/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// bring your own functionality
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);
	         // allocate memory for the impulse response
        h.resize(num_taps, 0.0);
        float norm_cutoff;
        norm_cutoff = Fc / (Fs / 2.0);
        for(int i = 0; i < num_taps; i++) {
                if (i == static_cast<int>((num_taps-1)/2)){
                        h[i] = norm_cutoff;
                }
                else{
                        float tempy = std::sin(PI*norm_cutoff*(i-((num_taps)/2)))/(PI*norm_cutoff*(i-((num_taps)/2)));
                        h[i] = norm_cutoff*tempy;
                }
                h[i] = h[i] * std::pow((std::sin((i*PI)/num_taps)), 2);
        }
}

//Impuse respone for the bandpass filter
void impulseResponseBPF(float Fb, float Fe, float Fs, int num_taps, std::vector<float> &h){
	//Initializations
	float norm_pass = (Fe-Fb)/(Fs/2);
	float n_half = (num_taps-1)/2;
	float norm_center = ((Fe+Fb)/2)/(Fs/2);

	h.resize(num_taps, 0.0);

	for(int i=0; i<num_taps; i++){
		if(i == n_half){
			h[i] = norm_pass;
		} else {
			
			//If curr != n_half, assign
			h[i] = norm_pass*sin(PI*(norm_pass/2)*(i-n_half))/(PI*(norm_pass/2)*(i-n_half));
		}
		h[i] = h[i] * cos(i*PI*norm_center);
        h[i] = h[i] * pow(sin((i*PI)/num_taps), 2);
	}
}

//RRC Impulse response
void impulseResponseRRC(const float &Fs, const int &num_taps, std::vector<float> &h)
{
	//Duration of each symbol
	float T_symbol = 1.0/2375.0;
	//Roll off factor
	float beta = 0.90;
	//The response vector
	h.resize(num_taps,0.0); 
	//Additional declarations
	float t;

	//Loop for the RRC
	for(int k = 0; k < num_taps; k++)
	{
		t = (float)(((float)(k-num_taps/2.0)/Fs));
		if(t == 0.0)
		{
			h[k] = 1.0 + beta*((4/PI)-1);

		}
		else if(t == -T_symbol/(4.0*beta) || t == T_symbol/(4.0*beta))
		{
			h[k] = (beta/sqrt(2.0))*(((1.0+2.0/PI)*(sin(PI/(4.0*beta))))+((1.0-2.0/PI)*(cos(PI/(4.0*beta)))));
		}
		else
		{
			h[k] = (sin(PI*t*(1.0-beta)/T_symbol)+4.0*beta*(t/T_symbol)*cos(PI*t*(1.0+beta)/T_symbol))/(PI*t*(1.0-(4.0*beta*t/T_symbol)*(4.0*beta*t/T_symbol))/T_symbol);
		}
	}

}
// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi)
{
	// bring your own functionality
	// allocate memory for the output (filtered) data
	y.resize(x.size()+h.size()-1, 0.0);
        //Set up an array to look at initial values. For single input its just zeros 
        //Count variable used to keep track of values 
        int count;
        for (auto n = 0; n < y.size(); n++){
                count = 0;
                for (auto k = 0; k < h.size(); k++)
		{
			if((n-k >= 0) && (n-k < x.size()))
			{
                                y[n] += x[n-k]*h[k];
			}
			else
			{
                                y[n] += zi[zi.size()-1-count]*h[k];
				count += 1;
			}
                }
        }
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = x[x.size()-zi.size()-1+i]; 
	}
}

//Downsamples and does convoloution at the same time in order to go speedy 
void convolveWithDecim(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num)
{
	//Creates vector for down sampled data
	y.resize(x.size()/decim_num);
	//Loops through vlaues to do convoloution
        int count;
        for (auto n = 0; n < y.size(); n++)
	{
                count = 0;
                for (auto k = 0; k < h.size(); k++)
		{
			if((decim_num*n-k >= 0) && (decim_num*n-k < x.size()))
			{
				//Multiplies n by the downsample number so not every value is used only the samples we need
                                y[n] += x[decim_num*n-k]*h[k];
			}
			//Previous state data
			else
			{
                                y[n] += zi[zi.size()-1-count]*h[k];
				count += 1;
			}
                }
        }
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = x[x.size()-zi.size()-1+i]; 
	}
}

//Downsamples and does convoloution at the same time in order to go speedy 
void convolveWithDecimPointer(std::vector<float> &y, float* &x, const unsigned int block_size, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num)
{
	//Creates vector for down sampled data
	y.resize(block_size/decim_num);
	//Loops through vlaues to do convoloution
        int count;
        for (auto n = 0; n < y.size(); n++)
	{
                count = 0;
                for (auto k = 0; k < h.size(); k++)
		{
			if((decim_num*n-k >= 0) && (decim_num*n-k < block_size))
			{
				//Multiplies n by the downsample number so not every value is used only the samples we need
                                y[n] += (float)(*(x+(decim_num*n-k)))*h[k];
			}
			//Previous state data
			else
			{
                                y[n] += zi[zi.size()-1-count]*h[k];
				count += 1;
			}
                }
        }
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = x[block_size-zi.size()-1+i]; 
	}
}
//Meant to combine the the convoloution of both I and Q samples. Hopefully faster but who knows in this economy 
void convolveWithDecimIQ(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi,std::vector<float> &y1, const std::vector<float> &x1, std::vector<float> &zi1, const int &decim_num)
{
//Creates vector for down sampled data
y.resize(x.size()/decim_num);
y1.resize(x1.size()/decim_num);
//Loops through vlaues to do convoloution
int count;
for (auto n = 0; n < y.size(); n++)
{
	count = 0;
	for (auto k = 0; k < h.size(); k++)
	{
		if((decim_num*n-k >= 0) && (decim_num*n-k < x.size()))
		{
			//Multiplies n by the downsample number so not every value is used only the samples we need
			y[n] += x[decim_num*n-k]*h[k];
			y1[n] += x1[decim_num*n-k]*h[k];
		}
		//Previous state data
		else
		{
			y[n] += zi[zi.size()-1-count]*h[k];
			y1[n] += zi1[zi1.size()-1-count]*h[k];
				count += 1;
			}
                }
        }
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = x[x.size()-zi.size()-1+i]; 
		zi1[i] = x1[x1.size()-zi.size()-1+i]; 
	}
}

//Convoloution in mode 1, takes upsamling and decimation as an argument 
void convolveWithDecimMode1(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num, const int &up_sample)
{
	//Creates vector for down sampled data
	y.resize(x.size()*up_sample/decim_num);
	//Loops through vlaues to do convoloution
        int count;
	//Iterator for the impulse response 
	int iter;
        for (auto n = 0; n < y.size(); n++)
	{
		//Reset the values
                count = 0;
		iter = 1;
                for (auto k = 0; k < h.size(); k= k+iter)
		{
			//Once you reach a non-zero start iterating 24 
			if(decim_num*n-k == 0 || (decim_num*n-k) %up_sample == 0) iter = up_sample; 
			else continue;
			
			//
			if((decim_num*n-k >= 0) && (decim_num*n-k < x.size()*up_sample))
			{
                                y[n] += x[(decim_num*n-k) /up_sample]*h[k];
			}
			//Previous state data
			else
			{
				y[n] += zi[(zi.size()-1-count)/up_sample]*h[k];
			}

			count += 1;
                }
        }
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = x[x.size()-zi.size()-1+i]; 
	}
}
//Pointer mode 1 
void convolveWithDecimMode1Pointer(std::vector<float> &y,float* &x, const unsigned int block_size, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num, const int &up_sample)
{
	//Creates vector for down sampled data
	y.resize(block_size*up_sample/decim_num);
	//Loops through vlaues to do convoloution
        int count;
	//Iterator for the impulse response 
	int iter;
        for (auto n = 0; n < y.size(); n++)
	{
		//Reset the values
                count = 0;
		iter = 1;
                for (auto k = 0; k < h.size(); k= k+iter)
		{
			//Once you reach a non-zero start iterating 24 
			if(decim_num*n-k == 0 || (decim_num*n-k) %up_sample == 0) iter = up_sample; 
			else continue;
			
			//
			if((decim_num*n-k >= 0) && (decim_num*n-k < block_size*up_sample))
			{
                                y[n] += (float)(*(x+((decim_num*n-k)/up_sample)))*h[k];
			}
			//Previous state data
			else
			{
				y[n] += zi[(zi.size()-1-count)/up_sample]*h[k];
			}

			count += 1;
                }
        }
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = x[block_size-zi.size()-1+i]; 
	}
}

void convolveWithDecimMode1RDS(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num, const int &up_sample)
{
	//Creates vector for down sampled data
	y.resize(x.size()*up_sample/decim_num);
	//Loops through vlaues to do convoloution
        int count;
	//Iterator for the impulse response 
	int iter;
        for (auto n = 0; n < y.size(); n++)
	{
		//Reset the values
                count = 0;
		iter = 1;
                for (auto k = 0; k < h.size(); k= k+iter)
		{
			//Once you reach a non-zero start iterating 24 
			if(decim_num*n-k == 0 || (decim_num*n-k) %up_sample == 0) iter = up_sample; 
			else continue;
			
			//
			if((decim_num*n-k >= 0) && (decim_num*n-k < x.size()*up_sample))
			{
                                y[n] += x[(decim_num*n-k) /up_sample]*h[k];
			}
			//Previous state data
			else
			{
				y[n] += zi[(zi.size()-1-count)/up_sample]*h[k];
			}

			count += 1;
                }
		y[n] = y[n]*up_sample;
        }
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = x[x.size()-zi.size()-1+i]; 
	}
}
//Downsamples and does convoloution at the same time in order to go speedy 
void convolveWithDecimSquare(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num)
{
	//Creates vector for down sampled data
	y.resize(x.size()/decim_num);
	//Loops through vlaues to do convoloution
        int count;
        for (auto n = 0; n < y.size(); n++)
	{
                count = 0;
                for (auto k = 0; k < h.size(); k++)
		{
			if((decim_num*n-k >= 0) && (decim_num*n-k < x.size()))
			{
				//Multiplies n by the downsample number so not every value is used only the samples we need
                                y[n] += pow(x[decim_num*n-k],2)*h[k];
			}
			//Previous state data
			else
			{
                                y[n] += zi[zi.size()-1-count]*h[k];
				count += 1;
			}
                }
        }
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = x[x.size()-zi.size()-1+i]; 
	}
}

//combines mixing and convoloutoo
void convolveWithDecimAndMixer(std::vector<float> &y, const std::vector<float> &x,const std::vector<float> &x1 , const std::vector<float> &h, std::vector<float> &zi, const int &decim_num)
{
	//Creates vector for down sampled data
	y.resize(x.size()/decim_num);
	//Loops through vlaues to do convoloution
        int count;
        for (auto n = 0; n < y.size(); n++)
	{
                count = 0;
                for (auto k = 0; k < h.size(); k++)
		{
			if((decim_num*n-k >= 0) && (decim_num*n-k < x.size()))
			{
				//Multiplies n by the downsample number so not every value is used only the samples we need
                                y[n] += x[decim_num*n-k]*x1[decim_num*n-k]*h[k]*2;
			}
			//Previous state data
			else
			{
                                y[n] += zi[zi.size()-1-count]*h[k];
				count += 1;
			}
                }
        }
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = x[x.size()-zi.size()-1+i]; 
	}
}
