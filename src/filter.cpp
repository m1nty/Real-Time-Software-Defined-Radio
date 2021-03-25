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

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi)
{
	// bring your own functionality
	// allocate memory for the output (filtered) data
	y.resize(x.size()+h.size()-1, 0.0);
        //Set up an array to look at inital values. For single input its just zeros 
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
