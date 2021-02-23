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

//Place for band pass filter 


//Function for block processing
void blockProcessing(std::vector<float> &y1,std::vector<float> &y2,  const std::vector<float> &audio_right, const std::vector<float> &audio_left,unsigned short int num_taps, const float &Fc, const float &Fs, const int &block_size)
{
	//Defining needed variables
	std::vector<float> zi_right, zi_left,h_right,h_left, x_temp1, x_temp2, y_temp1, y_temp2;
	int position = 0; 
	//Making sure require vectors for convoloution are set up 
	zi_right.resize(num_taps-1, 0.0);
	zi_left.resize(num_taps-1, 0.0);
	x_temp1.resize(block_size, 0.0);
	x_temp2.resize(block_size, 0.0);
	//Performs impulse reponse to get h values
	impulseResponseLPF(Fs, Fc, num_taps, h_right);
	impulseResponseLPF(Fs, Fc, num_taps, h_left);
	y1.resize(audio_right.size()+h_right.size()-1, 0.0);
        y2.resize(audio_left.size()+h_left.size()-1, 0.0);
	//Keep going through each block and perform convoloution
	while(true)
	{
		//First need for loop to assign values to x	
		for(auto i = 0; i < block_size; i++) 
		{
			x_temp1[i] = audio_right[position+i];
			x_temp2[i] = audio_left[position+i];
		}
		//Does the convoloution
		//printRealVector(zi_right);
		convolveFIR(y_temp1,x_temp1,h_right,zi_right);	
		convolveFIR(y_temp2,x_temp2,h_left,zi_left);	
		//Gets the y values 
		for(auto j = 0; j < y_temp1.size(); j++)
		{
			y1[position+j] = y_temp1[j];
			y2[position+j] = y_temp2[j];
		}
		//Resets all the values with zeros
		std::fill(y_temp1.begin(), y_temp1.end(), 0.0);
		std::fill(y_temp2.begin(), y_temp2.end(), 0.0);
		//Iterates the positon  
		position = position + block_size; 
		if (position > audio_right.size())
		{
			break;
		}
	}
}
