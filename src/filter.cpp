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
