/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
	for (auto m = 0; m < Xf.size(); m++) {
		for (auto k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
	// only the positive frequencies
	Xmag.resize(Xf.size(), static_cast<float>(0));
  	for (auto i = 0; i < Xf.size(); i++) {
    		Xmag[i] = std::abs(Xf[i])/Xf.size();
  	}	
}

// add your own code to estimate the PSD
void estimatePSD(const std::vector<float> &samples,const int &NFFT_1, const int &Fs, std::vector<float> &freq, std::vector<float> &psd_est)
{
	//Create some variables needed to estimate the PSD
	int freq_bin; 
	freq_bin = NFFT_1;
	float df;
	//Defines each frequency segment
	df = static_cast<float>(Fs)/static_cast<float>(freq_bin); 
	//std::cout << "The df value is " << df << "\n";
	//Resizes and set up the frequency vector
	freq.resize((Fs/2)/df);
	float number = 0.0;
	//Need to create that frequency vector 
	for(auto i =0; i < freq.size(); i++)
	{
		freq[i] = number;
		number += df; 
	}
	//Creates Hann window
	std::vector<float> hann;
	hann.resize(freq_bin);
	for(auto j = 0; j < hann.size(); j++)
	{
		hann[j] = std::pow(std::sin(j*PI/static_cast<float>(freq_bin)),2.0);
	}
	//Creates a empty vector for PSD samples
	std::vector<float> psd_list; 
	//Create number of segments
	int no_segments = std::floor(samples.size()/freq_bin);
	//Create vector for look
	std::vector<float> windowed_samples;
	//Define vectors needed for the loop
	std::vector<std::complex<float>> Xf,Xf_pos;
	Xf_pos.resize(freq_bin/2);
	std::vector<float> psd_seg;
	psd_seg.resize(Xf_pos.size());
	int position = 0; 	
	psd_list.resize(no_segments*(freq_bin/2));
	windowed_samples.resize(freq_bin); 
	//Now need to create a loop that iterates through all of these segments
	for(auto k=0; k<no_segments; k++)
	{
		//Multiplies by hann window 
		for(auto t=0; t < freq_bin; t++)
		{
			//windowed_samples[count] = samples[t]*hann[count];
			windowed_samples[t] = samples[position+t]*hann[t];
		}
		//Discrete Fourier Transfer
		DFT(windowed_samples, Xf);
		//Now, we want to keep only positive frequencies
		//We are going to try and combin the 2 with setting up psd segments, see how it goes
		for(auto l = 0; l < Xf_pos.size(); l++) 
		{
			//Lowkey may not need this. Will look later to see if can optimize
			Xf_pos[l] = Xf[l]; 
			psd_seg[l] = 2*((1/(static_cast<float>(Fs*freq_bin)/2.0))*std::pow(std::abs(Xf_pos[l]),2)); 
			//Find mag in dB
			psd_seg[l] = 10*(std::log10(psd_seg[l])); 
			psd_list[position+l] = psd_seg[l];
		}
		//Iterate the position that looks through the sample list
		position += freq_bin/2; 
	}
	//Now that we have all of the psd values, need to set them up to all the freq  
	psd_est.resize(freq_bin/2);
	for(auto n = 0; n < freq_bin/2; n++)
	{
		for(auto m = 0; m < no_segments; m++)
		{
			psd_est[n] += psd_list[n + m*(freq_bin/2)];
		}
		psd_est[n] = psd_est[n]/static_cast<float>(no_segments); 
	}
}
