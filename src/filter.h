/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h);
void impulseResponseBPF(float Fb, float Fe, float Fs, int num_taps, std::vector<float> &h);

void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi);

void convolveWithDecim(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num);

void convolveWithDecimIQ(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi,std::vector<float> &y1, const std::vector<float> &x1, std::vector<float> &zi1, const int &decim_num);

void convolveWithDecimMode1(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num, const int &up_sample);

void convolveWithDecimPointer(std::vector<float> &y, float* &x, const unsigned int block_size, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num); 

#endif // DY4_FILTER_H
