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
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void blockProcessing(std::vector<float> &y1,std::vector<float> &y2,  const std::vector<float> &audio_right, const std::vector<float> &audio_left,unsigned short int num_taps, const float &Fc, const float &Fs, const int &block_size);

#endif // DY4_FILTER_H
