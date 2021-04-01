/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_SUPPORT_H
#define DY4_SUPPORT_H

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

struct pll_state_type{
	float integrator, phaseEst, feedbackI, feedbackQ, trigOffset,ncoLast;
};

void fmPLL(std::vector<float> &ncoOut, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, pll_state_type &pll_state);

void fmPLLIQ(std::vector<float> &ncoOut,std::vector<float> &ncoOutQ, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth , pll_state_type &pll_state);

void pllCombine(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num,std::vector<float> &ncoOut, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth , pll_state_type &pll_state);
#endif
