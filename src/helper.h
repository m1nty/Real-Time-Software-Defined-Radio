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

struct pll_state{
	float integrator, phaseEst, feedbackI, feedbackQ, trigOffset,ncoLast;
};

void fmPLL(std::vector<float> &ncoOut, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, pll_state &state);

#endif