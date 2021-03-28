
/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_RF_MODULE_H
#define DY4_RF_MODULE_H

// add headers as needed
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

// declaration of a function prototypes
void fmDemodArctan(const std::vector<float> &I, const std::vector<float> &Q,std::vector<float> &prev_phase, std::vector<float> &fm_demod);
#endif // DY4_FILTER_H
