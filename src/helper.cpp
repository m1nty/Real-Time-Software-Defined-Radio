/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
#include "dy4.h"
#include "iofunc.h"


//fmPLL without state
void fmPLL(std::vector<float> &ncoOut, std::vector<float> pllIn, float freq, float Fs, float ncoScale = 1.0, float phaseAdjust = 0.0, float normBandwidth = 0.01){
    float Cp = 2.666;
	float Ci = 3.555;

	// gain for the proportional term
	float Kp = (normBandwidth)*Cp;
	// gain for the integrator term
	float Ki = (normBandwidth*normBandwidth)*Ci;

	// output array for the NCO
	ncoOut.resize(pllIn.size()+1);

	// initialize internal state
	float integrator = 0.0;
	float phaseEst = 0.0;
	float feedbackI = 1.0;
	float feedbackQ = 0.0;
	ncoOut[0] = 1.0;
	float trigOffset = 0;
    float errorI, errorQ, errorD;
    float trigArg;
	// note: state saving will be needed for block processing

	for (int k=0; k<pllIn.size(); k++)
    {

		// phase detector
		errorI = pllIn[k] * (+feedbackI);  // complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ) ; // feedback complex exponential

		// four-quadrant arctangent discriminator for phase error detection
		errorD = atan2(errorQ, errorI);

		// loop filter
		integrator = integrator + Ki*errorD;

		// update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator;

		// internal oscillator
		trigArg = 2*PI*(freq/Fs)*(trigOffset+k+1) + phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
    }
}