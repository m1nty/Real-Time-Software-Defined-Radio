/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
#include "dy4.h"
#include "iofunc.h"
#include "helper.h"

void fmPLL(std::vector<float> &ncoOut, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth , pll_state &state){
    
    float Cp = 2.666;
	float Ci = 3.555;

	// Gain Initilizations
    float Ki = (normBandwidth*normBandwidth)*Ci;
	float Kp = (normBandwidth)*Cp;

	ncoOut.resize(pllIn.size()+1);

    //PLL State Type Initializer
    float integrator = state.integrator;
    float phaseEst = state.phaseEst;
    float feedbackI = state.feedbackI;
    float feedbackQ = state.feedbackQ;
    ncoOut[0] = state.ncoLast;
    float trigOffset = state.trigOffset;

    // float phaseAdjust = 0.0;

	for (int k=0; k<pllIn.size(); k++)
    {
		float errorI = pllIn[k] * (+feedbackI);  // complex conjugate of the
		float errorQ = pllIn[k] * (-feedbackQ); // feedback complex exponential
    
    	// four-quadrant arctangent discriminator for phase error detection
		float errorD = atan2(errorQ, errorI);

		// loop filter
		integrator += Ki*errorD;

		// update phase estimate
		phaseEst += Kp*errorD + integrator;

		// internal oscillator
		float trigArg = 2*PI*(freq/Fs)*(trigOffset+k+1)+phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
    }

    //Update State Variables so they are saved to the struct object in main    
    state.integrator = integrator;
    state.phaseEst = phaseEst;
    state.feedbackI = feedbackI;
    state.feedbackQ = feedbackQ;
    state.ncoLast= ncoOut[ncoOut.size()-1];
    state.trigOffset= trigOffset + pllIn.size();

    //Resize to return 1:end of array
    ncoOut = std::vector<float>(ncoOut.begin(), ncoOut.end()-1);

    //block_data = std::vector<float>(block_data.begin()+1,block_data.end());

}