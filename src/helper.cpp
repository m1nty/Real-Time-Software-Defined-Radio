/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/
#include "dy4.h"
#include "iofunc.h"
#include "helper.h"

void fmPLL(std::vector<float> &ncoOut, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth , pll_state_type &pll_state){
    
    float Cp = 2.666;
	float Ci = 3.555;

	// Gain Initilizations
    float Ki = (normBandwidth*normBandwidth)*Ci;
	float Kp = (normBandwidth)*Cp;

	ncoOut.resize(pllIn.size()+1);

    //PLL State Type Initializer
    float integrator = pll_state.integrator;
    float phaseEst = pll_state.phaseEst;
    float feedbackI = pll_state.feedbackI;
    float feedbackQ = pll_state.feedbackQ;
    ncoOut[0] = pll_state.ncoLast;
    float trigOffset = pll_state.trigOffset;

	for (int k=0; k<pllIn.size(); k++)
    {
		float errorI = pllIn[k] * (+feedbackI);
		float errorQ = pllIn[k] * (-feedbackQ);
    	float errorD = atan2(errorQ, errorI);  //Phase Error Detector
		integrator += Ki*errorD; //Must increment by integrator gain
		phaseEst += Kp*errorD + integrator; //Phase Estimate Update

		//Oscillator Config
		float trigArg = 2*PI*(freq/Fs)*(trigOffset+k+1)+phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
    }

    //PLL State Type Edits  
    pll_state.integrator = integrator;
    pll_state.phaseEst = phaseEst;
    pll_state.feedbackI = feedbackI;
    pll_state.feedbackQ = feedbackQ;
    pll_state.ncoLast= ncoOut[ncoOut.size()-1];
    pll_state.trigOffset= trigOffset + pllIn.size();

    //ncOut Resize
    ncoOut = std::vector<float>(ncoOut.begin(), ncoOut.end()-1);
}
void fmPLLIQ(std::vector<float> &ncoOut,std::vector<float> &ncoOutQ, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth , pll_state_type &pll_state){
    
    float Cp = 2.666;
	float Ci = 3.555;

	// Gain Initilizations
    float Ki = (normBandwidth*normBandwidth)*Ci;
	float Kp = (normBandwidth)*Cp;

	ncoOut.resize(pllIn.size()+1);
	ncoOutQ.resize(pllIn.size()+1);

    //PLL State Type Initializer
    float integrator = pll_state.integrator;
    float phaseEst = pll_state.phaseEst;
    float feedbackI = pll_state.feedbackI;
    float feedbackQ = pll_state.feedbackQ;
    ncoOut[0] = pll_state.ncoLast;
    ncoOutQ[0] = pll_state.ncoLast;
    float trigOffset = pll_state.trigOffset;

	for (int k=0; k<pllIn.size(); k++)
    {
		float errorI = pllIn[k] * (+feedbackI);
		float errorQ = pllIn[k] * (-feedbackQ);
    	float errorD = atan2(errorQ, errorI);  //Phase Error Detector
		integrator += Ki*errorD; //Must increment by integrator gain
		phaseEst += Kp*errorD + integrator; //Phase Estimate Update

		//Oscillator Config
		float trigArg = 2*PI*(freq/Fs)*(trigOffset+k+1)+phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
		ncoOutQ[k+1] = sin(trigArg*ncoScale + phaseAdjust);
    }

    //PLL State Type Edits  
    pll_state.integrator = integrator;
    pll_state.phaseEst = phaseEst;
    pll_state.feedbackI = feedbackI;
    pll_state.feedbackQ = feedbackQ;
    pll_state.ncoLast= ncoOut[ncoOut.size()-1];
    pll_state.trigOffset= trigOffset + pllIn.size();

    //ncOut Resize
    ncoOut = std::vector<float>(ncoOut.begin(), ncoOut.end()-1);
}
