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
//MEGA function 
void pllCombine(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &decim_num,std::vector<float> &ncoOut, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth , pll_state_type &pll_state)
{
	float Cp = 2.666;
	float Ci = 3.555;

	// Gain Initilizations
	float Ki = (normBandwidth*normBandwidth)*Ci;
	float Kp = (normBandwidth)*Cp;

	ncoOut.resize(x.size()+1);

	//PLL State Type Initializer
	float integrator = pll_state.integrator;
	float phaseEst = pll_state.phaseEst;
	float feedbackI = pll_state.feedbackI;
	float feedbackQ = pll_state.feedbackQ;
	ncoOut[0] = pll_state.ncoLast;
	float trigOffset = pll_state.trigOffset;
	
	//Creates vector for down sampled data
	y.resize(x.size()/decim_num);
	//Loops through vlaues to do convoloution
        int count;
        for (auto n = 0; n < y.size(); n++)
	{
                count = 0;
                for (auto k = 0; k < h.size(); k++)
		{
			if((decim_num*n-k >= 0) && (decim_num*n-k < x.size()))
			{
				//Multiplies n by the downsample number so not every value is used only the samples we need
                                y[n] += pow(x[decim_num*n-k],2)*h[k];
			}
			//Previous state data
			else
			{
                                y[n] += zi[zi.size()-1-count]*h[k];
				count += 1;
			}
                }
		float errorI = y[n] * (+feedbackI);
		float errorQ = y[n] * (-feedbackQ);
		float errorD = atan2(errorQ, errorI);  //Phase Error Detector
		integrator += Ki*errorD; //Must increment by integrator gain
		phaseEst += Kp*errorD + integrator; //Phase Estimate Update

		//Oscillator Config
		float trigArg = 2*PI*(freq/Fs)*(trigOffset+n+1)+phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[n+1] = cos(trigArg*ncoScale + phaseAdjust);
	}
	//Assigns next zi value 
	for(auto i = 0; i < zi.size(); i++){
		zi[i] = pow(x[x.size()-zi.size()-1+i],2); 
	}
    	//PLL State Type Edits  
    	pll_state.integrator = integrator;
   	pll_state.phaseEst = phaseEst;
    	pll_state.feedbackI = feedbackI;
    	pll_state.feedbackQ = feedbackQ;
    	pll_state.ncoLast= ncoOut[ncoOut.size()-1];
    	pll_state.trigOffset= trigOffset + y.size();
}
