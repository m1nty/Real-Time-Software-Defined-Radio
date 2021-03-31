
/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Implementing fmSupportLib.py refactored in C++
*/


#include "dy4.h"
#include "genfunc.h"

void fmDemodArctan(const std::vector<float> &I, const std::vector<float> &Q,std::vector<float> &prev_phase, float* &queue_block) {
//Function to perform demodulation on the IQ samples 
	prev_phase = {0.0, 0.0};
	float prev_I = prev_phase[0];
	float prev_Q = prev_phase[1];

	for(unsigned int k = 0 ; k<I.size() ; k++)
	{
		if(pow(I[k],2) + pow(Q[k],2) == 0)
		{
			*(queue_block+k) = 0.0;
		} 
		else
		{
			*(queue_block+k) =(I[k]*(Q[k]-prev_Q)-Q[k]*(I[k]-prev_I))/(pow(I[k],2)+pow(Q[k],2));
		}
		prev_I = I[k];
        	prev_Q = Q[k];
		prev_phase[0] = prev_I;
		prev_phase[1] = prev_Q;
	}
}
