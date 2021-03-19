

/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Implementing fmMonoBasic.py
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "rf_module.h"

int main(int argc, char* argv[])
{
	//This is for mode 1 and 0
	//Need to improve when start working on mode 0 
	int mode = 0;
	if(argc < 2)
	{
		std::cerr << "Operating in mode 0" << std::endl;
	}
	else if (argc == 2)
	{
		mode = atoi(argv[1]);
		if(mode!=1)
		{
			std::cerr << "Wrong mode " <<mode <<std::endl;
			exit(1);
		}
	}
	else
	{
		std::cerr << "Usage " <<argv[0] <<std::endl;
	}

	//define values needed for processing
	//TODO May need to use a different data type
	int rf_Fs = 2400000;
	int rf_Fc = 100000;
	int rf_taps = 151;
	int rf_decim = 10;
	int audio_Fs = 48000;
	int audio_Fc = 16000;
	int audio_taps = 151; 
	int audio_decim = 5; 
	unsigned int block_size = 1024 * rf_decim * audio_decim * 2;
	int position = 0;
	unsigned int block_id = 0;
	std::vector<float> iq_data, i_data, q_data,iq_filter_coeff,i_inital, q_inital;
	std::vector<float> i_filter, q_filter, i_down, q_down;
	std::vector<float> prev_phase,demod_data; 
	std::vector<float> mono_coeff,audio_inital,audio_block, audio_filter;
	std::vector<short int> audio_data;
	std::vector<std::complex<float>> Xf;
	std::vector<float> vector_index;
	std::vector<float> Xmag;
	std::vector<float> psd_est, freq;
	std::vector<float> psd_est1, freq1;

	//Sets some inital values
	prev_phase.resize(2,0.0);
	audio_inital.resize(audio_taps-1,0.0);
	i_inital.resize(rf_taps-1,0.0);
	q_inital.resize(rf_taps-1,0.0);
	audio_inital.resize(audio_taps-1,0.0);
	
	iq_data.resize(block_size);
	i_data.resize(block_size/2);
	q_data.resize(block_size/2);
	while(true)
	{

		//Now that data is in, need to seperate the I and Q data into seperate vectors
		//Maybe this index is wrong
		readStdInBlock(block_size, block_id, iq_data);

		for(auto i = 0; i < (block_size)/2; i ++)
		{
			i_data[i] = iq_data[i*2];
			q_data[i] = iq_data[1+i*2];
		}

		//now generate the filter coefficents
		impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, iq_filter_coeff);

		//Filtering 
		convolveWithDecimIQ(i_filter, i_data, iq_filter_coeff, i_inital,q_filter, q_data, q_inital, rf_decim);
		//convolveWithDecim(q_filter, q_data, iq_filter_coeff, q_inital, rf_decim);

		//Demoadulate data
		fmDemodArctan(i_filter, q_filter, prev_phase, demod_data);
		
		//NOTE: To generate the PSD, it take a lil bit because of the Fourier transform. Uncoment this block to get PSD data
		//attempt to plot
		//if(position == 1024 * rf_decim * audio_decim * 2)
		//{
		//	std::vector<float> slice_data = std::vector<float>(demod_data.begin(), demod_data.begin() + 5*NFFT);
		//	DFT(slice_data, Xf);
		//	computeVectorMagnitude(Xf, Xmag);
		//	vector_index.clear();
		//	genIndexVector(vector_index, Xmag.size());
		//	logVector("demod_freq", vector_index, Xmag);
		//	int Fs1 = (2400000/10)/1000; //Divide by 1000 to make it easier to see
		//	estimatePSD(slice_data, NFFT, Fs1, freq, psd_est);
		//	logVector("demod_psd", freq, psd_est); 		
		//}
		
		//Now need use intermediate frequency to get audio.
		impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, mono_coeff);
		convolveWithDecim(audio_block, demod_data, mono_coeff, audio_inital, audio_decim);

		//NOTE: To generate the PSD, it take a lil bit because of the Fourier transform. Uncoment this block to get PSD data
		//if(position == 1024 * rf_decim * audio_decim * 2)
		//{
		//	std::vector<float> slice_data1 = std::vector<float>(audio_data.begin(), audio_data.begin() + 5*NFFT);
		//	int Fs2 = 48000/1000; //Divide by 1000 to make it easier to see
		//	estimatePSD(slice_data1, NFFT, Fs2, freq1, psd_est1);
		//	logVector("demod_audio", freq1, psd_est1); 
		//}

		//Write blocks to stdout
		audio_data.resize(audio_block.size());
		for(unsigned int l=0; l<audio_block.size(); l++)
		{
			//If block is null or a weird value, set it to zero
			if(std::isnan(audio_block[l])) 
			{
				audio_data[l] =0;
			}
			//Cast the value to short int for standard out and then scale it
			else 
			{
				audio_data[l] = static_cast<short int>(audio_block[l] *16384);
			}
		}
		
		//Write to standard out
		fwrite(&audio_data[0],sizeof(short int),audio_data.size(), stdout);

		//Fill the elements with zero to ensure when filtering happens no weird numbers are added
		std::fill(i_filter.begin(), i_filter.end(), 0);
		std::fill(q_filter.begin(), q_filter.end(), 0);
		std::fill(audio_block.begin(), audio_block.end(), 0);
	
		//Will keep iterating the block till there is nothing coming in from standard in 
		if((std::cin.rdstate()) != 0)
		{
			break;
		}
		else
		{
			block_id ++;
		}
	}
	std::cerr << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png"<< std::endl;

	return 0;
}
