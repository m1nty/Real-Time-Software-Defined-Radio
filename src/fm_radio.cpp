/*
Comp Eng 3DY4 (Computer Systems Integration Project)
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "rf_module.h"
#include "queue"
#include "mutex"
#include "condition_variable"
#include <thread>
#include "helper.h"

#include <iostream> 

#define QUEUE_BLOCKS 5 

//NOTE FOR WHOEVER TESTING
//To run in mode 0 type: cat ../data/my_samples_953.raw | ./experiment | aplay -c 1 -f S16_LE -r 48000
//To run in mode 1 type: cat ../data/my_samples_953.raw | ./experiment 1 | aplay -c 1 -f S16_LE -r 48000

//Rf thread
//TODO Eventuall add these functions to a seperate file so its less ugly
//Should be ready for mode 1. All that needs to be done is assign the sampling frequency 
void rf_thread(int &mode, std::queue<std::vector<float>> &sync_queue, std::mutex &radio_mutex, std::condition_variable &cvar) 
{
	//Need to set up conditional for this 
	int rf_Fs;
	//Depending on the mode sets the sampling frequency to the corresponding value
	if(mode == 1) rf_Fs = 2500000;
	else rf_Fs = 2400000;
	
	//These variables are the same no matter what
	int rf_Fc = 100000;
	int rf_taps = 151;
	int rf_decim = 10;
	int audio_decim = 5; 
	unsigned int block_id = 0;
	unsigned int block_size = 1024 * rf_decim * audio_decim * 2;
	//define nessisary vectors
	std::vector<float> iq_data, i_data, q_data,iq_filter_coeff,i_inital, q_inital;
	std::vector<float> i_filter, q_filter, i_down, q_down;
	std::vector<float> prev_phase,demod_data; 

	//Resize inital states used in convoloution
	prev_phase.resize(2,0.0);
	i_inital.resize(rf_taps-1,0.0);
	q_inital.resize(rf_taps-1,0.0);
	
	//Resizes blocks
	iq_data.resize(block_size);
	i_data.resize(block_size/2);
	q_data.resize(block_size/2);
	std::cerr << "rf_Fs = " << rf_Fs << std::endl;
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

		//If queue is full should wait till it is empty
		std::unique_lock<std::mutex> queue_lock(radio_mutex);
		if(sync_queue.size() == QUEUE_BLOCKS)
		{
			cvar.wait(queue_lock);
		}
		//push vector onto queue 
		sync_queue.push(demod_data);

		//Fills with zeros
		std::fill(i_filter.begin(), i_filter.end(), 0);
		std::fill(q_filter.begin(), q_filter.end(), 0);
		
		//iterate block id
		block_id ++;
		
		//Unlock and notify 
		queue_lock.unlock();
		cvar.notify_one();

		//If nothing at STD in it breaks
		if((std::cin.rdstate()) != 0)
		{
			break;
		}
	}
}

//Thread for mono and stero 
//Need to add all the shit for mode 0  
void mono_stero_thread(int &mode, std::queue<std::vector<float>> &sync_queue, std::mutex &radio_mutex, std::condition_variable &cvar) 
{
	//Depending on the mode sets the sampling frequency to the corresponding value
	int rf_decim = 10;
	int audio_Fs = 240000;
	int audio_Fc = 16000;
	int audio_taps = 151; 
	int audio_decim = 5; 
	int stereo_decim = 5; 
	unsigned int block_size = 1024*rf_decim*audio_decim*2;

	unsigned int block_id = 0;
	int audio_up = 1;

	pll_state statePLL;
	statePLL.integrator = 0.0;
	statePLL.phaseEst = 0.0;
	statePLL.feedbackI = 1.0;
	statePLL.feedbackQ = 0.0;
	statePLL.trigOffset = 0.0;
	statePLL.ncoLast = 1.0;
	
	//If mode 1, change som values, and define the up sampler value 
	if (mode == 1) 
	{
		audio_Fs = 6000000;
		audio_decim = 125; 
		audio_up = 24; 
		audio_taps = audio_taps*audio_up; 
	}

	//Sets up nessisary vectors 
	std::vector<float> mono_coeff,audio_inital,audio_block, audio_filter;
	std::vector<float> stereo_coeff,recovery_initial,stereo_initial,extraction_initial, recovery_coeff, extraction_coeff;
	std::vector<float> stereo_lr, stereo_filt, stereo_data, bpf_recovery, recovery_pll, bpf_extraction, mixed;

	std::vector<short int> audio_data;
	
	//Sets some inital values
	audio_inital.resize(audio_taps-1,0.0);

	recovery_initial.resize(audio_taps-1,0.0);
	stereo_initial.resize(audio_taps-1,0.0);
	extraction_initial.resize(audio_taps-1,0.0);

	std::vector<float> block99;

	mixed.resize(5120,0.0);
	// std::cerr << "audio_Fs = " << audio_Fs << std::endl;

	//Creates the filter coefficents and then
	//Need to scale this when with mono mode 1
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, mono_coeff);
	impulseResponseBPF(18.5e3, 19.5e3, audio_Fs, audio_taps, recovery_coeff);
	impulseResponseBPF(22e3, 54e3, audio_Fs, audio_taps, extraction_coeff);
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, stereo_coeff);

	//randome var
	int mult = 1;
	while (true) 
	{
		//Creates lock
		std::unique_lock<std::mutex> queue_lock(radio_mutex);
		//Waits until there is something in the queue
		if(sync_queue.empty())
		{
			cvar.wait(queue_lock);
		}

		//Mode 1, with all the upsamping and pull shit 
		if(mode == 1)
		{
			std::vector<float> demod_data = sync_queue.front();
			convolveWithDecimMode1(audio_block, demod_data, mono_coeff, audio_inital, audio_decim, audio_up);
			mult = 24;


			//-----------------------STEREO CARRIER RECOVERY-------------------------------
			convolveWithDecim(bpf_recovery, demod_data, recovery_coeff, recovery_initial, 1);
        	fmPLL(recovery_pll, bpf_recovery, 19e3, 240e3,2.0,0.0, 0.01,statePLL);

			//-----------------------STEREO CHANNEL EXTRACTION-------------------------------
			convolveWithDecim(bpf_extraction, demod_data, extraction_coeff, extraction_initial, 1);

			//-----------------------STEREO PROCESSING-------------------------------
			//Mixing

			for (int i = 0; i <bpf_extraction.size();i++){
				mixed[i] = bpf_extraction[i]*recovery_pll[i];
			}		
			//LPF
			convolveWithDecim(stereo_filt, mixed, stereo_coeff, stereo_initial,stereo_decim);


			//Combiner
			stereo_lr.resize(2*audio_block.size());
			for(int i = 0 ; i<audio_block.size(); i++){
				stereo_lr[2*i] = (audio_block[i] + stereo_filt[i])/2; // left channel
				stereo_lr[2*i+1] = (audio_block[i] - stereo_filt[i])/2; // right channel
			
			}
		}
		//if in mode 0 
		else
		{
			std::vector<float> demod_data = sync_queue.front();
			convolveWithDecim(audio_block, demod_data, mono_coeff, audio_inital, audio_decim);

			//-----------------------STEREO CARRIER RECOVERY-------------------------------
			convolveWithDecim(bpf_recovery, demod_data, recovery_coeff, recovery_initial, 1);
        	fmPLL(recovery_pll, bpf_recovery, 19e3, 240e3,2.0,0.0, 0.01,statePLL);

			//-----------------------STEREO CHANNEL EXTRACTION-------------------------------
			convolveWithDecim(bpf_extraction, demod_data, extraction_coeff, extraction_initial, 1);

			//-----------------------STEREO PROCESSING-------------------------------
			//Mixing

			for (int i = 0; i <bpf_extraction.size();i++){
				mixed[i] = bpf_extraction[i]*recovery_pll[i];
			}		
			//LPF
			convolveWithDecim(stereo_filt, mixed, stereo_coeff, stereo_initial,stereo_decim);


			//Combiner
			stereo_lr.resize(2*audio_block.size());
			for(int i = 0 ; i<audio_block.size(); i++){
				stereo_lr[2*i] = (audio_block[i] + stereo_filt[i])/2; // left channel
				stereo_lr[2*i+1] = (audio_block[i] - stereo_filt[i])/2; // right channel
			
			}
		}
		
		//Pops the value of the queue 
		sync_queue.pop();

		queue_lock.unlock();
		cvar.notify_one();
		
		//Write blocks to stdout
		audio_data.resize(stereo_lr.size());
		for(unsigned int l=0; l<stereo_lr.size(); l++)
		{
			//If block is null or a weird value, set it to zero
			if(std::isnan(stereo_lr[l])) 
			{
				audio_data[l] =0;
			}
			//Cast the value to short int for standard out and then scale it
			else 
			{
				audio_data[l] = static_cast<short int>(stereo_lr[l] *16384);
			}
		}
		
		//Write to standard out
		fwrite(&audio_data[0],sizeof(short int),audio_data.size(), stdout);

		//Fill the elements with zero to ensure when filtering happens no weird numbers are added
		std::fill(audio_block.begin(), audio_block.end(), 0);
		std::fill(stereo_filt.begin(), stereo_filt.end(), 0);
		mixed.clear();

		//iterate block id
		block_id ++;

		//Will keep iterating the block till there is nothing coming in from standard in 
		if((std::cin.rdstate()) != 0)
		{
			break;
		}
	}
}
	
//Need to get args working
int main(int argc, char* argv[])
{
	//This is for mode 1 and 0
	//Need to improve when start working on mode 0 
	int mode = 0;
	//If no input, mode is 0
	if(argc < 2)
	{
		std::cerr << "Operating in mode 0" << std::endl;
	}
	//If 1 is the input, argc becomes 2
	else if (argc == 2)
	{
		mode = atoi(argv[1]);
		//If something sus happens exits 
		if(mode!=1)
		{
			std::cerr << "Wrong mode " <<mode <<std::endl;
			exit(1);
		}
		//if the argument is 1, the mode is set to 1
		else
		{
			std::cerr << "Operating in mode " << mode <<  std::endl;
		}
	}
	//Incase somesithing really suspect happens
	else
	{
		std::cerr << "Usage " <<argv[0] <<std::endl;
		exit(1);
	}

	//Define stuff for threads
	std::queue<std::vector<float>> sync_queue;
	std::mutex radio_mutex;
	std::condition_variable cvar; 

	//THIS IS STUFF FOR PLOTTING IN CASE WE NEEDA DO THAT 
	//define values needed for processing
	std::vector<std::complex<float>> Xf;
	std::vector<float> vector_index;
	std::vector<float> Xmag;
	std::vector<float> psd_est, freq;
	std::vector<float> psd_est1, freq1;
	
	//Creates threads
	std::thread rf = std::thread(rf_thread, std::ref(mode), std::ref(sync_queue), std::ref(radio_mutex), std::ref(cvar));
	std::thread mono_stero = std::thread(mono_stero_thread, std::ref(mode), std::ref(sync_queue), std::ref(radio_mutex), std::ref(cvar));

	//Once all standard in blocks have been read, threads are joined
	rf.join();
	mono_stero.join(); 

	//Print this case I always forget the command
	std::cerr << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png"<< std::endl;

	return 0;
}