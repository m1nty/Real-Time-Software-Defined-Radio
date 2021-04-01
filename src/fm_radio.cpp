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
#include <math.h>
#include <chrono>

#include <iostream> 

#define QUEUE_BLOCKS 5 
#define BLOCK_SIZE 307200 

//Rf thread where IQ samples read in from rf dongle. This is where they are filtered and demodulated 
void rf_thread(int &mode, std::queue<void *> &sync_queue,std::queue<void *> &rds_queue, std::mutex &radio_mutex,std::mutex &rds_mutex, std::condition_variable &cvar,std::condition_variable &cvar1) 
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
	unsigned int block_id = 0;
	unsigned int block_size = BLOCK_SIZE;
	//define nessisary vectors
	std::vector<float> iq_data, i_data, q_data,iq_filter_coeff,i_initial, q_initial;
	std::vector<float> i_filter, q_filter, i_down, q_down;
	std::vector<float> prev_phase; 

	//Stuff for pointer implementation
	static float queue_block[QUEUE_BLOCKS][BLOCK_SIZE];
	//Resize initial states used in convoloution
	prev_phase.resize(2,0.0);
	i_initial.resize(rf_taps-1,0.0);
	q_initial.resize(rf_taps-1,0.0);
	
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
		convolveWithDecimIQ(i_filter, i_data, iq_filter_coeff, i_initial,q_filter, q_data, q_initial, rf_decim);
		//convolveWithDecim(q_filter, q_data, iq_filter_coeff, q_initial, rf_decim);

		unsigned int queue_entry = block_id % QUEUE_BLOCKS;
		//Demoadulate data
		float * pointer_block = &queue_block[queue_entry][0];
		fmDemodArctan(i_filter, q_filter, prev_phase, pointer_block);
		
		std::unique_lock<std::mutex> queue_lock(radio_mutex);
		//Probs not right but issa attempt
		if(mode == 1)
		{
			if(sync_queue.size() == QUEUE_BLOCKS-1)
			{
				cvar.wait(queue_lock);
			}
			sync_queue.push((void *)&queue_block[queue_entry][0]);

			//Fills with zeros
			std::fill(i_filter.begin(), i_filter.end(), 0);
			std::fill(q_filter.begin(), q_filter.end(), 0);
			
			//iterate block id
			block_id ++;
			
			//Unlock and notify 
			queue_lock.unlock();
			cvar.notify_one();
			if((std::cin.rdstate()) != 0)
			{
				break;
			}
		}
		else
		{
			std::unique_lock<std::mutex> queue1_lock(rds_mutex);
			if(sync_queue.size() == QUEUE_BLOCKS-1||rds_queue.size() == QUEUE_BLOCKS-1) 
			{
				//std::cerr << "Issue so need to lock" << std::endl;
				if(sync_queue.size() == QUEUE_BLOCKS-1)
					cvar.wait(queue_lock);
				if(rds_queue.size() == QUEUE_BLOCKS-1)
					cvar1.wait(queue1_lock);
			}
			//push pointers to address onto queue 
			sync_queue.push((void *)&queue_block[queue_entry][0]);
			rds_queue.push((void *)&queue_block[queue_entry][0]);

			//Fills with zeros
			std::fill(i_filter.begin(), i_filter.end(), 0);
			std::fill(q_filter.begin(), q_filter.end(), 0);
			
			//iterate block id
			block_id ++;
			
			//Unlock and notify 
			queue_lock.unlock();
			cvar.notify_one();
			cvar1.notify_one();
			if((std::cin.rdstate()) != 0)
			{
				break;
			}
		}
	}
}

//Thread for mono and stero. This is where all of the audio data is processed  
void mono_stero_thread(int &mode, std::queue<void *> &sync_queue, std::mutex &radio_mutex, std::condition_variable &cvar) 
{
	//Depending on the mode sets the sampling frequency to the corresponding value
	int audio_Fs = 240000;
	int audio_Fc = 16000;
	int audio_taps = 151; 
	int audio_taps_1 = 151; 
	int audio_decim = 5; 
	int stereo_decim = 5; 
	unsigned int block_size = BLOCK_SIZE;

	unsigned int block_id = 0;
	int audio_up = 1;

	//For Pll block processing
	pll_state_type pll_state;
	pll_state.integrator = 0.0;
	pll_state.phaseEst = 0.0;
	pll_state.feedbackI = 1.0;
	pll_state.feedbackQ = 0.0;
	pll_state.trigOffset = 0.0;
	pll_state.ncoLast = 1.0;
	
	//If mode 1, change som values, and define the up sampler value 
	if (mode == 1) 
	{
		audio_Fs = 6000000;
		audio_decim = 125; 
		audio_up = 24; 
		audio_taps = audio_taps*audio_up; 
	}

	//Sets up vectors 
	std::vector<float> mono_coeff,audio_initial,audio_block, audio_filter;
	std::vector<float> stereo_coeff,recovery_initial,stereo_initial,extraction_initial, recovery_coeff, extraction_coeff;
	std::vector<float> stereo_lr, stereo_filt, stereo_data, bpf_recovery, recovery_pll, bpf_extraction, mixed;
	std::vector<short int> audio_data;
	
	//Sets some initial values
	audio_initial.resize(audio_taps-1,0.0);

	recovery_initial.resize(audio_taps-1,0.0);
	stereo_initial.resize(audio_taps-1,0.0);
	extraction_initial.resize(audio_taps-1,0.0);

	mixed.resize(block_size/20,0.0);
	// std::cerr << "audio_Fs = " << audio_Fs << std::endl;

	//Creates the filter coefficents and then
	//Need to scale this when with mono mode 1
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, mono_coeff);
	impulseResponseBPF(18.5e3, 19.5e3, audio_Fs, audio_taps_1, recovery_coeff);
	impulseResponseBPF(22e3, 54e3, audio_Fs, audio_taps_1, extraction_coeff);
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, stereo_coeff);

	//for parcebels theorm
	int mult = 1;
	
	//Loop where processes occur
	while (true) 
	{
		//Creates lock
		std::unique_lock<std::mutex> queue_lock(radio_mutex);
		//Waits until there is something in the queue
		if(sync_queue.empty())
		{
			cvar.wait(queue_lock);
		}
		//Assigns front of quene to vector
		float *ptr_block = (float *)sync_queue.front();
		//Pops the value of the queue 
		sync_queue.pop();
		queue_lock.unlock();
		cvar.notify_one();

		//Mode 1, with Upsample/Pulls
		if(mode == 1)
		{
			convolveWithDecimMode1Pointer(audio_block, ptr_block, block_size/20, mono_coeff, audio_initial, audio_decim, audio_up);
			mult = 24;

			//-----------------------STEREO CARRIER RECOVERY-------------------------------
			convolveWithDecimPointer(bpf_recovery, ptr_block,BLOCK_SIZE/20, recovery_coeff, recovery_initial, 1);
			fmPLL(recovery_pll, bpf_recovery, 19e3, 240e3,2.0,0.0, 0.01,pll_state);

			//-----------------------STEREO CHANNEL EXTRACTION-------------------------------
			convolveWithDecimPointer(bpf_extraction,ptr_block,BLOCK_SIZE/20 , extraction_coeff, extraction_initial, 1);

			//-----------------------STEREO PROCESSING-------------------------------
			//Mixing
			for (int i = 0; i < bpf_extraction.size();i++)
			{
				mixed[i] = bpf_extraction[i]*recovery_pll[i];
			}		
			//LPF
			convolveWithDecimMode1(stereo_filt, mixed, stereo_coeff, stereo_initial,stereo_decim,audio_up);
			//Combiner
			stereo_lr.resize(2*audio_block.size());
			for(int i = 0 ; i < audio_block.size(); i++)
			{
				stereo_lr[2*i] = (audio_block[i] + stereo_filt[i])/2; // left channel
				stereo_lr[2*i+1] = (audio_block[i] - stereo_filt[i])/2; // right channel
			}
		}
		//if in mode 0 
		else
		{
			//Gets mono Data
			convolveWithDecimPointer(audio_block, ptr_block,BLOCK_SIZE/20, mono_coeff, audio_initial, audio_decim);

			//-----------------------STEREO CARRIER RECOVERY-------------------------------
			convolveWithDecimPointer(bpf_recovery, ptr_block,BLOCK_SIZE/20, recovery_coeff, recovery_initial, 1);
			fmPLL(recovery_pll, bpf_recovery, 19e3, 240e3,2.0,0.0, 0.01,pll_state);

			//-----------------------STEREO CHANNEL EXTRACTION-------------------------------
			convolveWithDecimPointer(bpf_extraction,ptr_block,BLOCK_SIZE/20 , extraction_coeff, extraction_initial, 1);

			//-----------------------STEREO PROCESSING-------------------------------
			//Mixing
			for (int i = 0; i <bpf_extraction.size();i++)
			{
				mixed[i] = bpf_extraction[i]*recovery_pll[i];
			}		
			//LPF
			convolveWithDecim(stereo_filt, mixed, stereo_coeff, stereo_initial,stereo_decim);

			//Combiner
			stereo_lr.resize(2*audio_block.size());
			for(int i = 0 ; i<audio_block.size(); i++)
			{
				stereo_lr[2*i] = (audio_block[i] + stereo_filt[i])/2; // left channel
				stereo_lr[2*i+1] = (audio_block[i] - stereo_filt[i])/2; // right channel
			}
		}
		
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
				audio_data[l] = static_cast<short int>(stereo_lr[l] *16384*mult);
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
		if((std::cin.rdstate()) != 0&&sync_queue.empty())
		{
			break;
		}
	}
}

//Inital RDS Thread which does all of the filtering
void rds_thread(int &mode, std::queue<void *> &rds_queue, std::mutex &rds_mutex, std::condition_variable &cvar, std::queue<std::vector<float>> &frame_queue, std::mutex &frame_mutex, std::condition_variable &cvarframe) 
{
	//Will only perform these process if in mode 0 
	if(mode == 0)
	{
		//Defining the vectors from python
		std::vector<float> extract_RDS_coeff, pre_state_extract, square_coeff, square_state, lpf_coeff_rds, lpf_3k_state, anti_img_coeff, anti_img_state, rrc_coeff, rrc_state;
		//Other vectors
		std::vector<float> extract_rds,extract_rds_squared,pre_Pll_rds, post_Pll,mixed,lpf_filt_rds, resample_rds, rrc_rds;
		//Defining constants
		int num_taps = 151; 
		unsigned int block_id = 0;
		//For first BPF
		float initial_RDS_lower_freq = 54000;
		float initial_RDS_higher_freq = 60000;
		float Fs = 240000;
		//For second BPF
		float squared_lower_freq = 113500;
		float squared_higher_freq = 114500;
		//PLL
		float freq_centered = 114000;
		float phase_adj = PI/3.3-PI/1.5;
		pll_state_type pll_state;
		pll_state.integrator = 0.0;
		pll_state.phaseEst = 0.0;
		pll_state.feedbackI = 1.0;
		pll_state.feedbackQ = 0.0;
		pll_state.trigOffset = 0.0;
		pll_state.ncoLast = 1.0;
		//LPF 
		float cutoff_LPF = 3000; 
		//Rational Resampler
		int upsample_val = 19;
		int downsample_val = 80;
		float cutoff_anti_img = 57000/2;
		//Values for clock recovery
		std::vector<float> symbols_I;

		//Define initial states
		pre_state_extract.resize(num_taps-1); 
		square_state.resize(num_taps-1);
		lpf_3k_state.resize(num_taps-1);
		anti_img_state.resize(num_taps*19-1); 
		rrc_state.resize(num_taps-1);
		//Get coeeficents
		impulseResponseBPF(initial_RDS_lower_freq, initial_RDS_higher_freq, Fs, num_taps,extract_RDS_coeff);
		impulseResponseBPF(squared_lower_freq, squared_higher_freq, Fs, num_taps, square_coeff);
		impulseResponseLPF(Fs, cutoff_LPF, num_taps, lpf_coeff_rds);
		impulseResponseLPF(Fs*19, cutoff_anti_img, num_taps*19, anti_img_coeff);
		impulseResponseRRC(57000, num_taps, rrc_coeff);
	
		//Loop to calculate all of it 
		while(true)
		{
			//Creates lock
			std::unique_lock<std::mutex> queue_lock(rds_mutex);
			//Waits until there is something in the queue
			if(rds_queue.empty())
			{
				cvar.wait(queue_lock);
			}
			//Assigns front of quene to vector
			float *ptr_block = (float *)rds_queue.front();
			//Pops the value of the queue 
			rds_queue.pop();
			queue_lock.unlock();
			cvar.notify_one();

			// ****************************************************************** 
			// -----------------------RDS Data Processing------------------------ 
			// ******************************************************************
				
			// ------------------------Extraction--------------------------------
			// Performs convoloution to extract the data 
			convolveWithDecimPointer(extract_rds, ptr_block,BLOCK_SIZE/20, extract_RDS_coeff, pre_state_extract, 1.0);

			// ---------------------Carrier Recovery-----------------------------
			//Combined the squaring, BPF and Pll for a bit of a time save 
			//Phase_adj was determined by looking at the constalation diagrams. But all of the Quadrature calculation functions were removed to increase run time
			pllCombine(pre_Pll_rds, extract_rds, square_coeff,square_state, 1,post_Pll, freq_centered, 240000, 0.5,phase_adj-PI/1.4 , 0.001 , pll_state);
			
			// ---------------------Demodulation-mixed----------------------------
			//Low pass filter combined with mixer
			convolveWithDecimAndMixer(lpf_filt_rds, post_Pll, extract_rds, lpf_coeff_rds, lpf_3k_state, 1);

			//Resampler
			//Resamples and passes through the anti-imaging filter at the same time
			convolveWithDecimMode1RDS(resample_rds, lpf_filt_rds, anti_img_coeff,anti_img_state,downsample_val,upsample_val);
			
			//RRC Filter
			convolveWithDecim(rrc_rds, resample_rds, rrc_coeff, rrc_state, 1);

			//Push to thread which dedicated to the rest of data processing in RDS
			std::unique_lock<std::mutex> frame_lock(frame_mutex);
			if(frame_queue.size() == QUEUE_BLOCKS-1)
			{
				cvarframe.wait(frame_lock);
			}
			//Push onto queue 
			frame_queue.push(rrc_rds);
			//std::cerr << "Size of diff_bits = " << diff_bits.size() <<std::endl; 
			frame_lock.unlock();
			cvarframe.notify_one();
			//Fills these vectors with zeros
			std::fill(extract_rds.begin(), extract_rds.end(), 0);
			std::fill(pre_Pll_rds.begin(), pre_Pll_rds.end(), 0);
			//std::fill(post_Pll.begin(), post_Pll.end(), 0);
			std::fill(lpf_filt_rds.begin(), lpf_filt_rds.end(), 0);
			std::fill(resample_rds.begin(), resample_rds.end(), 0);
			std::fill(rrc_rds.begin(), rrc_rds.end(), 0);
			//iterate block id
			block_id ++;

			//Will keep iterating the block till there is nothing coming in from standard in 
			if((std::cin.rdstate()) != 0&&rds_queue.empty())
			{
				break;
			}
		}
	}
}
	
//Thread which does all of the clock recovery, manchester decoding and frame sync
void frame_thread(int &mode, std::queue<std::vector<float>> &frame_queue,std::queue<void *> &rds_queue, std::mutex &frame_mutex, std::condition_variable &cvarframe)
{
	if(mode == 0)
	{
		//Block id
		int block_id = 0;
		//Stuff for clock sync 
		unsigned int initial_offset = 0; 
		//variables for manchester decoding
		int count_0_pos =0;
		int count_1_pos = 0;
		unsigned int start_pos = 0;
		//Some block processing values 
		float lonely_bit = 0;
		int front_bit = 0; 
		//Values for differential decoding
		int prebit = 0; 
		unsigned int offset = 0;
		//Vecotrs needed for these processes 
		std::vector<float> rrc_rds;
		std::vector<int> bit_stream; 
		std::vector<int> diff_bits; 
		std::vector<float> symbols_I;
		std::vector<int> potential_syndrome,block,prev_sync_bits;
		potential_syndrome.resize(10);
		block.resize(26);
		prev_sync_bits.resize(27);
		//Poistion variables for frame sync
		unsigned int position = 0;
		unsigned int printposition = 0; 
		int last_position = -1;
		//The parity matrix
		std::vector<std::vector<int>> H {{1,0,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,1},{1,0,1,1,0,1,1,1,0,0},{0,1,0,1,1,0,1,1,1,0},{0,0,1,0,1,1,0,1,1,1},{1,0,1,0,0,0,0,1,1,1},{1,1,1,0,0,1,1,1,1,1},{1,1,0,0,0,1,0,0,1,1},{1,1,0,1,0,1,0,1,0,1},{1,1,0,1,1,1,0,1,1,0},{0,1,1,0,1,1,1,0,1,1},{1,0,0,0,0,0,0,0,0,1},{1,1,1,1,0,1,1,1,0,0}, {0,1,1,1,1,0,1,1,1,0},{0,0,1,1,1,1,0,1,1,1},{1,0,1,0,1,0,0,1,1,1},{1,1,1,0,0,0,1,1,1,1}, {1,1,0,0,0,1,1,0,1,1}};
		//Defining the syndromes
		std::vector<int> syndrome_A {1,1,1,1,0,1,1,0,0,0};
		std::vector<int> syndrome_B {1,1,1,1,0,1,0,1,0,0};
		std::vector<int> syndrome_C {1,0,0,1,0,1,1,1,0,0};
		std::vector<int> syndrome_D {1,0,0,1,0,1,1,0,0,0};
		//Loop in which all of the processes take place 
		while(true) 
		{
			//Check for reading in from the queue 
			std::unique_lock<std::mutex> frame_lock(frame_mutex);
			//Waits until there is something in the queue
			if(frame_queue.empty())
			{
				cvarframe.wait(frame_lock);
			}
			//Assigns front of quene to vector
			rrc_rds = frame_queue.front();
			//std::cerr << "diff_bits size " << diff_bits.size()<<std::endl;
			//Pops the value of the queue 
			frame_queue.pop();
			frame_lock.unlock();
			cvarframe.notify_one();

			//If its the first block find the inital offset of the first 24 bits to find the peak then take every 24th bit 
			if(block_id == 0)
			{
				//finds the index of the max
				//float placehold = abs(rrc_rds[0]);
				float placehold = abs(rrc_rds[0]);
				for(unsigned int i = 1; i < 24;i++)
				{
					if(std::abs(rrc_rds[i]) > placehold)
					{
						placehold = std::abs(rrc_rds[i]); 
						initial_offset = i; 
					}
				}
				std::cerr << "initial offset for clock recovery = " << initial_offset <<std::endl; 
			}
			//symbols_I.resize(std::floor((rrc_rds.size())/24));
			symbols_I.resize((int)((rrc_rds.size())/24));
			
			//Gets the symbols from the rrc array. Gets every 24th sample from the offset
			for(unsigned int k=0; k < symbols_I.size();k++)
			{
				symbols_I[k] = rrc_rds[24*k+initial_offset];
				//std::cerr << symbols_I[k] <<std::endl;
			}

			//Then figures out what the index will need to be for the next block
			//NOTE: For some reason, in python this is needed. But we get more syndromes when this is commented out in C++
			//It makes the constalations look alot nicer in C++ but for some reason we loss all the syndromes in block 1
			//for(unsigned int j = 0; j < 24; j++)
			//{	
			//	//std::cerr << rrc_rds[rrc_rds.size()-24+j] << std::endl;
			//	if(rrc_rds[rrc_rds.size()-24+j] == symbols_I[symbols_I.size()-1])
			//	{
			//		initial_offset = 24-j; 
			//	}
			//}

			// ---------------------RDS Data Processing----------------------------
			//initial screening to determine if we start at bit 0 or 1 when doing manchester decoding
			if(block_id == 0)
			{
				//Loops through small chunck of symbols vector 
				for(unsigned int j; j < symbols_I.size()/4; j++)
				{
					//Checks for doubles. Which ever start point has the least doubles is where will start from 
					if((symbols_I[2*j] > 0 && symbols_I[2*j+1] > 0) || (symbols_I[2*j] < 0 && symbols_I[2*j+1]<0))
						count_0_pos += 1;
					else if((symbols_I[2*j+1] > 0 && symbols_I[2*j+2] > 0) || (symbols_I[2*j+1] < 0 && symbols_I[2*j+2]<0))
						count_1_pos += 1;
				}
				//Depending on what start position has more doubles, determine the appropriate start position 
				if(count_0_pos > count_1_pos)
					start_pos = 1;
				else if(count_1_pos > count_0_pos) 
					start_pos = 0;
			}
			//Now the recovery 
			bit_stream.resize((int)(symbols_I.size()/2)-start_pos,0);
			
			//Uses prev bit just in case
			//In the case that the start position is 1, we need to use the last unused bit from the prev block and the first bit of the current block,
			//in order to get the bit that will be appended to the front 
			if(start_pos == 1 && block_id != 0)
			{
				//Which ever symbol is bigger will indicate if we go from a high to a low or vice versa 
				if(lonely_bit > symbols_I[0])
					front_bit = 1 ;
				else if(symbols_I[0] > lonely_bit)
					front_bit = 0;
			}
			//Figures out what every bit is
			for(unsigned int k =0; k<bit_stream.size(); k++) 
			{
				//Break if we exceed the vector bounds
				if(start_pos+2*k+1 > symbols_I.size()-1)
					break;
				//Which ever symbol is bigger will indicate if we go from a high to a low or vice versa 
				//Then we can set the value of bit stream to the corresponding bit
				if(symbols_I[2*k+start_pos] > symbols_I[2*k+1+start_pos])
					bit_stream[k] = 1;
				else if(symbols_I[2*k+start_pos] < symbols_I[2*k+1+start_pos])
					bit_stream[k] = 0;
			}
			//If start pos 1 append bit to front and set value for the bit of the end 
			if(start_pos == 1) 
			{
				bit_stream.insert(bit_stream.begin(), front_bit);
				//Set last unused symbol to the value of lonely bit 
				lonely_bit = symbols_I[symbols_I.size()-1]; 
			}

			//Differential decoding
			//Set initial pre bit and offset
			if(block_id == 0)
			{
				//Initall we need to use the first bit of the bit stream for differential decoding
				prebit = bit_stream[0];
				offset = 1;
			}
			else
			{
				//Afterwards just use the last bit from the previous bit	
				offset = 0; 
			}
			//Perform XOR on bits to perform differential decoding
			diff_bits.resize(bit_stream.size()-offset,0);
			for(unsigned int t = 0; t< diff_bits.size(); t++)
			{
				diff_bits[t] = prebit ^ bit_stream[t+offset];
				//Set current bit to the bit that will be used in the next loop iteration
				prebit = bit_stream[t+offset];
			}
			//Just in case set last bit to the prebit for the next block 
			prebit = bit_stream[bit_stream.size()-1];

			//Prints block Number
			std::cerr << " " << std::endl;
			std::cerr << "****************Prcoessing Block: " << block_id<< "****************" << std::endl;

			//Frame Sync
			//insert remaining bits to front of vector
			//Insert bits not used to the beginning of this vecotr so no syndromes are lost 
			if(block_id != 0){
				diff_bits.insert(diff_bits.begin(), prev_sync_bits.begin(), prev_sync_bits.end());
			}
			
			//Sets inital postion
			position = 0;
			while(true){
				//Creates a block of size, and get a new block each cycle shifted by position
				for(unsigned int y = 0; y < 26; y++)
				{
					block[y] = diff_bits[y+position];
				}
				//Matrix multiplication in order to get a vector with 10 elements to check if its a syndrome
				std::fill(potential_syndrome.begin(), potential_syndrome.end(), 0);
				for(unsigned int i = 0 ; i<potential_syndrome.size() ; i++)
				{
					for(unsigned int j = 0 ; j<26 ; j++)
					{
						//XOR instead of add, AND instad of multiply
						potential_syndrome[i] = potential_syndrome[i] ^ (block[j] && H[j][i]);
					}
				}
				//The following conditionals see if the caclulated value corresponds to a syndrome
				//Checks if syndrome A
				if(potential_syndrome == syndrome_A){ 
					if(last_position == -1 || printposition-last_position == 26){ 
						last_position = printposition;
						std::cerr << "Syndrome A at position " << printposition << std::endl;
						last_position = printposition;
					}
					else{
						std::cerr << "False positive Syndrome A at position " << printposition << std::endl;
					}
				}
				//Checks if syndrome B
				else if(potential_syndrome == syndrome_B){ 
					if(last_position == -1 || printposition-last_position == 26){ 
						std::cerr << "Syndrome B at position " << printposition << std::endl;
						last_position = printposition;
					}
					else{
						std::cerr << "False positive Syndrome B at position " << printposition << std::endl;
					}
				}
				//Checks if syndrome C
				else if(potential_syndrome == syndrome_C){ 
					if(last_position == -1 || printposition-last_position == 26){ 
						std::cerr << "Syndrome C at position " << printposition << std::endl;
						last_position = printposition;
					}
					else{
						std::cerr << "False positive Syndrome C at position " << printposition << std::endl;
					}
				}
				//Checks if syndrome D
				else if(potential_syndrome == syndrome_D){
					if(last_position == -1 || printposition-last_position == 26){ 
						std::cerr << "Syndrome D at position " << printposition << std::endl;
						last_position = printposition;
					}
					else{
						std::cerr << "False positive Syndrome D at position " << printposition << std::endl;
					}
				}
				//Breaks once it reaches the end
				//Iterate the posiition variabel
				position += 1;
				if(position+26 > diff_bits.size()-1)
				{
					//std::cerr<< "Size of diff_bit = " << diff_bits.size() << std::endl; 
					//std::cerr <<"Position left at " << position << std::endl;
					break;
				}
				printposition+=1;
			}
			//Creates vector  of bits not used to be used in the next block 
			for(unsigned int g = 0; g < prev_sync_bits.size(); g++)
			{
				prev_sync_bits[g] = diff_bits[position-1+g];
			}

			//Iterate the block ID
			block_id ++;
			//Runs until both queues are empty and there is nothing at the standard in 
			if((std::cin.rdstate()) != 0 && frame_queue.empty()&&rds_queue.empty())
			{
				break;
			}
		}
	}
}

//Main function which also takes arguments to determine the mode in which the radio will operate in
int main(int argc, char* argv[])
{
	//This is for mode 1 and 0
	//Need to improve when start working on mode 0 
	int mode = 0;
	//If no input, mode is 0
	std::cerr << argc <<std::endl;
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
	std::queue<void *> sync_queue;
	std::queue<void *> rds_queue;
	std::mutex radio_mutex,rds_mutex;
	std::condition_variable cvar,cvar1,cvarframe; 

	//THIS IS STUFF FOR PLOTTING IN CASE WE NEEDA DO THAT 
	//define values needed for processing
	std::vector<std::complex<float>> Xf;
	std::vector<float> vector_index;
	std::vector<float> Xmag;
	std::vector<float> psd_est, freq;
	std::vector<float> psd_est1, freq1;
	std::queue<std::vector<float>> frame_queue;
	std::mutex frame_mutex;
	
	//Creates threads
	std::thread rf = std::thread(rf_thread, std::ref(mode), std::ref(sync_queue), std::ref(rds_queue), std::ref(radio_mutex),std::ref(rds_mutex), std::ref(cvar),std::ref(cvar1));
	std::thread mono_stero = std::thread(mono_stero_thread, std::ref(mode), std::ref(sync_queue), std::ref(radio_mutex), std::ref(cvar));
	std::thread rds  = std::thread(rds_thread, std::ref(mode), std::ref(rds_queue), std::ref(rds_mutex), std::ref(cvar1),std::ref(frame_queue), std::ref(frame_mutex), std::ref(cvarframe));
	std::thread frame= std::thread(frame_thread, std::ref(mode), std::ref(frame_queue),std::ref(rds_queue), std::ref(frame_mutex), std::ref(cvarframe));

	//Once all standard in blocks have been read, threads are joined
	rf.join();
	mono_stero.join(); 
	rds.join();
	frame.join();

	//Print this case I always forget the command
	std::cerr << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png"<< std::endl;

	return 0;
}
