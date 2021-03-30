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

#include <iostream> 

#define QUEUE_BLOCKS 5 
#define BLOCK_SIZE 307200 

//NOTE FOR WHOEVER TESTING
//To run in mode 0 type: cat ../data/my_samples_953.raw | ./experiment | aplay -c 1 -f S16_LE -r 48000
//To run in mode 1 type: cat ../data/my_samples_953.raw | ./experiment 1 | aplay -c 1 -f S16_LE -r 48000

//Rf thread
//TODO Eventuall add these functions to a seperate file so its less ugly
void rf_thread(int &mode, std::queue<void *> &sync_queue,std::queue<void *> &rds_queue, std::mutex &radio_mutex, std::condition_variable &cvar,std::condition_variable &cvar1) 
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
	unsigned int block_size = BLOCK_SIZE;
	//define nessisary vectors
	std::vector<float> iq_data, i_data, q_data,iq_filter_coeff,i_inital, q_inital;
	std::vector<float> i_filter, q_filter, i_down, q_down;
	std::vector<float> prev_phase,demod_data; 

	//Stuff for pointer implementation
	static float queue_block[QUEUE_BLOCKS][BLOCK_SIZE];
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

		unsigned int queue_entry = block_id % QUEUE_BLOCKS;
		//Demoadulate data
		float * pointer_block = &queue_block[queue_entry][0];
		fmDemodArctan(i_filter, q_filter, prev_phase, demod_data, pointer_block);
		
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
			if(sync_queue.size() == QUEUE_BLOCKS-1||rds_queue.size() == QUEUE_BLOCKS-1) //|| sync_queue.size() != rds_queue.size())
			{
				//std::cerr << "Issue so need to lock" << std::endl;
				if(sync_queue.size() == QUEUE_BLOCKS-1)
					cvar.wait(queue_lock);
				if(rds_queue.size() == QUEUE_BLOCKS-1)
					cvar1.wait(queue_lock);
			}
			//if(sync_queue.size() == QUEUE_BLOCKS-1)
			//{
			//	cvar.wait(queue_lock);
			//}
			//if(rds_queue.size() == QUEUE_BLOCKS-1)
			//{
			//	cvar1.wait(queue_lock);
			//}
			//push vector onto queue 
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

//Thread for mono and stero 
//Need to add all the shit for mode 0  
void mono_stero_thread(int &mode, std::queue<void *> &sync_queue, std::mutex &radio_mutex, std::condition_variable &cvar) 
{
	//Depending on the mode sets the sampling frequency to the corresponding value
	int rf_decim = 10;
	int audio_Fs = 240000;
	int audio_Fc = 16000;
	int audio_taps = 151; 
	int audio_taps_1 = 151; 
	int audio_decim = 5; 
	int stereo_decim = 5; 
	unsigned int block_size = BLOCK_SIZE;

	unsigned int block_id = 0;
	int audio_up = 1;

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

	//Sets up nessisary vectors 
	std::vector<float> mono_coeff,audio_inital,audio_block, audio_filter;
	std::vector<float> stereo_coeff,recovery_initial,stereo_initial,extraction_initial, recovery_coeff, extraction_coeff;
	std::vector<float> stereo_lr, stereo_filt, stereo_data, bpf_recovery, recovery_pll, bpf_extraction, mixed,demod_data;

	std::vector<short int> audio_data;
	
	//Sets some inital values
	audio_inital.resize(audio_taps-1,0.0);

	recovery_initial.resize(audio_taps-1,0.0);
	stereo_initial.resize(audio_taps-1,0.0);
	extraction_initial.resize(audio_taps-1,0.0);
	demod_data.resize(block_size/20);

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

		//Mode 1, with all the upsamping and pull shit 
		if(mode == 1)
		{
			convolveWithDecimMode1Pointer(audio_block, ptr_block, block_size/20, mono_coeff, audio_inital, audio_decim, audio_up);
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
			convolveWithDecimPointer(audio_block, ptr_block,BLOCK_SIZE/20, mono_coeff, audio_inital, audio_decim);

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

//RDS Thread
void rds_thread(int &mode, std::queue<void *> &rds_queue, std::mutex &radio_mutex, std::condition_variable &cvar) 
{
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
		float inital_RDS_lower_freq = 54000;
		float inital_RDS_higher_freq = 60000;
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
		//Values for rational resampler
		float rrc_Fs = 57000; 
		unsigned int rrc_taps = 151;
		//Values for clock recovery
		unsigned int inital_offset = 0; 
		int final_symb = 0; 
		//Differential decoding 
		int count_0_pos =0;
		int count_1_pos = 0;
		unsigned int start_pos = 0;
		float lonely_bit = 0;
		int front_bit = 0; 
		int prebit = 0; 
		unsigned int offset = 0;
		//Frame sync 
		unsigned int position = 0;
		unsigned int printposition = 0; 
		int last_position = -1;
		std::vector<int> potential_syndrome,block,prev_sync_bits;
		potential_syndrome.resize(10);
		block.resize(26);
		prev_sync_bits.resize(27);
		std::vector<std::vector<int>> H {{1,0,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,1,0,0,0,0,0},{0,0,0,0,0,1,0,0,0,0},{0,0,0,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,1},{1,0,1,1,0,1,1,1,0,0},{0,1,0,1,1,0,1,1,1,0},{0,0,1,0,1,1,0,1,1,1},{1,0,1,0,0,0,0,1,1,1},{1,1,1,0,0,1,1,1,1,1},{1,1,0,0,0,1,0,0,1,1},{1,1,0,1,0,1,0,1,0,1},{1,1,0,1,1,1,0,1,1,0},{0,1,1,0,1,1,1,0,1,1},{1,0,0,0,0,0,0,0,0,1},{1,1,1,1,0,1,1,1,0,0}, {0,1,1,1,1,0,1,1,1,0},{0,0,1,1,1,1,0,1,1,1},{1,0,1,0,1,0,0,1,1,1},{1,1,1,0,0,0,1,1,1,1}, {1,1,0,0,0,1,1,0,1,1}};
		std::vector<int> syndrome_A {1,1,1,1,0,1,1,0,0,0};
		std::vector<int> syndrome_B {1,1,1,1,0,1,0,1,0,0};
		std::vector<int> syndrome_C {1,0,0,1,0,1,1,1,0,0};
		std::vector<int> syndrome_D {1,0,0,1,0,1,1,0,0,0};

		//Define inital states
		pre_state_extract.resize(num_taps-1); 
		square_state.resize(num_taps-1);
		lpf_3k_state.resize(num_taps-1);
		anti_img_state.resize(num_taps*19-1); 
		rrc_state.resize(num_taps-1);
		//Get coeeficents
		impulseResponseBPF(inital_RDS_lower_freq, inital_RDS_higher_freq, Fs, num_taps,extract_RDS_coeff);
		impulseResponseBPF(squared_lower_freq, squared_higher_freq, Fs, num_taps, square_coeff);
		impulseResponseLPF(Fs, cutoff_LPF, num_taps, lpf_coeff_rds);
		impulseResponseLPF(Fs*19, cutoff_anti_img, num_taps*19, anti_img_coeff);
		impulseResponseRRC(57000, num_taps, rrc_coeff);
	
		//Loop to calculate all of it 
		while(true)
		{
			//Creates lock
			std::unique_lock<std::mutex> queue_lock(radio_mutex);
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
			std::cerr << " " << std::endl;
			std::cerr << "****************Prcoessing Block: " << block_id<< "****************" << std::endl;
				
			// ------------------------Extraction--------------------------------
			// Performs convoloution to extract the data 
			convolveWithDecimPointer(extract_rds, ptr_block,BLOCK_SIZE/20, extract_RDS_coeff, pre_state_extract, 1.0);

			// ---------------------Carrier Recovery-----------------------------
			//Squares the elements
			extract_rds_squared.resize(extract_rds.size());
			for(unsigned n = 0; n < extract_rds.size(); n++) 
			{
				extract_rds_squared[n] = extract_rds[n]*extract_rds[n]; 
			}
			//Second BPF
			convolveWithDecim(pre_Pll_rds, extract_rds_squared, square_coeff, square_state, 1);
			//Pll
			fmPLL(post_Pll, pre_Pll_rds, freq_centered, 240e3,0.5,phase_adj, 0.001,pll_state);

			// ---------------------Demodulation-mixed----------------------------
			//mixing 
			mixed.resize(post_Pll.size());
			for(unsigned m = 0; m < post_Pll.size(); m++)
			{
				mixed[m] = post_Pll[m] * extract_rds[m]*2;
				//std::cerr << mixed[m] << std::endl;
			}
			//LPF
			convolveWithDecim(lpf_filt_rds, mixed, lpf_coeff_rds, lpf_3k_state, 1);

			//Need to use the concoloution from mode 1 for upsampling
			convolveWithDecimMode1(resample_rds, lpf_filt_rds, anti_img_coeff,anti_img_state,downsample_val,upsample_val);
			for(unsigned int x = 0; x < resample_rds.size(); x++)
			{
				resample_rds[x] = resample_rds[x] * upsample_val;
				//std::cerr << resample_rds[x] << std::endl;
			}
			
			//Lets try some sus shit

			//RRC Filter
			convolveWithDecim(rrc_rds, resample_rds, rrc_coeff, rrc_state, 1);
			if(block_id == 3)
			{
				std::vector<float> time;
				genIndexVector(time, rrc_rds.size());
				std::cerr << "VectorLog time"<<std::endl;
				logVector("rrc", time, rrc_rds);
			}

			//Clock and data recovery
			//Determines where to initally sample 
			if(block_id == 0)
			{
				if(std::abs(rrc_rds[0]) < std::abs(rrc_rds[12]) || std::abs(rrc_rds[24]) < std::abs(rrc_rds[12]))
					inital_offset = 12;
				else
					inital_offset =0 ;
				std::cerr << "inital offset for clock recovery = " << inital_offset <<std::endl; 
			}
			std::vector<float> symbols_I;
			symbols_I.resize(std::floor((rrc_rds.size()-inital_offset)/24));
			//Gets the symbols from the rrc array 
			for(unsigned int k=0; k < symbols_I.size();k++)
			{
				symbols_I[k] = rrc_rds[24*k+inital_offset];
				//std::cerr << symbols_I[k] <<std::endl;
			}

			//NOTE May need to add something for block procssing here 
			
			// ---------------------RDS Data Processing----------------------------
			//Inital screening 
			if(block_id == 0)
			{
				//Loops through small chunck of symbols vector 
				for(unsigned int j; j < symbols_I.size()/4; j++)
				{
					//Checks for doubles 
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
				std::cerr << "Doubles when start pos is 0="<< count_0_pos << " Doubles when start pos is 1="<<count_1_pos << std::endl;
			}
			//Now the recovery 
			std::vector<int> bit_stream; 
			bit_stream.resize((int)((symbols_I.size()/2)-start_pos));
			
			//Uses prev bit just in case
			if(start_pos == 1 && block_id != 0)
			{
				if(lonely_bit > symbols_I[0])
					front_bit = 1 ;
				else if(symbols_I[0] > lonely_bit)
					front_bit = 0;
			}
			//Figures out what every bit is
			for(unsigned int k =0; k<bit_stream.size(); k++) 
			{
				if(start_pos+2*k+1 > symbols_I.size()-1)
					break;
				if(symbols_I[2*k+start_pos] > symbols_I[2*k+1+start_pos])
					bit_stream[k] = 1;
				else if(symbols_I[2*k+start_pos] < symbols_I[2*k+1+start_pos])
					bit_stream[k] = 0;
			}
			//If start pos 1 append bit to front and set value for the bit of the end 
			if(start_pos == 1) 
			{
				bit_stream.insert(bit_stream.begin(), front_bit);
				lonely_bit = symbols_I[-1]; 
			}
			//std::cerr << "Size of bit_stream = " << bit_stream.size() <<std::endl; 

			//Differential decoding
			//Set inital pre bit and offset
			if(block_id == 0)
			{
				prebit = bit_stream[0];
				offset = 1;
			}
			else
			{
				offset = 0; 
			}
			//Perform XOR on bits
			std::vector<int> diff_bits;
			diff_bits.resize(bit_stream.size()-offset);
			for(unsigned int t = 0; t< diff_bits.size(); t++)
			{
				diff_bits[t] = prebit ^ bit_stream[t+offset];
				prebit = bit_stream[t+offset];
			}
			prebit = bit_stream[-1];
			//std::cerr << "Size of diff_bits = " << diff_bits.size() <<std::endl; 

			//From Sync
			//insert remaining bits to front of vector
			if(block_id != 0){
				diff_bits.insert(diff_bits.begin(), prev_sync_bits.begin(), prev_sync_bits.end());
			}

			//std::cerr << "Size of diff_bits = " << diff_bits.size() <<std::endl; 

			position = 0;
			while(true){
				//Creates a block of size 26
				for(unsigned int y = position; y < position+26; y++)
				{
					block[y-position] = diff_bits[y];
				}
				//Matric multiplication 
				for(unsigned int i = 0 ; i<potential_syndrome.size() ; i++)
				{
					for(unsigned int j = 0 ; j<26 ; j++)
					{
						potential_syndrome[i] = (potential_syndrome[i] && !(block[j] && H[j][i])) || (!potential_syndrome[i] && (block[j] && H[j][i]));
					}
				}
				//convert to int
				//potential_syndrome = static_cast<short int>(potential_syndrome);
				//std::cerr << "Potential syndrome " << potential_syndrome[0] <<potential_syndrome[1]<< potential_syndrome[2]<<potential_syndrome[3]<<potential_syndrome[4]<<potential_syndrome[5]<<potential_syndrome[6]<<potential_syndrome[7]<<potential_syndrome[8]<<potential_syndrome[9]<<std::endl;
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
				position += 1;
				if(position+26 > diff_bits.size()-1){
					break;
				}
				printposition+=1;
			}
			//Creates list of bits not used 
			for(unsigned int g = 0; g < prev_sync_bits.size(); g++)
				prev_sync_bits[g] = diff_bits[position-1+g];

			
			//END
			std::fill(extract_rds.begin(), extract_rds.end(), 0);
			std::fill(pre_Pll_rds.begin(), pre_Pll_rds.end(), 0);
			std::fill(post_Pll.begin(), post_Pll.end(), 0);
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
	std::queue<void *> sync_queue;
	std::queue<void *> rds_queue;
	std::mutex radio_mutex;
	std::condition_variable cvar,cvar1; 

	//THIS IS STUFF FOR PLOTTING IN CASE WE NEEDA DO THAT 
	//define values needed for processing
	std::vector<std::complex<float>> Xf;
	std::vector<float> vector_index;
	std::vector<float> Xmag;
	std::vector<float> psd_est, freq;
	std::vector<float> psd_est1, freq1;
	
	//Creates threads
	std::thread rf = std::thread(rf_thread, std::ref(mode), std::ref(sync_queue), std::ref(rds_queue), std::ref(radio_mutex), std::ref(cvar),std::ref(cvar1));
	std::thread mono_stero = std::thread(mono_stero_thread, std::ref(mode), std::ref(sync_queue), std::ref(radio_mutex), std::ref(cvar));
	std::thread rds  = std::thread(rds_thread, std::ref(mode), std::ref(rds_queue), std::ref(radio_mutex), std::ref(cvar1));

	//Once all standard in blocks have been read, threads are joined
	std::cerr << "do we get here" << std::endl;
	rf.join();
	mono_stero.join(); 
	rds.join();

	//Print this case I always forget the command
	std::cerr << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png"<< std::endl;

	return 0;
}
