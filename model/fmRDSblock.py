
#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use "custom" fmDemodArctan
from fmSupportLib import fmDemodArctan
from fmSupportLib import my_filterImpulseResponse
from fmSupportLib import my_convoloution 
from fmPll import fmPll
from fmRRC import impulseResponseRootRaisedCosine

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 240000
audio_decim = 5
# add other settings for audio, like filter taps, ...
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5

#Parameters for RDS
#For the first BPF
inital_RDS_lower_freq = 54000
inital_RDS_higher_freq = 60000
#Second BPF
squared_lower_freq = 113500
squared_higher_freq = 114500
#PLL Parameters 
freq_centered = 114000
phase_shift = math.pi/4.5
#LPF
cutoff_LPF = 3000

#Parity Matrix
H = np.matrix([[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,1,1,0,1,1,1,0,0],[0,1,0,1,1,0,1,1,1,0],[0,0,1,0,1,1,0,1,1,1],[1,0,1,0,0,0,0,1,1,1],[1,1,1,0,0,1,1,1,1,1],[1,1,0,0,0,1,0,0,1,1],[1,1,0,1,0,1,0,1,0,1],[1,1,0,1,1,1,0,1,1,0],[0,1,1,0,1,1,1,0,1,1],[1,0,0,0,0,0,0,0,0,1],[1,1,1,1,0,1,1,1,0,0], [0,1,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,0,1,1,1],[1,0,1,0,1,0,0,1,1,1],[1,1,1,0,0,0,1,1,1,1], [1,1,0,0,0,1,1,0,1,1]]) 

if __name__ == "__main__":

    # read the raw IQ data from the recorded file
    # IQ data is normalized between -1 and +1 and interleaved
    # in_fname = "../data/iq_samples.raw"
    in_fname = "../data/samples_rds_1029.raw"
    iq_data = np.fromfile(in_fname, dtype='uint8')
    iq_data = (iq_data -128.0)/128.0
    print("Read raw RF data from \"" + in_fname + "\" in uint8 format. Block size is ", len(iq_data))
    iq_data = iq_data[0:8*307200]

    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, \
                                                    rf_Fc/(rf_Fs/2), \
                                                    window=('hann'))
    # coefficients for the filter to extract mono audio
    audio_coeff = signal.firwin(audio_taps, audio_Fc/(audio_Fs/2), window=('hann'))

    # set up drawing

    # select a block_size that is in KB and
    # a multiple of decimation factors
    block_size = 307200
    block_count = 0

    # states needed for continuity in block processing
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    state_phase = 0
    # add state as needed for the mono channel filter

    # ****************************************************************** 
    # -----------------------RDS States and Coeff-----------------------
    # ******************************************************************
    
    #Coefficents for extracting RDS data 
    extract_RDS_coeff = signal.firwin(rf_taps, [inital_RDS_lower_freq/(audio_Fs/2), inital_RDS_higher_freq/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
    pre_state_extract = np.zeros(rf_taps-1) 
    #BPF coefficents after squaring non-linearity
    square_coeff = signal.firwin(rf_taps, [squared_lower_freq/(audio_Fs/2), squared_higher_freq/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
    square_state =np.zeros(rf_taps-1)
    #Pll values 
    freq_centered =114000
    phase_adj = math.pi/3.3-math.pi/1.5
    state_Pll =[0.0, 0.0, 1.0, 0.0, 1.0, 0.0]

    #Coefficents for the LPF of fc 3000Hz
    lpf_coeff_rds = signal.firwin(rf_taps, cutoff_LPF/(audio_Fs/2), window=('hann'))
    lpf_3k_state = np.zeros(rf_taps-1) 
    lpf_3k_state_Q = np.zeros(rf_taps-1) 
    #Values for rational resampler
    upsample_val = 19
    downsample_val = 80 
    anti_img_coeff = signal.firwin(rf_taps, (57000/2)/((240000*19)/2), window=('hann'))
    anti_img_state = np.zeros(rf_taps-1)
    anti_img_state_Q = np.zeros(rf_taps-1)
    #Values for RRC
    rrc_Fs = 57000
    rrc_taps = 151 
    rrc_coeff = impulseResponseRootRaisedCosine(rrc_Fs, rrc_taps)
    rrc_state = np.zeros(rrc_taps -1) 
    rrc_state_Q = np.zeros(rrc_taps -1) 
    #Values for clock recoverey
    inital_offset = 0
    final_symb = 0 
    #differential decoding
    lonely_bit = 0 
    front_bit = 0
    remain_symbol = 0.0
    #Frame sync
    printposition = 0 
    prev_sync_bits = np.zeros(0) 
    last_position = -1
    # if the number of samples in the last block is less than the block size
    # it is fine to ignore the last few samples from the raw IQ file
    while (block_count+1)*block_size < len(iq_data):
        # if you wish to have shorter runtimes while troubleshooting
        # you can control the above loop exit condition as you see fit
        print("")
        print('Processing block ' + str(block_count))
        # filter to extract the FM channel (I samples are even, Q samples are odd)
        i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                        iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
                        zi=state_i_lpf_100k)
        q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                        iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
                        zi=state_q_lpf_100k)
        
        # downsample the FM channel
        i_ds = i_filt[::rf_decim]
        q_ds = q_filt[::rf_decim]

        # FM demodulator
        fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)

        # ****************************************************************** 
        # -----------------------RDS Data Processing------------------------ 
        # ******************************************************************
        
        # ------------------------Extraction--------------------------------
        # Performs convoloution to extract the data
        extract_rds, pre_state_extract = signal.lfilter(extract_RDS_coeff,1.0,fm_demod,zi=pre_state_extract)
        if block_count == 0: 
            print(extract_rds)

        # ---------------------Carrier Recovery-----------------------------
        #Squaring Nonolinearity
        #All this means is that we need to point why multiple each element by itself 
        squared_rds = np.square(extract_rds)

        #BPF
        pre_Pll_rds,square_state = signal.lfilter(square_coeff,1.0,squared_rds,zi=square_state)
        
        #PLL 
        post_Pll, post_Pll_Q, state_Pll =  fmPll(pre_Pll_rds, freq_centered, 240000, state_Pll, ncoScale = 0.5, phaseAdjust =phase_adj , normBandwidth = 0.001)

        # -----------------------Demodulation-------------------------------
        # Mixer 
        #Just a mixer which is just some good old point wise multiplication 
        #I component
        mixed_rds = np.multiply(extract_rds, post_Pll[0:len(extract_rds):1])*2
        #Q component
        mixed_rds_Q = np.multiply(extract_rds, post_Pll_Q[0:len(extract_rds):1])*2
        #print(mixed_rds)

        #LPF
        #I Compent 
        lpf_filt_rds,lpf_3k_state = signal.lfilter(lpf_coeff_rds,1.0, mixed_rds,zi=lpf_3k_state)
        #Q Compent 
        lpf_filt_rds_Q,lpf_3k_state_Q = signal.lfilter(lpf_coeff_rds,1.0, mixed_rds_Q,zi=lpf_3k_state_Q)
        #print(lpf_filt_rds)
        upsample_rds = np.zeros(len(lpf_filt_rds)*19) #Creates a list of empty zeros 
        upsample_rds_Q = np.zeros(len(lpf_filt_rds)*19) #Creates a list of empty zeros 

        #Rational Resampler
        #Upsamples by 19
        for i in range (len(lpf_filt_rds)):
            upsample_rds[i*upsample_val] = lpf_filt_rds[i]
            upsample_rds_Q[i*upsample_val] = lpf_filt_rds_Q[i]
        #Downsample by 19 in order to to get a frequency of 57kHz  
        #I Compent 
        anti_img,anti_img_state= signal.lfilter(anti_img_coeff,1.0, upsample_rds,zi=anti_img_state)
        #Q Compent 
        anti_img_Q,anti_img_state_Q= signal.lfilter(anti_img_coeff,1.0, upsample_rds_Q,zi=anti_img_state_Q)
        #Downsample by 19 in order to to get a frequency of 57kHz  
        resample_rds = anti_img[::downsample_val]*upsample_val
        resample_rds_Q = anti_img_Q[::downsample_val]*upsample_val

        #RRC Filter
        rrc_rds,rrc_state = signal.lfilter(rrc_coeff, 1.0, resample_rds,zi=rrc_state)
        #Q Component
        rrc_rds_Q, rrc_state_Q= signal.lfilter(rrc_coeff, 1.0, resample_rds_Q,zi=rrc_state_Q)

        #Clock and data recovery
        if block_count ==0:
            int_offset = (np.where(rrc_rds[0:24] == np.max(rrc_rds[0:24])))[0][0]
        #    if(abs(rrc_rds[0]) > abs(rrc_rds[12])):
        #        int_offset = 0;
        #    else: 
        #        int_offset = 12;
            print("Inital offset for clock recovery ", int_offset)

        #Go to every 24th sample 
        symbols_I = rrc_rds[int_offset::24]
        #print(len(symbols_I)) 
        symbols_Q = rrc_rds_Q[int_offset::24]
        #block processing the offset for the next block 
        int_offset = 24 - np.where(rrc_rds[len(rrc_rds)-24::] == symbols_I[-1])[0][0] 
        
        #Plotting
        if block_count == 0: 
            fig, (p_adjust1) = plt.subplots(nrows=1)
            fig, (rrc) = plt.subplots(nrows=1)
            fig.subplots_adjust(hspace = 1.0)
            p_adjust1.scatter(symbols_I, symbols_Q, s=10)
            p_adjust1.set_ylim(-1.25, 1.25)
            rrc.plot(rrc_rds[0:512], c = 'r') 
            rrc.plot(rrc_rds_Q[0:512], c = 'b') 

        # ---------------------RDS Data Processing--------------------------
        #Screening only needs to happen at block 0
        if block_count == 0:
            count_0_pos = 0
            count_1_pos = 0 
            for m in range(int(len(symbols_I)/4)): 
                #Counts the Amount of doubles when the start position is 0
                if((symbols_I[2*m] > 0 and symbols_I[2*m+1] > 0) or (symbols_I[2*m] < 0 and symbols_I[2*m+1]<0)):
                    count_0_pos += 1
                #Counts the Amount of doubles when the start position is 1 
                elif((symbols_I[2*m+1] > 0 and symbols_I[2*m+2] > 0) or (symbols_I[2*m+1] < 0 and symbols_I[2*m+2]<0)):
                    count_1_pos += 1
            
            print("Amount of doub when start 0 ", count_0_pos, " Amount of doub when 1 ", count_1_pos)
            # Now decide which point is 
            if(count_0_pos > count_1_pos): 
                start_pos = 1
            elif(count_1_pos > count_0_pos): 
                start_pos = 0
            print("Start position ", start_pos)
        
        bit_stream = np.zeros(int(len(symbols_I)/2)-start_pos)
        #Converts the bits
        flag_bit = 0 
        #Uses the last left over bit and first bit to get the proper value 
        if(start_pos == 1 and block_count != 0):
            if(lonely_bit > symbols_I[0]): 
                front_bit = 1
            elif(lonely_bit < symbols_I[0]): 
                front_bit = 0
        #Figures out what every bit is     
        for k in range(len(bit_stream)):
            #Now figure out what each bit is 
            if(start_pos+2*k+1 > len(symbols_I)-1):
                break
            #Should assign the bits 
            if(symbols_I[2*k+start_pos] > symbols_I[2*k+1+start_pos]): 
                bit_stream[k] = 1
            elif(symbols_I[2*k+start_pos] < symbols_I[2*k+1+start_pos]): 
                bit_stream[k] = 0
        #Sees if bit left over
        if(start_pos == 1): 
            #Should append to the front of array 
            #When this is added it breaks alot. which is weird 
            #It is appending correctly but for some reason it works v bad when we do append 
            bit_stream = np.insert(bit_stream, 0, front_bit, axis=0)
            lonely_bit = symbols_I[-1]
        #print("Len bit stream ",len(bit_stream))

        #Differential decoding
        if block_count == 0:
            prebit = bit_stream[0] 
            offset = 1
        else:
            offset = 0 
        diff_bits = np.zeros(len(bit_stream)-offset) 
        #Does XOR on the bits
        for t in range(len(diff_bits)): 
            diff_bits[t] = (prebit and not bit_stream[t+offset]) or (not prebit and bit_stream[t+offset])
            prebit = bit_stream[t+offset] 
        #Set first bit to last bit to then have XOR performed on it 
        prebit = bit_stream[-1]

        #Frame Sync and Error detection 
        #Need to check for sydromes: 
        if block_count != 0:
            diff_bits = np.insert(diff_bits, 0, prev_sync_bits, axis=0)
        #print("Len of diff bits ", len(diff_bits))
        position = 0 
        while True:
            block = diff_bits[position:position+26]
            #Below line for testing
            potential_syndrome = np.zeros(10) 
            #Does the binary matrix multiplication that uses XORs and ANDs
            for i in range(len(potential_syndrome)):
                    for j in range(26): 
                        mult = block[j] and H[j,i]
                        potential_syndrome[i] = (potential_syndrome[i] and not mult) or (not potential_syndrome[i] and mult)
            #convert to int
            potential_syndrome = potential_syndrome.astype(int)
            print(potential_syndrome)
            #Checks if syndrome A
            if ((potential_syndrome).tolist() == [1,1,1,1,0,1,1,0,0,0]):
                if(last_position == -1 or printposition-last_position == 26): 
                    last_position = printposition
                    print("Syndrome A at position ", printposition)
                    last_position = printposition
                else:
                    print("False positive Syndrome A at position ", printposition)
            #Checks if syndrome B
            elif ((potential_syndrome).tolist() == [1,1,1,1,0,1,0,1,0,0]):
                if(last_position == -1 or printposition-last_position == 26): 
                    print("Syndrome B at position ", printposition)
                    last_position = printposition
                else:
                    print("False positive Syndrome B at position ", printposition)
            #Checks if syndrome C
            elif ((potential_syndrome).tolist() == [1,0,0,1,0,1,1,1,0,0]):
                if(last_position == -1 or printposition-last_position == 26): 
                    print("Syndrome C at position ", printposition)
                    last_position = printposition
                else:
                    print("False positive Syndrome C at position ", printposition)
            #Checks if syndrome D
            elif ((potential_syndrome).tolist() == [1,0,0,1,0,1,1,0,0,0]):
                if(last_position == -1 or printposition-last_position == 26): 
                    print("Syndrome D at position ", printposition)
                    last_position = printposition
                else:
                    print("False positive Syndrome D at position ", printposition)
            #Breaks once it reaches the end
            position += 1 
            if( position+26 > len(diff_bits)-1):
                break
            printposition += 1 
        #Creates list of bits not used 
        prev_sync_bits = diff_bits[position-1::]
        #print(len(prev_sync_bits))

        #Iterates through the blocks 
        block_count += 1

# uncomment assuming you wish to show some plots
plt.show()
