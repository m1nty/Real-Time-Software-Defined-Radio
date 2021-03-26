
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
    in_fname = "../data/test6.raw"
    iq_data = np.fromfile(in_fname, dtype='uint8')
    iq_data = (iq_data -128.0)/128.0
    print("Read raw RF data from \"" + in_fname + "\" in float32 format. Block size is ", len(iq_data))

    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, \
                                                    rf_Fc/(rf_Fs/2), \
                                                    window=('hann'))
    # coefficients for the filter to extract mono audio
    audio_coeff = signal.firwin(audio_taps, audio_Fc/(audio_Fs/2), window=('hann'))

    # set up drawing
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
    fig.subplots_adjust(hspace = 1.0)

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
    state_Pll = np.zeros(6)
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

    #Frame sync
    
    # if the number of samples in the last block is less than the block size
    # it is fine to ignore the last few samples from the raw IQ file
    while (block_count+1)*block_size < len(iq_data):
        # if you wish to have shorter runtimes while troubleshooting
        # you can control the above loop exit condition as you see fit

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
        extract_rds, pre_state_extract = signal.lfilter(extract_RDS_coeff,1.0,fm_demod,pre_state_extract)

        # ---------------------Carrier Recovery-----------------------------
        #Squaring Nonolinearity
        #All this means is that we need to point why multiple each element by itself 
        squared_rds = np.square(extract_rds)

        #BPF
        pre_Pll_rds,square_state = signal.lfilter(square_coeff,1.0,squared_rds,square_state)
        
        #PLL 
        #NOTE Needa ass shit
        post_Pll, post_Pll_Q =  fmPll(pre_Pll_rds, freq_centered, 240000, ncoScale = 0.5, phaseAdjust =phase_adj , normBandwidth = 0.001)

        # -----------------------Demodulation-------------------------------
        # Mixer 
        #Just a mixer which is just some good old point wise multiplication 
        #I component
        mixed_rds = np.multiply(extract_rds, post_Pll[0:len(extract_rds):1])*2
        #Q component
        mixed_rds_Q = np.multiply(extract_rds, post_Pll_Q[0:len(extract_rds):1])*2

        #LPF
        #I Compent 
        lpf_filt_rds,lpf_3k_state = signal.lfilter(lpf_coeff_rds,1.0, mixed_rds,lpf_3k_state)
        #Q Compent 
        lpf_filt_rds_Q,lpf_3k_state_Q = signal.lfilter(lpf_coeff_rds,1.0, mixed_rds_Q,lpf_3k_state_Q)

        upsample_rds = np.zeros(len(lpf_filt_rds)*19) #Creates a list of empty zeros 
        upsample_rds_Q = np.zeros(len(lpf_filt_rds)*19) #Creates a list of empty zeros 

        #Rational Resampler
        #Upsamples by 19
        for i in range (len(lpf_filt_rds)):
            upsample_rds[i*upsample_val] = lpf_filt_rds[i]
            upsample_rds_Q[i*upsample_val] = lpf_filt_rds_Q[i]
        #Downsample by 19 in order to to get a frequency of 57kHz  
        #I Compent 
        anti_img,anti_img_state= signal.lfilter(anti_img_coeff,1.0, upsample_rds,anti_img_state)
        #Q Compent 
        anti_img_Q,anti_img_state_Q= signal.lfilter(anti_img_coeff,1.0, upsample_rds_Q,anti_img_state_Q)
        #Downsample by 19 in order to to get a frequency of 57kHz  
        resample_rds = anti_img[::downsample_val]*upsample_val
        resample_rds_Q = anti_img_Q[::downsample_val]*upsample_val

        #RRC Filter
        rrc_rds,rrc_state = signal.lfilter(rrc_coeff, 1.0, resample_rds,rrc_state)
        #Q Component
        rrc_rds_Q, rrc_state_Q= signal.lfilter(rrc_coeff, 1.0, resample_rds_Q,rrc_state_Q)

        #Clock and data recovery
        #TODO add some fancy shit for block processing 
        #Need to sample the shit
        #Check each 24 samples to try and identify the symbols
        #Pretty sure this method is good 
        if block_count == 0: 
            int_offset = (np.where(rrc_rds[0:24] == np.max(rrc_rds[0:24])))[0][0]
        #Go to every 24th sample 
        symbols_I = rrc_rds[int_offset::24]
        symbols_Q = rrc_rds_Q[int_offset::24]

        block_count += 1



# uncomment assuming you wish to show some plots
# plt.show()
