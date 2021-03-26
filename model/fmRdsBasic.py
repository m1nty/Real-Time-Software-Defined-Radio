#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math

# use "custom" fmDemodArctan
from fmSupportLib import fmDemodArctan
from fmSupportLib import my_convoloution
from fmSupportLib import my_filterImpulseResponse
from fmPll import fmPll
from fmRRC import impulseResponseRootRaisedCosine
# the radio-frequency (RF) sampling rate
# this sampling rate is either configured on RF hardware
# or documented when a raw file with IQ samples is provided
rf_Fs = 2.4e6

# the cutoff frequency to extract the FM channel from raw IQ data
rf_Fc = 100e3

# the number of taps for the low-pass filter to extract the FM channel
# this default value for the width of the impulse response should be changed
# depending on some target objectives, like the width of the transition band
# and/or the minimum expected attenuation from the pass to the stop band
rf_taps = 151

# the decimation rate when reducing the front end sampling rate (i.e., RF)
# to a smaller samping rate at the intermediate frequency (IF) where
# the demodulated data will be split into the mono/stereo/radio data channels
rf_decim = 10

# audio sampling rate (we assume audio will be at 48 KSamples/sec)
audio_Fs = 48000

# complete your own settings for the mono channel
# (cutoff freq, audio taps, decimation rate, ...)
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5
#Parity matrix
H = np.matrix([[1,0,0,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,1],[1,0,1,1,0,1,1,1,0,0],[0,1,0,1,1,0,1,1,1,0],[0,0,1,0,1,1,0,1,1,1],[1,0,1,0,0,0,0,1,1,1],[1,1,1,0,0,1,1,1,1,1],[1,1,0,0,0,1,0,0,1,1],[1,1,0,1,0,1,0,1,0,1],[1,1,0,1,1,1,0,1,1,0],[0,1,1,0,1,1,1,0,1,1],[1,0,0,0,0,0,0,0,0,1],[1,1,1,1,0,1,1,1,0,0], [0,1,1,1,1,0,1,1,1,0],[0,0,1,1,1,1,0,1,1,1],[1,0,1,0,1,0,0,1,1,1],[1,1,1,0,0,0,1,1,1,1], [1,1,0,0,0,1,1,0,1,1]]) 

if __name__ == "__main__":
    # read the raw IQ data from the recorded file
    # IQ data is normalized between -1 and +1 and interleaved
    # in_fname = "../data/iq_samples.raw"
    in_fname = "../data/test6.raw"
    iq_data = np.fromfile(in_fname, dtype='uint8')
    iq_data = (iq_data -128.0)/128.0
    print("Read raw RF data from \"" + in_fname + "\" in float32 format. Block size is ", len(iq_data))

    iq_data=iq_data[:8*307200]
    # Additional params needed for our own functions
    i_pre = np.zeros(rf_taps-1) 
    q_pre = np.zeros(rf_taps-1)
    audio_pre = np.zeros(audio_taps-1)

    #These are the in lab methods they are comments out in order to implement our own functions
    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

    # filter to extract the FM channel (I samples are even, Q samples are odd)
    i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
    q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])

    # downsample the FM channel
    i_ds = i_filt[::rf_decim]
    q_ds = q_filt[::rf_decim]

    # FM demodulator (check the library)
    fm_demod, dummy = fmDemodArctan(i_ds, q_ds)
    # we use a dummy because there is no state for this single-pass model

    # set up drawing
    fig, (ax0, ax1, ax2,ax4) = plt.subplots(nrows=4)
    fig.subplots_adjust(hspace = 1.0)

    # PSD after FM demodulation
    ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax0.set_ylabel('PSD (db/Hz)')
    ax0.set_title('Demodulated FM')

    # ****************************************************************** 
    # -----------------------------RDS----------------------------------
    # *****************************************************************

    # ------------------------Extraction--------------------------------
    #BPF filter to extract inital data
    #Use lfilter from scipy 
    inital_RDS_lower_freq = 54000
    inital_RDS_higher_freq = 60000
    # Intermidiate frequency sample rate 
    if_Fs = 240000

    # Creates BPF filter coefficents
    extract_RDS_coeff = signal.firwin(rf_taps, [inital_RDS_lower_freq/(if_Fs/2), inital_RDS_higher_freq/(if_Fs/2)], window=('hann'), pass_zero="bandpass")

    # Performs convoloution to extract the data
    extract_rds = signal.lfilter(extract_RDS_coeff,1.0,fm_demod)

    # ---------------------Carrier Recovery-----------------------------
    #Squaring Nonolinearity
    #All this means is that we need to point why multiple each element by itself 
    squared_rds = np.square(extract_rds)

    #TODO eventually try to combine this with the PLL to increase speed
    #BPF
    #Now that values
    squared_lower_freq = 113500
    squared_higher_freq = 114500
    square_Fs = 240000
    square_coeff = signal.firwin(rf_taps, [squared_lower_freq/(square_Fs/2), squared_higher_freq/(square_Fs/2)], window=('hann'), pass_zero="bandpass")
    pre_Pll_rds = signal.lfilter(square_coeff,1.0,squared_rds)

    #PLL
    #For RDS we need to output a bit more extra stuff
    #Know it it centered arounf 
    freq_centered = 114000
    # note, we also need the Q component to properly tune the PLL using constalation diagrams 
    # This works with what nicolici gave us 
    state_Pll =[0.0, 0.0, 1.0, 0.0, 1.0, 0.0]
    phase_adj = math.pi/3.3-math.pi/1.5
    post_Pll, post_Pll_Q, state_Pll =  fmPll(pre_Pll_rds, freq_centered, 240000, state_Pll, ncoScale = 0.5, phaseAdjust =phase_adj , normBandwidth = 0.001)

    # -----------------------Demodulation-------------------------------
    # Mixer 
    #Just a mixer which is just some good old point wise multiplication 
    #I component
    mixed_rds = np.multiply(extract_rds, post_Pll[0:len(extract_rds):1])*2
    #Q component
    mixed_rds_Q = np.multiply(extract_rds, post_Pll_Q[0:len(extract_rds):1])*2

    #LPF
    cutoff_LPF = 3000
    lpf_coeff_rds = signal.firwin(rf_taps, cutoff_LPF/(square_Fs/2), window=('hann'))
    #I Compent 
    lpf_filt_rds = signal.lfilter(lpf_coeff_rds,1.0, mixed_rds)
    #Q Compent 
    lpf_filt_rds_Q = signal.lfilter(lpf_coeff_rds,1.0, mixed_rds_Q)

    #Rational Resampler
    # to get 240000 to 57000 need to multiply by 19 / 80
    upsample_rds = np.zeros(len(lpf_filt_rds)*19) #Creates a list of empty zeros 
    upsample_rds_Q = np.zeros(len(lpf_filt_rds)*19) #Creates a list of empty zeros 

    #Upsamples by 19
    for i in range (len(lpf_filt_rds)):
        upsample_rds[i*19] = lpf_filt_rds[i]
        upsample_rds_Q[i*19] = lpf_filt_rds_Q[i]

    #Downsample by 19 in order to to get a frequency of 57kHz  
    #TODO Add anti imagining filter `
    anti_img_coeff = signal.firwin(rf_taps, (57000/2)/((240000*19)/2), window=('hann'))
    #I Compent 
    anti_img = signal.lfilter(anti_img_coeff,1.0, upsample_rds)
    #Q Compent 
    anti_img_Q = signal.lfilter(anti_img_coeff,1.0, upsample_rds_Q)
    #Downsample by 19 in order to to get a frequency of 57kHz  
    resample_rds = anti_img[::80]*19
    resample_rds_Q = anti_img_Q[::80]*19
    #After this both signals should be 57kHz

    #Root raised cosine filter
    #Black box in order to get the coefficents
    rrc_Fs = 57000
    rrc_taps = 151 
    rrc_coeff = impulseResponseRootRaisedCosine(rrc_Fs, rrc_taps)
    #I Component
    rrc_rds = signal.lfilter(rrc_coeff, 1.0, resample_rds)
    #Q Component
    rrc_rds_Q = signal.lfilter(rrc_coeff, 1.0, resample_rds_Q)
    print("RRC output this long ", len(rrc_rds))

    #Clock and data recovery
    #Need to sample the shit
    #Check each 24 samples to try and identify the symbols
    int_offset = 0 
    
    #Pretty sure this method is good 
    int_offset = (np.where(rrc_rds[0:24] == np.max(rrc_rds[0:24])))[0][0]

    print("int offset ", int_offset)
    #Go to every 24th sample 
    symbols_I = rrc_rds[int_offset::24]
    symbols_Q = rrc_rds_Q[int_offset::24]

    # ---------------------RDS Data Processing--------------------------
    #Decoding 
    #Now that we have the symbols for high and low, we can then figure out the bits each represent
    bit_stream = np.zeros(int(len(symbols_I)/2))
    count_0_pos = 0
    count_1_pos = 0 

    #Screening
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

    #Sets the values 
    print("Lenth of symbols ",len(symbols_I) ," Length bit stream", len(bit_stream))
    for k in range(len(bit_stream)):
        #Now figure out what each bit is 
        if(start_pos+2*k+1 > len(symbols_I)-1):
            break
        #Should assign the bits 
        if(symbols_I[2*k+start_pos] > symbols_I[2*k+1+start_pos]): 
            bit_stream[k] = 1
        elif(symbols_I[2*k+start_pos] < symbols_I[2*k+1+start_pos]): 
            bit_stream[k] = 0
        #print(bit_stream[k], " and ", k)

    #Differential decoding
    diff_bits = np.zeros(len(bit_stream)-1) 
    prebit = bit_stream[0] 

    #Does XOR on the bits
    for t in range(len(diff_bits)): 
        diff_bits[t] = (prebit and not bit_stream[t+1]) or (not prebit and bit_stream[t+1])
        prebit = bit_stream[t+1] 
    print("Size of diff bit ", len(diff_bits))

    #NOTE This shit slaps, problem is not here 
    #Frame Sync and Error detection 
    #Need to check for sydromes: 
    position = 0 
    while True:
        block = diff_bits[position:position+26]
        #Below line used to test if this shit works correctly, and it do
        potential_syndrome = np.zeros(10) 
        #Does the binary matrix multiplication that uses XORs and ANDs
        for i in range(len(potential_syndrome)):
                for j in range(26): 
                    mult = block[j] and H[j,i]
                    potential_syndrome[i] = (potential_syndrome[i] and not mult) or (not potential_syndrome[i] and mult)
        #convert to int
        potential_syndrome = potential_syndrome.astype(int)
        #print(potential_syndrome)
        if ((potential_syndrome).tolist() == [1,1,1,1,0,1,1,0,0,0]):
            print("Syndrome A at position ", position)
        elif ((potential_syndrome).tolist() == [1,1,1,1,0,1,0,1,0,0]):
            print("Syndrome B at position ", position)
        elif ((potential_syndrome).tolist() == [1,0,0,1,0,1,1,1,0,0]):
            print("Syndrome C at position ", position)
        elif ((potential_syndrome).tolist() == [1,0,0,1,0,1,1,0,0,0]):
            print("Syndrome D at position ", position)
        #Breaks once it reaches the end
        position += 1 
        if( position+26 > len(diff_bits)):
            break

    #Retrieve Radio text
    #Once it is all synced up, we can then start to recieve radio text from the blocks 

    # save PSD plots
    ax1.psd(pre_Pll_rds, NFFT=512, Fs=(square_Fs/rf_decim)/1e3)
    ax1.set_ylabel('PSD (db/Hz)')
    ax1.set_title('Pre Pll')

    ax2.psd(post_Pll, NFFT=512, Fs=(square_Fs/rf_decim)/1e3)
    ax2.set_ylabel('PSD (db/Hz)')
    ax2.set_title('Post Pll')

    ax4.psd(rrc_rds, NFFT=512, Fs=(rrc_Fs/1e3))
    ax4.set_ylabel('PSD (db/Hz)')
    ax4.set_title('Post RRC Filter')

    fig.savefig("../data/fmRdsBasic.png")
    plt.show()

    #Plot scatter plots in order to tune the PLL 
    fig, (p_adjust1) = plt.subplots(nrows=1)
    fig.subplots_adjust(hspace = 1.0)
    p_adjust1.scatter(symbols_I, symbols_Q, s=10)
    p_adjust1.set_ylim(-1.25, 1.25)
    plt.show()

    #plt.plot(10*pre_Pll_rds[10180:10200], c='b') 
    #plt.plot(post_Pll[10180:10200], c = 'r') 
    #plt.show()
    plt.plot(rrc_rds[10000:10240], c = 'r') 
    plt.plot(rrc_rds_Q[10000:10240], c = 'b') 
    plt.show()
    # write audio data to file (assumes audio_data samples are -1 to +1)
    #wavfile.write("../data/fmMonoBasic.wav", int(audio_Fs), np.int16((audio_data/2)*32767))
    # during FM transmission audio samples in the mono channel will contain
    # the sum of the left and right audio channels; hence, we first
    # divide by two the audio sample value and then we rescale to fit
    # in the range offered by 16-bit signed int representation
