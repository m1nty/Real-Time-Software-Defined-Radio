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
audio_Fs = 48e3

# complete your own settings for the mono channel
# (cutoff freq, audio taps, decimation rate, ...)
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5

if __name__ == "__main__":
    # read the raw IQ data from the recorded file
    # IQ data is normalized between -1 and +1 and interleaved
    # in_fname = "../data/iq_samples.raw"
    in_fname = "../data/test4.raw"
    iq_data = np.fromfile(in_fname, dtype='float32')
    print("Read raw RF data from \"" + in_fname + "\" in float32 format")

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
    #TODO mess with these values to see what works best 
    # note, we also need the Q component to properly tune the PLL using constalation diagrams 
    post_Pll, post_Pll_Q =  fmPll(pre_Pll_rds, freq_centered, square_Fs, ncoScale = 0.5, phaseAdjust = math.pi/2, normBandwidth = 0.01)

    # -----------------------Demodulation-------------------------------
    # Mixer 
    #Just a mixer which is just some good old point wise multiplication 
    #I component
    mixed_rds = np.multiply(extract_rds, post_Pll[1::])
    #Q component
    mixed_rds_Q = np.multiply(extract_rds, post_Pll_Q[1::])

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
    resample_rds = upsample_rds[::80]
    resample_rds_Q = upsample_rds_Q[::80]
    #After this both signals should be 57kHz

    #Root raised cosine filter
    #Black box in order to get the coefficents
    rrc_Fs = square_Fs * (19/80)
    rrc_taps = 151 
    rrc_coeff = impulseResponseRootRaisedCosine(rrc_Fs, rrc_taps)
    #I Component
    rrc_rds = signal.lfilter(rrc_coeff, 1.0, resample_rds)
    #Q Component
    rrc_rds_Q = signal.lfilter(rrc_coeff, 1.0, resample_rds_Q)

    #Clock and data recovery
    #Need to sample the shit
    #We assume 24 bits per symbol, look at each 24 samples to determine what the symbol is
    sample_per_symbol = 24 
    symbols_I = np.zeros(int(len(rrc_rds)/sample_per_symbol)) 
    symbols_Q = np.zeros(int(len(rrc_rds_Q)/sample_per_symbol)) 
    position = 0

    #Check each 24 samples to try and identify the symbols
    #TODO is 24 samples per symbol the H or the 1???
    symbols_I = rrc_rds[::24]
    symbols_Q = rrc_rds_Q[::24]

    #Plot scatter plots in order to tune the PLL 

    # ---------------------RDS Data Processing--------------------------
    #Deconding 

    #Frame Sync and Error detection 

    #RDS Application Layer

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
    plt.show()

    #plt.plot(10*pre_Pll_rds[10180:10200], c='b') 
    #plt.plot(post_Pll[10180:10200], c = 'r') 
    #plt.show()
    #plt.plot(rrc_rds, c = 'r') 
    #plt.show()
    # write audio data to file (assumes audio_data samples are -1 to +1)
    #wavfile.write("../data/fmMonoBasic.wav", int(audio_Fs), np.int16((audio_data/2)*32767))
    # during FM transmission audio samples in the mono channel will contain
    # the sum of the left and right audio channels; hence, we first
    # divide by two the audio sample value and then we rescale to fit
    # in the range offered by 16-bit signed int representation
