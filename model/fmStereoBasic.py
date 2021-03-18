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

from fmPll import fmPll

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

    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))
    print("Completed Creating RF Coeff")

    # filter to extract the FM channel (I samples are even, Q samples are odd)
    i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
    print("Completed i conv")
    q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])
    print("Completed q conv")


    # downsample the FM channel
    i_ds = i_filt[::rf_decim]
    q_ds = q_filt[::rf_decim]

    # FM demodulator (check the library)
    fm_demod, dummy = fmDemodArctan(i_ds, q_ds)
    # we use a dummy because there is no state for this single-pass model

    # set up drawing
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
    fig.subplots_adjust(hspace = 1.0)

    # PSD after FM demodulation
    ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax0.set_ylabel('PSD (db/Hz)')
    ax0.set_title('Demodulated FM')



    #-----------------------STEREO CARRIER RECOVERY-------------------------------
	bandPassCoeff = signal.firwin(rf_taps, [18.5e3, 19.5e3], window=('hann'), pass_zero="bandpass")
    bpf_recovery = signal.lfilter(rf_coeff, 1.0, fm_demod)
    print("Completed Recovery BPF Convolution")

    # PSD after extracting Pilot Tone
    ax1.psd(bpf_recovery, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax1.set_ylabel('PSD (db/Hz)')
    ax1.set_title('Extracted Stereo Pilot Tone')
    
    #PLL/NCO
    #def fmPll(pllIn, freq, Fs, ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01):
    recovery_pll = fmPll(bpf_recovery, 19e3, 2400e3)

    
    #-----------------------STEREO CHANNEL EXTRACTION-------------------------------
	bandPassCoeff = signal.firwin(rf_taps, [22e3,54e3], window=('hann'), pass_zero="bandpass")
    bpf_extraction = signal.lfilter(bandPassCoeff, 1.0, fm_demod)
    print("Completed Extraction BPF Convolution")

    # PSD after extracting Stereo Channel
    ax1.psd(bpf_extraction, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax1.set_ylabel('PSD (db/Hz)')
    ax1.set_title('Extracted Stereo Channel')

    #-----------------------STEREO PROCESSING-------------------------------
    #Mixing
    for x in range(len(recovery_pll)):
        mixed = recovery_pll[x]*bpf_extraction[x]
    #LPF
    stereo_coeff = signal.firwin(rf_taps, 16e3, window=('hann'))
    stereo_filt = signal.lfilter(stereo_coeff, 1.0, mixed)  
    #Decimation
    stereo_data = stereo_filt[::5] 

    # PSD after extracting stereo channel
    ax1.psd(stereo_data, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax1.set_ylabel('PSD (db/Hz)')
    ax1.set_title('Extracted Stereo')
    

    

    
    







    # write audio data to file (assumes audio_data samples are -1 to +1)
    wavfile.write("../data/fmMonoStereo.wav", int(audio_Fs), np.int16((audio_data/2)*32767))
    # during FM transmission audio samples in the mono channel will contain
    # the sum of the left and right audio channels; hence, we first
    # divide by two the audio sample value and then we rescale to fit
    # in the range offered by 16-bit signed int representation
