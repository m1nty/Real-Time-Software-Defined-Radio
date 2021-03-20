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
audio_Fs = 240e3

# complete your own settings for the mono channel
# (cutoff freq, audio taps, decimation rate, ...)
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5
block_size = 1024*rf_decim*audio_decim*2


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
    # i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
    # q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])

    i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0:block_size:2])
    q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1:block_size:2])

    print("Done all coeff/iq")


    #-----------------TAKE HOME EXERCISE------------------
    # rf_coeff = my_filterImpulseResponse(rf_Fc, rf_Fs, rf_taps)
    # print("Completed Creating RF Coeff")
    # i_filt, dummy_i = my_convoloution(iq_data[0::2], rf_coeff, rf_taps, i_pre) 
    # print("Completed i conv")
    # q_filt, dummy_q = my_convoloution(iq_data[1::2], rf_coeff, rf_taps, q_pre) 
    # print("Completed q conv")

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

    #-------------------IN LAB EXPERIMENT--------------------------
    # coefficients for the filter to extract mono audio
    # For in lab experiment, the scipy firwin function can be utilized 
    audio_coeff = signal.firwin(audio_taps, audio_Fc/(audio_Fs/2), window=('hann'))

    # extract the mono audtio data through filtering
    # For the in lab experiment, the scipy function lfilter can be used
    audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)  

    #-----------------TAKE HOME EXERCISE------------------
    # audio_coeff = my_filterImpulseResponse(audio_Fc, audio_Fs, audio_taps)
    # print("Completed creating audio coeff") 
    # audio_filt, dummy_audio = my_convoloution(fm_demod, audio_coeff, audio_taps, audio_pre) 
    # print("Completed audio conv") 

    # PSD after extracting mono audio
    ax1.psd(audio_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax1.set_ylabel('PSD (db/Hz)')
    ax1.set_title('Extracted Mono')

    # downsample audio data
    audio_data = audio_filt[::audio_decim] 

    # PSD after decimating mono audio
    ax2.psd(audio_data, NFFT=512, Fs=audio_Fs/1e3)
    ax2.set_ylabel('PSD (db/Hz)')
    ax2.set_title('Mono Audio')

    # save PSD plots
    fig.savefig("../data/fmMonoBasic.png")
    # write audio data to file (assumes audio_data samples are -1 to +1)

    # wavfile.write("../data/fmMonoBasic.wav", int(audio_Fs), np.int16((audio_data/2)*32767))

    # during FM transmission audio samples in the mono channel will contain
    # the sum of the left and right audio channels; hence, we first
    # divide by two the audio sample value and then we rescale to fit
    # in the range offered by 16-bit signed int representation
    
    print("Stereo")

    #-----------------------STEREO CARRIER RECOVERY-------------------------------
    bandPassCoeff = signal.firwin(rf_taps, [18.5e3/(audio_Fs/2), 19.5e3/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
    bpf_recovery = signal.lfilter(bandPassCoeff, 1.0, fm_demod)
    print("Completed Recovery BPF Convolution")
    recovery_pll = fmPll(bpf_recovery, 19e3, 240e3,2)

    #-----------------------STEREO CHANNEL EXTRACTION-------------------------------
    bandPassCoeff = signal.firwin(rf_taps,[22e3/(audio_Fs/2), 54e3/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
    bpf_extraction = signal.lfilter(bandPassCoeff, 1.0, fm_demod)
    print("Completed Extraction BPF Convolution")

    #-----------------------STEREO PROCESSING-------------------------------
    #Mixing
    mixed = np.multiply(recovery_pll[0:len(bpf_extraction):1],bpf_extraction)
    #LPF
    stereo_coeff = signal.firwin(rf_taps, 16e3/(audio_Fs/2), window=('hann'))
    stereo_filt = signal.lfilter(stereo_coeff, 1.0, mixed)  
    #Decimation
    stereo_data = stereo_filt[::5] 

    #Combiner
    combined_l, combined_r = audio_data, audio_data

    for i in range(len(audio_data)):
        combined_l[i] = (audio_data[i]+stereo_data[i])/2
        combined_r[i] = (audio_data[i]-stereo_data[i])/2


    fig, (ax0,ax1,ax2,ax3,ax4) = plt.subplots(nrows=5)
    fig.subplots_adjust(hspace = 1.0)

    ax0.psd(bpf_recovery, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax0.set_ylabel('PSD (db/Hz)')
    ax0.set_title('Extracted Stereo Pilot Tone')

    ax1.psd(bpf_extraction, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax1.set_ylabel('PSD (db/Hz)')
    ax1.set_title('Extracted Stereo Channel')

    ax2.psd(mixed, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax2.set_ylabel('PSD (db/Hz)')
    ax2.set_title('cos(α - β) and cos(α + β) Components')

    ax3.psd(stereo_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax3.set_ylabel('PSD (db/Hz)')
    ax3.set_title('LPF cos(α - β) Extraction')

    ax4.psd(audio_data, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax4.psd(stereo_data, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax4.set_ylabel('PSD (db/Hz)')
    ax4.set_title('Mono Audio vs Stereo Audio')

    fig, (ax5,ax6) = plt.subplots(nrows=2)
    fig.subplots_adjust(hspace = 1.0)
    ax5.psd(combined_l, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax5.set_ylabel('PSD (db/Hz)')
    ax5.set_title('Left Audio Channel')

    ax6.psd(combined_r, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
    ax6.set_ylabel('PSD (db/Hz)')
    ax6.set_title('Right Audio Channel')

    # combined_2D = np.array( audio_data.tolist(), audio_data.tolist())
    # wavfile.write("../data/fmStereo.wav", int(audio_Fs), combined_2D)

    plt.show()

    # # save PSD plots
    # fig.savefig("../data/fmStereoBasic.png")
    # plt.show()
