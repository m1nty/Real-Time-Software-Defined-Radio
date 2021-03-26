
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
    in_fname = "../data/test5.raw"
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
    freq_centered = 114000
    block_count = 0

    # states needed for continuity in block processing
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    state_phase = 0
    # add state as needed for the mono channel filter

    #Coefficents for extracting RDS data 
    extract_RDS_coeff = signal.firwin(rf_taps, [inital_RDS_lower_freq/(audio_Fs/2), inital_RDS_higher_freq/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
    #BPF coefficents after squaring non-linearity
    square_coeff = signal.firwin(rf_taps, [squared_lower_freq/(audio_Fs/2), squared_higher_freq/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
    #Coefficents for the LPF of fc 3000Hz
    lpf_coeff_rds = signal.firwin(rf_taps, cutoff_LPF/(audio_Fs/2), window=('hann'))

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

        block_count += 1

    print('Finished processing the raw I/Q samples')


# uncomment assuming you wish to show some plots
# plt.show()
