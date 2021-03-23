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

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 240e3
audio_decim = 5
# add other settings for audio, like filter taps, ...
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5

if __name__ == "__main__":

    # read the raw IQ data from the recorded file
    # IQ data is normalized between -1 and +1 and interleaved
    in_fname = "../data/test1.raw"
    iq_data = np.fromfile(in_fname, dtype='float32')
    print("Read raw RF data from \"" + in_fname + "\" in float32 format")

    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))
    # coefficients for the filter to extract mono audio
    audio_coeff = signal.firwin(audio_taps, audio_Fc/(audio_Fs/2), window=('hann'))

    # set up drawing
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
    fig.subplots_adjust(hspace = 1.0)

    # select a block_size that is in KB and
    # a multiple of decimation factors
    block_size = 1024 * rf_decim * audio_decim * 2
    block_count = 0

    # states needed for continuity in block processing
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    state_phase = 0

    state_recovery = np.zeros(rf_taps-1)
    state_extraction = np.zeros(rf_taps-1)
    # add state as needed for the mono channel filter

    # audio buffer that stores all the audio blocks
    audio_data = np.array([]) # to be updated by you during in-lab
    audio_pre_state = np.zeros(audio_taps-1)
    stereo_data = np.array([])
    stereo_pre_state = np.zeros(audio_taps-1)

    # INITIALIZE STATE VARIABLES FOR FIRST BLOCK
    # recovery_state;
    # [integrator, phaseEst, feedbackI, feedbackQ, ncoOut[0], trigOffset]
    recovery_state = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0]

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

        # extract the mono audtio data through filtering
        audio_filt, audio_pre_state = signal.lfilter(audio_coeff, 1.0, fm_demod,zi=audio_pre_state)

        # downsample audio data
        # audio_block = ... change as needed
        audio_block = audio_filt[::audio_decim]

        # concatanete most recently processed audio_block
        # to the previous blocks stored in audio_data
        audio_data = np.concatenate((audio_data, audio_block))



        #-----------------------STEREO CARRIER RECOVERY-------------------------------
        # 1. firwin
        bpcoeff_recovery = signal.firwin(rf_taps, [18.5e3/(audio_Fs/2), 19.5e3/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
        # 2. lfilter
        bpf_recovery, state_recovery = signal.lfilter(bpcoeff_recovery, 1.0, fm_demod, zi=state_recovery)
        # 3. pll
        recovery_pll, recovery_state = fmPll(bpf_recovery, 19e3, 240e3, recovery_state, 2)

        if block_count == 3:
            final_prev = recovery_pll[-10:len(recovery_pll)]

        if block_count == 4:
            indices = np.arange(0, 20)
            first_next = recovery_pll[0:10]

            win, (td_plot) = plt.subplots(1)
            win.suptitle("Time Domain Plot of End + Start Window")

            xs = indices
            print(final_prev.shape)
            print(first_next.shape)
            ys = np.concatenate([final_prev, first_next])

            print(xs, ys)
            print(recovery_pll)
            print(type(recovery_pll))

            td_plot.plot(xs, ys)
            plt.show()

            td_plot.set(xlabel="Sample #", ylabel="Amplitude Value")

        #recovery PLL block code
        print("Completed Recovery BPF Convolution")


        #-----------------------STEREO CHANNEL EXTRACTION-------------------------------
        bpcoeff_extraction = signal.firwin(rf_taps,[22e3/(audio_Fs/2), 54e3/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
        bpf_extraction, state_extraction = signal.lfilter(bpcoeff_extraction, 1.0, fm_demod,zi=state_extraction)

        #-----------------------STEREO PROCESSING-------------------------------
        #Mixing
        mixed = np.multiply(recovery_pll[0:len(bpf_extraction):1],bpf_extraction)
        mixed = mixed*2

        #LPF
        stereo_coeff = signal.firwin(rf_taps, 16e3/(audio_Fs/2), window=('hann'))
        stereo_filt, stereo_pre_state = signal.lfilter(stereo_coeff, 1.0, fm_demod,zi=stereo_pre_state)
        #Decimation
        stereo_block = stereo_filt[::5]
        stereo_data = np.concatenate((stereo_data, stereo_block))

        #Combiner
        combined_l, combined_r = audio_data, audio_data

        for i in range(len(audio_data)):
            combined_l[i] = (audio_data[i]+stereo_data[i])/2
            combined_r[i] = (audio_data[i]-stereo_data[i])/2

        # to save runtime select the range of blocks to log iq_data
        # this includes both saving binary files as well plotting PSD
        # below we assume we want to plot for graphs for blocks 10 and 11
        if block_count >= 10 and block_count < 12:
            # PSD after FM demodulation
            ax0.clear()
            ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
            ax0.set_ylabel('PSD (dB/Hz)')
            ax0.set_xlabel('Freq (kHz)')
            ax0.set_title('Demodulated FM (block ' + str(block_count) + ')')
            # output binary file name (where samples are written from Python)
            fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
            # create binary file where each sample is a 32-bit float
            fm_demod.astype('float32').tofile(fm_demod_fname)

            # PSD after extracting mono audio
            ax1.clear()
            ax1.psd(audio_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
            ax1.set_ylabel('PSD (dB/Hz)')
            ax1.set_xlabel('Freq (kHz)')
            ax1.set_title('Extracted Mono')

            # PSD after decimating mono audio
            ax2.clear()
            ax2.psd(audio_block, NFFT=512, Fs=audio_Fs/1e3)
            ax2.set_ylabel('PSD (dB/Hz)')
            ax2.set_xlabel('Freq (kHz)')
            ax2.set_title('Mono Audio')

            # save figure to file
            fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

            #STEREO PLOTS
            fig, (ax0,ax1,ax2,ax3,ax4) = plt.subplots(nrows=5)
            fig.subplots_adjust(hspace = 1.0)

            ax0.psd(bpf_recovery, NFFT=512, Fs=(audio_Fs/audio_decim)/1e3)
            ax0.set_ylabel('PSD (db/Hz)')
            ax0.set_title('Extracted Stereo Pilot Tone')

            ax1.psd(bpf_extraction, NFFT=512, Fs=(audio_Fs/audio_decim)/1e3)
            ax1.set_ylabel('PSD (db/Hz)')
            ax1.set_title('Extracted Stereo Channel')

            ax2.psd(mixed, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
            ax2.set_ylabel('PSD (db/Hz)')
            ax2.set_title('cos(α - β) and cos(α + β) Components')

            ax3.psd(stereo_filt, NFFT=512, Fs=(audio_Fs/audio_decim)/1e3)
            ax3.set_ylabel('PSD (db/Hz)')
            ax3.set_title('LPF cos(α - β) Extraction')

            ax4.psd(audio_data, NFFT=512, Fs=(audio_Fs/audio_decim)/1e3)
            ax4.psd(stereo_data, NFFT=512, Fs=(audio_Fs/audio_decim)/1e3)
            ax4.set_ylabel('PSD (db/Hz)')
            ax4.set_title('Mono Audio vs Stereo Audio')

            fig, (ax5,ax6) = plt.subplots(nrows=2)
            fig.subplots_adjust(hspace = 1.0)
            ax5.psd(combined_l, NFFT=512, Fs=(audio_Fs/audio_decim)/1e3)
            ax5.set_ylabel('PSD (db/Hz)')
            ax5.set_title('Left Audio Channel')

            ax6.psd(combined_r, NFFT=512, Fs=(audio_Fs/audio_decim)/1e3)
            ax6.set_ylabel('PSD (db/Hz)')
            ax6.set_title('Right Audio Channel')

            print(rf_Fs/rf_decim)

        block_count += 1

    print('Finished processing the raw I/Q samples')

# write audio data to a .wav file (assumes audio_data samples are -1 to +1)

    combined_l = np.int16((combined_l)*32767)
    combined_r = np.int16((combined_r)*32767)
    stereo = np.array([combined_l, combined_r]).transpose()

    print("Writing Wav File")
    wavfile.write("../data/fmMonoBlock.wav", int(48e3), stereo)

# uncomment assuming you wish to show some plots
plt.show()
