import numpy as np
import math

def fmPll(pllIn, freq, Fs, recovery_state, ncoScale = 1.0, phaseAdjust = 0.0, \
    normBandwidth = 0.01):

    Cp = 2.666
    Ci = 3.555
    Kp = (normBandwidth)*Cp
    Ki = (normBandwidth*normBandwidth)*Ci

    ncoOut = np.empty(len(pllIn)+1)
    ncoOutQ = np.empty(len(pllIn)+1)

    integrator = recovery_state[0]
    phaseEst = recovery_state[1]
    feedbackI = recovery_state[2]
    feedbackQ = recovery_state[3]
    ncoOut[0] = recovery_state[4]
    trigOffset = recovery_state[5]

    for k in range(len(pllIn)):
        # phase detector
        errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
        errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential
        # four-quadrant arctangent discriminator for phase error detection
        errorD = math.atan2(errorQ, errorI)
        # loop filter
        integrator = integrator + Ki*errorD
        # update phase estimate
        phaseEst = phaseEst + Kp*errorD + integrator
        # internal oscillator
        trigArg = 2*math.pi*(freq/Fs)*(trigOffset+k+1) + phaseEst
        feedbackI = math.cos(trigArg)
        feedbackQ = math.sin(trigArg)
        ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
        ncoOutQ[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)

    recovery_state[0] = integrator
    recovery_state[1] = phaseEst
    recovery_state[2] = feedbackI
    recovery_state[3] = feedbackQ
    recovery_state[4] = ncoOut[-1]
    recovery_state[5] = trigOffset + len(pllIn)

    return ncoOut, ncoOutQ,recovery_state
    # for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state
    # for RDS add also the quadrature NCO component to the output

if __name__ == "__main__":
    pass
