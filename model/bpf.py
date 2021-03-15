import math

def bpf(fb, fe, fs, n_taps):
    n_center = ((fe + fb)/2)/(fs/2)
    n_pass = ((fe - fb))/(fs/2)
    for x in range(n_taps):
        if i==(n_taps-1/2):
            h[i] = n_pass
        else:
            h[i] = n_pass*((math.sin(math.pi*(n_pass/2)*(i-(n_taps-1)/2)))/(pi*(n_pass/2)*(i-(n_taps-1)/2)))
        h[i] = h[i]*math.cos(i*math.pi*n_center)
        h[i] = h[i]*math.pow(math.sin(i*pi/n_taps),2)

