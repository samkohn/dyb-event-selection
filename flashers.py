'''
This module computes flasher cuts.

'''
from math import log10, pow

def fID(fMax, fQuad):
    if fMax == 0 and fQuad == 0:
        return 3
    else:
        return log10(pow(fQuad, 2) + pow(fMax/0.45, 2))

def fPSD(f_t1, f_t2):
    if f_t1 == 1 and f_t2 == 1:
        return 3
    else:
        t1_term = 4*pow(1-f_t1, 2)
        t2_term = 1.8*pow(1-f_t2, 2)
        return log10(t1_term + t2_term)

def isFlasher(event_fID, event_fPSD):
    return (event_fID >= 0 or event_fPSD >= 0)
