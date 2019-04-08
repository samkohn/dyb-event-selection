'''
This module finds delayed-like events.

'''
from common import *

_EMIN = 6
_EMAX = 12

def isDelayedLike(event_detector, event_energy):
    return (event_detector in AD_DETECTORS
            and event_energy > _EMIN
            and event_energy < _EMAX)
