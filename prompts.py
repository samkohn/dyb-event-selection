'''
This module finds prompt-like events.

'''
from common import *

_EMIN = 0.7
_EMAX = 12

def isPromptLike(event_detector, event_energy):
    return (event_detector in AD_DETECTORS
            and event_energy > _EMIN
            and event_energy < _EMAX)

_NH_EMIN = 1.5
_NH_EMAX = 12

def isPromptLike_nH(event_detector, event_energy):
    return (event_detector in AD_DETECTORS
            and event_energy > _NH_EMIN
            and event_energy < _NH_EMAX)
