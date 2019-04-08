'''
This module finds delayed-like events.

'''

_EMIN = 6
_EMAX = 12
_AD_DETECTORS = [1, 2, 3, 4, 5, 6, 7, 8]

def isDelayedLike(event_detector, event_energy):
    return (event_detector in _AD_DETECTORS
            and event_energy > _EMIN
            and event_energy < _EMAX)
