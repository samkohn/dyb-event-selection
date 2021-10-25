'''
This module finds AD events.

'''
from common import *

_EMIN = 0.65

def isAnyADEvent(event_detector):
    return event_detector in AD_DETECTORS

def isADEvent(event_detector, event_energy):
    return (event_detector in AD_DETECTORS
            and event_energy >= _EMIN)

_EMIN_THU = 1.5
_EMAX_THU = 12
def isADEvent_THU(event_detector, event_energy):
    return (event_detector in AD_DETECTORS
            and event_energy >= _EMIN_THU)

_EMIN_LOW_THU = 0.7
def isADEvent_THU_lowenergy(event_detector, event_energy):
    return (event_detector in AD_DETECTORS
            and event_energy >= _EMIN_LOW_THU)

_COINCIDENCE_TIME_NS = int(9e5)
def hasCoincidence(event_isADEvent, dt_last_ADevent, dt_next_ADevent,
        event_isFlasher):
    return (event_isFlasher in (0, 2)
            and event_isADEvent
            and (dt_last_ADevent < _COINCIDENCE_TIME_NS
                or dt_next_ADevent < _COINCIDENCE_TIME_NS))

_HIGH_EMIN = 12
def hasHighEnergy(event_isADEvent, event_energy):
    return event_isADEvent and (event_energy > _HIGH_EMIN)
