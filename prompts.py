'''
This module finds prompt-like events.

'''
from common import *

_EMIN = 0.7
_EMAX = 12

def isPromptLike(event_detector, event_energy, x, y, z):
    correct_detector = event_detector in AD_DETECTORS
    correct_energy = (event_energy > _EMIN
            and event_energy < _EMAX)
    correct_position = (x**2 + y**2 < IAV_INNER_RADIUS**2
            and z**2 < IAV_INNER_HALF_HEIGHT**2)
    return correct_detector and correct_energy and correct_position

_NH_EMIN = 1.5
_NH_EMAX = 12

def isPromptLike_nH(event_detector, event_energy, x, y, z):
    correct_detector = event_detector in AD_DETECTORS
    correct_energy = (event_energy > _NH_EMIN
            and event_energy < _NH_EMAX)
    correct_position = (x**2 + y**2 < OAV_INNER_RADIUS**2
            and z**2 < OAV_INNER_HALF_HEIGHT**2)
    return correct_detector and correct_energy and correct_position
