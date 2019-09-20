'''
This module finds delayed-like events.

'''
from common import *

_EMIN = 6
_EMAX = 12
_IBD_DT_MIN = int(1e3)
_IBD_DT_MAX = int(200e3)
_MULT_POST_MIN = int(200e3)

def isDelayedLike(event_detector, event_energy, x, y, z):
    correct_detector = event_detector in AD_DETECTORS
    correct_energy = (event_energy > _EMIN
            and event_energy < _EMAX)
    correct_position = (x**2 + y**2 < IAV_INNER_RADIUS**2
            and z**2 < IAV_INNER_HALF_HEIGHT**2)
    return correct_detector and correct_energy and correct_position

def isIBDDelayed(tag_DelayedLike,
        dt_previous_PromptLike, num_PromptLikes_400us,
        dt_next_DelayedLike, tag_WSMuonVeto, tag_ADMuonVeto,
        tag_ShowerMuonVeto, tag_flasher):
    if not tag_DelayedLike:
        return False
    if tag_flasher:
        return False
    if tag_WSMuonVeto or tag_ADMuonVeto or tag_ShowerMuonVeto:
        return False
    passesDtCut = (dt_previous_PromptLike > _IBD_DT_MIN
            and dt_previous_PromptLike < _IBD_DT_MAX)
    passesMultPre = num_PromptLikes_400us == 1
    passesMultPost = (dt_next_DelayedLike > _MULT_POST_MIN
            or dt_next_DelayedLike == -1)
    return (passesDtCut and passesMultPre and passesMultPost)

_NH_EMIN = 1.9
_NH_EMAX = 2.74
_NH_IBD_DT_MIN = int(1e3)
_NH_IBD_DT_MAX = int(400e3)
_NH_MULT_POST_MIN = int(400e3)

def isDelayedLike_nH(event_detector, event_energy, x, y, z):
    correct_detector = event_detector in AD_DETECTORS
    correct_energy = (event_energy > _NH_EMIN
            and event_energy < _NH_EMAX)
    correct_position = (x**2 + y**2 < OAV_INNER_RADIUS**2
            and z**2 < OAV_INNER_HALF_HEIGHT**2)
    return correct_detector and correct_energy and correct_position

def isIBDDelayed_nH(tag_DelayedLike,
        dt_previous_PromptLike, num_PromptLikes_800us,
        dt_next_DelayedLike, tag_WSMuonVeto, tag_ADMuonVeto,
        tag_ShowerMuonVeto, tag_flasher):
    if not tag_DelayedLike:
        return False
    if tag_flasher:
        return False
    if tag_WSMuonVeto or tag_ADMuonVeto or tag_ShowerMuonVeto:
        return False
    passesDtCut = (dt_previous_PromptLike > _NH_IBD_DT_MIN
            and dt_previous_PromptLike < _NH_IBD_DT_MAX)
    passesMultPre = num_PromptLikes_800us == 1
    passesMultPost = (dt_next_DelayedLike > _NH_MULT_POST_MIN
            or dt_next_DelayedLike == -1)
    return (passesDtCut and passesMultPre and passesMultPost)
