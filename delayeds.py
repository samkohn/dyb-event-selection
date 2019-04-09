'''
This module finds delayed-like events.

'''
from common import *

_EMIN = 6
_EMAX = 12
_IBD_DT_MIN = int(1e3)
_IBD_DT_MAX = int(200e3)
_MULT_POST_MIN = int(200e3)

def isDelayedLike(event_detector, event_energy):
    return (event_detector in AD_DETECTORS
            and event_energy > _EMIN
            and event_energy < _EMAX)

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
    passesMultPost = dt_next_DelayedLike > _MULT_POST_MIN
    return (passesDtCut and passesMultPre and passesMultPost)

