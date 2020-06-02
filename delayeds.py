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
        dt_previous_PromptLike,
        dr_previous_PromptLike,
        muonVeto_previous_PromptLike,
        num_PromptLikes_400us,
        dts_recent_PromptLikes,
        dt_next_ADevent,
        dt_next_WSMuon,
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

_NH_DMC_EMIN = 1.9
_NH_DMC_EMAX = 2.74
_NH_DMC_IBD_DT_MIN = int(1e3)
_NH_DMC_IBD_DT_MAX = int(400e3)
_NH_DMC_MULT_POST_MIN = int(400e3)

def isDelayedLike_nH_DMC(event_detector, event_energy):
    return (event_detector in AD_DETECTORS
            and event_energy > _NH_DMC_EMIN
            and event_energy < _NH_DMC_EMAX)

def isIBDDelayed_nH_DMC(tag_DelayedLike,
        dt_previous_PromptLike,
        dr_previous_PromptLike,
        muonVeto_previous_PromptLike,
        num_PromptLikes_800us,
        dts_PromptLikes_800us,
        dt_next_ADevent,
        dt_next_WSMuon,
        dt_next_DelayedLike, tag_WSMuonVeto, tag_ADMuonVeto,
        tag_ShowerMuonVeto, tag_flasher):
    if not tag_DelayedLike:
        return False
    if tag_flasher:
        return False
    if tag_WSMuonVeto or tag_ADMuonVeto or tag_ShowerMuonVeto:
        return False
    passesDtCut = (dt_previous_PromptLike > _NH_DMC_IBD_DT_MIN
            and dt_previous_PromptLike < _NH_DMC_IBD_DT_MAX)
    passesMultPre = num_PromptLikes_800us == 1
    passesMultPost = (dt_next_DelayedLike > _NH_DMC_MULT_POST_MIN
            or dt_next_DelayedLike == -1)
    return (passesDtCut and passesMultPre and passesMultPost)

_NH_THU_MULT_PRE_MIN = int(400e3)
_NH_THU_MULT_POST_MIN = int(400e3)
_NH_THU_DR_MAX = 500
_NH_THU_MAX_TIME = int(1500e3)
_NH_THU_MIN_TIME = int(1e3)
_NH_THU_DIST_TIME_MAX = 800  # Known as phi or DT cut
_NH_THU_DIST_TIME_CONST = 1000/600e3  # Constant for converting time to distance
_NH_THU_DIST_TIME_STR = (f'(dr_to_prompt[1] + {_NH_THU_DIST_TIME_CONST} *'
    'dt_to_prompt[1])')
_NH_THU_DIST_TIME_CUT_STR = f'({_NH_THU_DIST_TIME_STR} < {_NH_THU_DIST_TIME_MAX})'

def nH_THU_DT(dr_mm, dt_ns):
    return dr_mm + dt_ns * _NH_THU_DIST_TIME_CONST


def isDelayedLike_nH_THU(event_detector, event_energy):
    return isDelayedLike_nH_DMC(event_detector, event_energy)

def isIBDDelayed_nH_THU(tag_DelayedLike,
        dt_previous_PromptLike,
        dr_previous_PromptLike,
        muonVeto_previous_PromptLike,
        num_ADevents_800us,
        dts_ADevents_800us,
        dt_next_ADevent,
        dt_next_WSMuon,
        dt_next_DelayedLike, tag_WSMuonVeto, tag_ADMuonVeto,
        tag_ShowerMuonVeto, tag_flasher):
    if not tag_DelayedLike:
        return False
    if tag_flasher not in (0, 2):
        return False
    if (tag_WSMuonVeto
            or tag_ADMuonVeto
            or tag_ShowerMuonVeto
            or muonVeto_previous_PromptLike):
        return False
    passesDtCut = (dt_previous_PromptLike > _NH_THU_DT_MIN
            and dt_previous_PromptLike < _NH_THU_DT_MAX)
    passesDrCut = dr_previous_PromptLike < _NH_THU_DR_MAX
    passesMultPre = (num_ADevents_800us > 0
            and all(dt == dt_previous_PromptLike
                or dt > dt_previous_PromptLike + _NH_THU_MULT_PRE_MIN
                for dt in dts_ADevents_800us[:num_ADevents_800us]))
    passesMultPost = (((dt_next_ADevent + dt_previous_PromptLike
            > _NH_THU_MULT_POST_MIN) or dt_next_ADevent == -1)
        and ((dt_next_WSMuon + dt_previous_PromptLike
            > _NH_THU_MULT_POST_MIN) or dt_next_WSMuon == -1))
    return passesDtCut and passesDrCut and passesMultPre and passesMultPost

