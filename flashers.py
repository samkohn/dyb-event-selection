'''
This module computes flasher cuts.

'''
from math import log10, pow, sqrt

from common import *

_NOT_RELEVANT = 3

def fID(fMax, fQuad):
    if fMax == 0 and fQuad == 0:
        return _NOT_RELEVANT
    else:
        return log10(pow(fQuad, 2) + pow(fMax/0.45, 2))

def fPSD(f_t1, f_t2):
    if f_t1 == 1 and f_t2 == 1:
        return _NOT_RELEVANT
    else:
        t1_term = 4*pow(1-f_t1, 2)
        t2_term = 1.8*pow(1-f_t2, 2)
        return log10(t1_term + t2_term)

def isFlasher(event_fID, event_fPSD, f2inch_maxQ, event_detector):
    if event_detector not in AD_DETECTORS:
        return False
    if event_fID == _NOT_RELEVANT and event_fPSD == _NOT_RELEVANT:
        return False
    fID_cut = int(event_fID >= 0)
    fPSD_cut = int(event_fPSD >= 0) * 2
    f2inch_cut = int(f2inch_maxQ > 100) * 4
    return fID_cut + fPSD_cut + f2inch_cut

def isFlasher_nH(event_fID, event_fPSD, f2inch_maxQ, event_detector):
    if event_detector not in AD_DETECTORS:
        return False
    if event_fID == _NOT_RELEVANT:
        return False
    fID_cut = int(event_fID >= 0)
    f2inch_cut = int(f2inch_maxQ > 100) * 4
    return fID_cut + f2inch_cut

RESID_FLASHER_TOP_RING_Z_MM = 2400
RESID_FLASHER_TOP_RING_R2_MM2 = 250000
RESID_FLASHER_QUN_R_MM = 2200
RESID_FLASHER_Q1Q2 = 0.6
RESID_FLASHER_FID = -0.3
RESID_FLASHER_TRIANGLE_Q1Q2 = 0.5
RESID_FLASHER_TRIANGLE_CONST = -0.8
def isFlasher_nH_no_resid(
    event_fID,
    event_fPSD,
    f2inch_maxQ,
    event_detector,
    event_Q1,
    event_Q2,
    event_x,
    event_y,
    event_z,
):
    is_nominal_flasher = isFlasher_nH(event_fID, event_fPSD, f2inch_maxQ, event_detector)
    event_r2 = event_x * event_x + event_y * event_y
    event_r = sqrt(event_r2)
    event_Q1_by_Q2 = event_Q1/event_Q2
    is_top_ring_flasher = (
        event_z > RESID_FLASHER_TOP_RING_Z_MM
        and event_r2 > RESID_FLASHER_TOP_RING_R2_MM2
    )
    is_qun_flasher = event_r > RESID_FLASHER_QUN_R_MM
    is_other_resid_flasher = (
        event_fID > RESID_FLASHER_FID
        and event_Q1_by_Q2 > RESID_FLASHER_Q1Q2
        and (
            event_Q1_by_Q2 * RESID_FLASHER_Q1Q2
            + RESID_FLASHER_TRIANGLE_CONST
            - event_fID
        ) < 0
    )
    return (
        is_nominal_flasher
        + is_top_ring_flasher * 8
        + is_qun_flasher * 16
        + is_other_resid_flasher * 32
    )
