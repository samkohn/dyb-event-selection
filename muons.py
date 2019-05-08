'''
This module computes muon cuts.

'''
from common import *

_WSMUON_NHIT_THRESHOLD = 12
_ADMUON_THRESHOLD = 3000
_SHOWER_MUON_THRESHOLD = int(3e5)

_WSMUON_VETO_NEXT_NS = int(2e3)
_WSMUON_VETO_LAST_NS = int(600e3)
_ADMUON_VETO_LAST_NS = int(1.4e6)
_SHOWER_MUON_VETO_LAST_NS = int(0.4e9)

def isWSMuon(event_detector, event_nHit, event_triggerType):
    return (event_detector in WP_DETECTORS and
            event_nHit > _WSMUON_NHIT_THRESHOLD and
            (event_triggerType & 0x2) == 0)

def isADMuon(event_detector, event_charge):
    return (event_detector in AD_DETECTORS
            and event_charge > _ADMUON_THRESHOLD
            and event_charge <= _SHOWER_MUON_THRESHOLD)

def isShowerMuon(event_detector, event_charge):
    return (event_detector in AD_DETECTORS
            and event_charge > _SHOWER_MUON_THRESHOLD)

def isVetoedByWSMuon(dt_last_ws, dt_next_ws):
    return (dt_last_ws < _WSMUON_VETO_LAST_NS
            or dt_next_ws < _WSMUON_VETO_NEXT_NS)

def isVetoedByADMuon(dt_last_ad):
    return (dt_last_ad < _ADMUON_VETO_LAST_NS)

def isVetoedByShowerMuon(dt_last_shower):
    return (dt_last_shower < _SHOWER_MUON_VETO_LAST_NS)

