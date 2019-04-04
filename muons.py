'''
This module computes muon cuts.

'''

_WSMUON_NHIT_THRESHOLD = 12
_ADMUON_THRESHOLD = 3000
_SHOWER_MUON_THRESHOLD = 300000

_WSMUON_VETO_NEXT_NS = 2000
_WSMUON_VETO_LAST_NS = 600000
_ADMUON_VETO_LAST_NS = 1400000
_SHOWER_MUON_VETO_LAST_NS = 400000000

def isWSMuon(event_detector, event_nHit):
    IWS = 5
    OWS = 6
    return (event_detector in (IWS, OWS) and
            event_nHit >= _WSMUON_NHIT_THRESHOLD)

def isADMuon(event_charge):
    return (event_charge > _ADMUON_THRESHOLD and
            event_charge <= _SHOWER_MUON_THRESHOLD)

def isShowerMuon(event_charge):
    return event_charge > _SHOWER_MUON_THRESHOLD

def isVetoedByWSMuon(dt_last_ws, dt_next_ws):
    return (dt_last_ws > _WSMUON_VETO_LAST_NS
            and dt_next_ws > _WSMUON_VETO_NEXT_NS)

def isVetoedByADMuon(dt_last_ad):
    return (dt_last_ad > _ADMUON_VETO_LAST_NS)

def isVetoedByShowerMuon(dt_last_shower):
    return (dt_last_shower > _SHOWER_MUON_VETO_LAST_NS)

