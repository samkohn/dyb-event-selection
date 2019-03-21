'''
This module computes muon cuts.

'''

_WSMUON_NHIT_THRESHOLD = 12
_ADMUON_THRESHOLD = 3000
_SHOWER_MUON_THRESHOLD = 300000

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

