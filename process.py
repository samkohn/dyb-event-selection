'''
Process the basic data by computing various cuts, quantities and tags.

'''
from __future__ import print_function

from collections import deque

from ROOT import TFile, TTree
from flashers import fID, fPSD, isFlasher
from muons import isWSMuon, isADMuon, isShowerMuon
from translate import (TreeBuffer, float_value, assign_value,
        fetch_value, unsigned_int_value, long_value)

def main():

    filename = 'out.root'
    infile = TFile(filename, 'UPDATE')
    indata = infile.Get('data')
    outdata = TTree('computed', 'Computed quantities by Sam Kohn')

    buf = TreeBuffer()
    buf.fID = float_value()
    buf.fPSD = float_value()
    buf.tag_flasher = unsigned_int_value()
    buf.tag_WSMuon = unsigned_int_value()
    buf.tag_ADMuon = unsigned_int_value()
    buf.tag_ShowerMuon = unsigned_int_value()
    buf.dt_previous_WSMuon = long_value()
    buf.dt_previous_ADMuon = long_value()
    buf.dt_previous_ShowerMuon = long_value()
    buf.num_ShowerMuons_5sec = unsigned_int_value()
    buf.dts_ShowerMuons_5sec = long_value(20)

    outdata.Branch('fID', buf.fID, 'fID/F')
    outdata.Branch('fPSD', buf.fPSD, 'fPSD/F')
    outdata.Branch('tag_flasher', buf.tag_flasher, 'tag_flasher/i')
    outdata.Branch('tag_WSMuon', buf.tag_WSMuon, 'tag_WSMuon/i')
    outdata.Branch('tag_ADMuon', buf.tag_ADMuon, 'tag_ADMuon/i')
    outdata.Branch('tag_ShowerMuon', buf.tag_ShowerMuon,
            'tag_ShowerMuon/i')
    outdata.Branch('dt_previous_WSMuon', buf.dt_previous_WSMuon,
            'dt_previous_WSMuon/L')
    outdata.Branch('dt_previous_ADMuon', buf.dt_previous_ADMuon,
            'dt_previous_ADMuon/L')
    outdata.Branch('dt_previous_ShowerMuon', buf.dt_previous_ShowerMuon,
            'dt_previous_ShowerMuon/L')
    outdata.Branch('num_ShowerMuons_5sec', buf.num_ShowerMuons_5sec,
            'num_ShowerMuons_5sec/i')
    outdata.Branch('dts_ShowerMuons_5sec', buf.dts_ShowerMuons_5sec,
            'dts_ShowerMuons_5sec[num_ShowerMuons_5sec]/L')

    last_WSMuon_time = 0
    last_ADMuon_time = {n:0 for n in range(9)}
    last_ShowerMuon_time = {n:0 for n in range(9)}
    MUON_COUNT_TIME = 5*10**9  # 5 seconds, in nanoseconds
    recent_shower_muons = {n:deque() for n in range(9)}

    for event_number in xrange(indata.GetEntries()):
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)

        timestamp = fetch_value(indata, 'timeStamp', int)
        detector = fetch_value(indata, 'detector', int)
        nHit = fetch_value(indata, 'nHit', int)
        charge = fetch_value(indata, 'charge', float)
        fMax = fetch_value(indata, 'fMax', float)
        fQuad = fetch_value(indata, 'fQuad', float)
        fPSD_t1 = fetch_value(indata, 'fPSD_t1', float)
        fPSD_t2 = fetch_value(indata, 'fPSD_t2', float)

        try:
            event_fID = fID(fMax, fQuad)
            assign_value(buf.fID, event_fID)
            event_fPSD = fPSD(fPSD_t1, fPSD_t2)
            assign_value(buf.fPSD, event_fPSD)
            event_isFlasher = isFlasher(event_fID, event_fPSD)
            assign_value(buf.tag_flasher, event_isFlasher)
            event_isWSMuon = isWSMuon(detector, nHit)
            assign_value(buf.tag_WSMuon, event_isWSMuon)
            event_isADMuon = isADMuon(charge)
            assign_value(buf.tag_ADMuon, event_isADMuon)
            event_isShowerMuon = isShowerMuon(charge)
            assign_value(buf.tag_ShowerMuon, event_isShowerMuon)
        except:
            print(fMax, fQuad, fPSD_t1, fPSD_t2)
            raise

        if event_isWSMuon:
            last_WSMuon_time = timestamp
        if event_isADMuon:
            last_ADMuon_time[detector] = timestamp
        if event_isShowerMuon:
            last_ShowerMuon_time[detector] = timestamp
            recent_shower_muons[detector].append(timestamp)

        # Remove muons that happened greater than MUON_COUNT_TIME ago
        while len(recent_shower_muons[detector]) > 0 and (timestamp
                - recent_shower_muons[detector][0] > MUON_COUNT_TIME):
            recent_shower_muons[detector].popleft()

        recent_dts = [timestamp - muon_time for muon_time in
                recent_shower_muons[detector]]

        assign_value(buf.dt_previous_WSMuon, timestamp -
                last_WSMuon_time)
        assign_value(buf.dt_previous_ADMuon, timestamp -
                last_ADMuon_time[detector])
        assign_value(buf.dt_previous_ShowerMuon, timestamp -
                last_ShowerMuon_time[detector])
        assign_value(buf.num_ShowerMuons_5sec,
                len(recent_shower_muons[detector]))
        for i, dt in enumerate(recent_dts):
            assign_value(buf.dts_ShowerMuons_5sec, dt, i)

        outdata.Fill()

    infile.Write()
    infile.Close()

if __name__ == '__main__':
    main()
