'''
Process the basic data by computing various cuts, quantities and tags.

'''
from __future__ import print_function

from collections import deque

from ROOT import TFile, TTree
from flashers import fID, fPSD, isFlasher
import muons
from translate import (TreeBuffer, float_value, assign_value,
        fetch_value, int_value, unsigned_int_value, long_value)

def done_with_cache(buf):
    if buf.dt_next_WSMuon != 0:
        return True
    return False

def main():

    filename = 'out.root'
    infile = TFile(filename, 'UPDATE')
    indata = infile.Get('data')
    outdata = TTree('computed', 'Computed quantities by Sam Kohn')

    fill_buf = TreeBuffer()
    fill_buf.noTree_timestamp = long_value()
    fill_buf.noTree_detector = int_value()
    fill_buf.noTree_site = int_value()
    fill_buf.fID = float_value()
    fill_buf.fPSD = float_value()
    fill_buf.tag_flasher = unsigned_int_value()
    fill_buf.tag_WSMuon = unsigned_int_value()
    fill_buf.tag_ADMuon = unsigned_int_value()
    fill_buf.tag_ShowerMuon = unsigned_int_value()
    fill_buf.tag_WSMuonVeto = unsigned_int_value()
    fill_buf.tag_ADMuonVeto = unsigned_int_value()
    fill_buf.tag_ShowerMuonVeto = unsigned_int_value()
    fill_buf.dt_previous_WSMuon = long_value()
    fill_buf.dt_next_WSMuon = long_value()
    fill_buf.dt_previous_ADMuon = long_value()
    fill_buf.dt_previous_ShowerMuon = long_value()
    fill_buf.num_ShowerMuons_5sec = unsigned_int_value()
    fill_buf.dts_ShowerMuons_5sec = long_value(20)

    outdata.Branch('fID', fill_buf.fID, 'fID/F')
    outdata.Branch('fPSD', fill_buf.fPSD, 'fPSD/F')
    outdata.Branch('tag_flasher', fill_buf.tag_flasher, 'tag_flasher/i')
    outdata.Branch('tag_WSMuon', fill_buf.tag_WSMuon, 'tag_WSMuon/i')
    outdata.Branch('tag_ADMuon', fill_buf.tag_ADMuon, 'tag_ADMuon/i')
    outdata.Branch('tag_ShowerMuon', fill_buf.tag_ShowerMuon,
            'tag_ShowerMuon/i')
    outdata.Branch('tag_WSMuonVeto', fill_buf.tag_WSMuonVeto,
            'tag_WSMuonVeto/i')
    outdata.Branch('tag_ADMuonVeto', fill_buf.tag_ADMuonVeto,
            'tag_ADMuonVeto/i')
    outdata.Branch('tag_ShowerMuonVeto', fill_buf.tag_ShowerMuonVeto,
            'tag_ShowerMuonVeto/i')
    outdata.Branch('dt_previous_WSMuon', fill_buf.dt_previous_WSMuon,
            'dt_previous_WSMuon/L')
    outdata.Branch('dt_next_WSMuon', fill_buf.dt_next_WSMuon,
            'dt_next_WSMuon/L')
    outdata.Branch('dt_previous_ADMuon', fill_buf.dt_previous_ADMuon,
            'dt_previous_ADMuon/L')
    outdata.Branch('dt_previous_ShowerMuon', fill_buf.dt_previous_ShowerMuon,
            'dt_previous_ShowerMuon/L')
    outdata.Branch('num_ShowerMuons_5sec', fill_buf.num_ShowerMuons_5sec,
            'num_ShowerMuons_5sec/i')
    outdata.Branch('dts_ShowerMuons_5sec', fill_buf.dts_ShowerMuons_5sec,
            'dts_ShowerMuons_5sec[num_ShowerMuons_5sec]/L')

    event_cache = deque()
    last_WSMuon_time = 0
    last_ADMuon_time = {n:0 for n in range(9)}
    last_ShowerMuon_time = {n:0 for n in range(9)}
    MUON_COUNT_TIME = 5*10**9  # 5 seconds, in nanoseconds
    recent_shower_muons = {n:deque() for n in range(9)}

    for event_number in xrange(indata.GetEntries()):
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)
        buf = fill_buf.clone_type()

        timestamp = fetch_value(indata, 'timeStamp', int)
        assign_value(buf.noTree_timestamp, timestamp)
        detector = fetch_value(indata, 'detector', int)
        assign_value(buf.noTree_detector, detector)
        site = fetch_value(indata, 'site', int)
        assign_value(buf.noTree_site, site)
        nHit = fetch_value(indata, 'nHit', int)
        charge = fetch_value(indata, 'charge', float)
        fMax = fetch_value(indata, 'fMax', float)
        fQuad = fetch_value(indata, 'fQuad', float)
        fPSD_t1 = fetch_value(indata, 'fPSD_t1', float)
        fPSD_t2 = fetch_value(indata, 'fPSD_t2', float)
        f2inch_maxQ = fetch_value(indata, 'f2inch_maxQ', float)

        try:
            event_fID = fID(fMax, fQuad)
            assign_value(buf.fID, event_fID)
            event_fPSD = fPSD(fPSD_t1, fPSD_t2)
            assign_value(buf.fPSD, event_fPSD)
            event_isFlasher = isFlasher(event_fID, event_fPSD, f2inch_maxQ)
            assign_value(buf.tag_flasher, event_isFlasher)
            event_isWSMuon = muons.isWSMuon(detector, nHit)
            assign_value(buf.tag_WSMuon, event_isWSMuon)
            event_isADMuon = muons.isADMuon(charge)
            assign_value(buf.tag_ADMuon, event_isADMuon)
            event_isShowerMuon = muons.isShowerMuon(charge)
            assign_value(buf.tag_ShowerMuon, event_isShowerMuon)
        except:
            print(fMax, fQuad, fPSD_t1, fPSD_t2)
            raise

        if event_isWSMuon:
            last_WSMuon_time = timestamp
            for cached_event in event_cache:
                if cached_event.noTree_site[0] == site:
                    assign_value(cached_event.dt_next_WSMuon,
                            timestamp - cached_event.noTree_timestamp[0])
                    assign_value(cached_event.tag_WSMuonVeto,
                            muons.isVetoedByWSMuon(cached_event.dt_previous_WSMuon[0],
                                cached_event.dt_next_WSMuon[0]))
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

        # Compute muon vetoes
        assign_value(buf.tag_ADMuonVeto,
                muons.isVetoedByADMuon(buf.dt_previous_ADMuon[0]))
        assign_value(buf.tag_ShowerMuonVeto,
                muons.isVetoedByShowerMuon(buf.dt_previous_ShowerMuon[0]))

        # Determine which of the oldest events are ready to go into the
        # new TTree. It is possible (due to different ADs) that the
        # oldest event might not be ready but another event might be. In
        # that case, to preserve event order, we will wait for the
        # oldest event to be ready anyways.
        num_to_delete = 0
        all_done_with_cache = True
        while all_done_with_cache:
            if done_with_cache(cached_event):
                num_to_delete += 1
            else:
                all_done_with_cache = False
        # Remove the oldest events from the cache and fill them into the
        # new TTree
        for _ in range(num_to_delete):
            cached_event = event_cache.popleft()
            cached_event.copyTo(fill_buf)
            outdata.Fill()

        event_cache.append(buf)
        #outdata.Fill()
    for cached_event in event_cache:
        assign_value(cached_event.dt_next_WSMuon, -1)
        assign_value(cached_event.tag_WSMuonVeto, 2)
        cached_event.copyTo(fill_buf)
        outdata.Fill()
    event_cache = []

    infile.Write()
    infile.Close()

if __name__ == '__main__':
    main()
