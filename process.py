'''
Process the basic data by computing various cuts, quantities and tags.

'''
from __future__ import print_function

from collections import deque
import time
import argparse
import logging

from ROOT import TFile, TTree
from common import *
from flashers import fID, fPSD, isFlasher
import muons
from prompts import isPromptLike
from delayeds import isDelayedLike
from translate import (TreeBuffer, float_value, assign_value,
        fetch_value, int_value, unsigned_int_value, long_value)

def done_with_cache(buf):
    flasher_done = buf.tag_flasher[0] != 0
    found_next_WSMuon = buf.dt_next_WSMuon[0] != 0
    detector = buf.noTree_detector[0]
    found_next_DelayedLike = (detector not in AD_DETECTORS
            or buf.dt_next_DelayedLike[0] != 0)
    return (flasher_done
            or (found_next_WSMuon
                and found_next_DelayedLike))

def main(debug):

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
    fill_buf.tag_PromptLike = unsigned_int_value()
    fill_buf.tag_DelayedLike = unsigned_int_value()
    fill_buf.dt_previous_WSMuon = long_value()
    fill_buf.dt_next_WSMuon = long_value()
    fill_buf.dt_previous_ADMuon = long_value()
    fill_buf.dt_previous_ShowerMuon = long_value()
    fill_buf.dt_previous_PromptLike = long_value()
    fill_buf.dt_next_DelayedLike = long_value()
    fill_buf.num_ShowerMuons_5sec = unsigned_int_value()
    fill_buf.dts_ShowerMuons_5sec = long_value(20)
    fill_buf.num_PromptLikes_400us = unsigned_int_value()
    fill_buf.dts_PromptLikes_400us = long_value(200)

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
    outdata.Branch('tag_PromptLike', fill_buf.tag_PromptLike,
            'tag_PromptLike/i')
    outdata.Branch('tag_DelayedLike', fill_buf.tag_DelayedLike,
            'tag_DelayedLike/i')
    outdata.Branch('dt_previous_WSMuon', fill_buf.dt_previous_WSMuon,
            'dt_previous_WSMuon/L')
    outdata.Branch('dt_next_WSMuon', fill_buf.dt_next_WSMuon,
            'dt_next_WSMuon/L')
    outdata.Branch('dt_previous_ADMuon', fill_buf.dt_previous_ADMuon,
            'dt_previous_ADMuon/L')
    outdata.Branch('dt_previous_ShowerMuon', fill_buf.dt_previous_ShowerMuon,
            'dt_previous_ShowerMuon/L')
    outdata.Branch('dt_previous_PromptLike',
            fill_buf.dt_previous_PromptLike, 'dt_previous_PromptLike/L')
    outdata.Branch('dt_next_DelayedLike', fill_buf.dt_next_DelayedLike,
            'dt_next_DelayedLike/L')
    outdata.Branch('num_ShowerMuons_5sec', fill_buf.num_ShowerMuons_5sec,
            'num_ShowerMuons_5sec/i')
    outdata.Branch('dts_ShowerMuons_5sec', fill_buf.dts_ShowerMuons_5sec,
            'dts_ShowerMuons_5sec[num_ShowerMuons_5sec]/L')
    outdata.Branch('num_PromptLikes_400us',
            fill_buf.num_PromptLikes_400us, 'num_PromptLikes_400us/i')
    outdata.Branch('dts_PromptLikes_400us',
            fill_buf.dts_PromptLikes_400us,
            'dts_PromptLikes_400us[num_PromptLikes_400us]/L')

    event_cache = deque()
    cache_assigned_WS = 0
    cache_assigned_delayed = {n: 0 for n in range(9)}
    last_WSMuon_time = 0
    last_ADMuon_time = {n:0 for n in range(9)}
    last_ShowerMuon_time = {n:0 for n in range(9)}
    MUON_COUNT_TIME = 5*10**9  # 5 seconds, in nanoseconds
    recent_shower_muons = {n:deque() for n in range(9)}
    last_PromptLike_time = {n:0 for n in range(9)}
    PROMPT_COUNT_TIME = int(400e3)  # 400 us, in nanoseconds
    recent_promptlikes = {n:deque() for n in range(9)}

    for event_number in xrange(indata.GetEntries()):
        logging.debug('event cache size: %d', len(event_cache))
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
        energy = fetch_value(indata, 'energy', float)

        event_fID = fID(fMax, fQuad)
        assign_value(buf.fID, event_fID)
        event_fPSD = fPSD(fPSD_t1, fPSD_t2)
        assign_value(buf.fPSD, event_fPSD)
        event_isFlasher = isFlasher(event_fID, event_fPSD, f2inch_maxQ)
        assign_value(buf.tag_flasher, event_isFlasher)
        if event_isFlasher:
            event_isWSMuon = False
            event_isADMuon = False
            event_isShowerMuon = False
        else:
            event_isWSMuon = muons.isWSMuon(detector, nHit)
            assign_value(buf.tag_WSMuon, event_isWSMuon)
            event_isADMuon = muons.isADMuon(detector, charge)
            assign_value(buf.tag_ADMuon, event_isADMuon)
            event_isShowerMuon = muons.isShowerMuon(detector, charge)
            assign_value(buf.tag_ShowerMuon, event_isShowerMuon)
        event_isPromptLike = isPromptLike(detector, energy)
        assign_value(buf.tag_PromptLike, event_isPromptLike)
        event_isDelayedLike = isDelayedLike(detector, energy)
        assign_value(buf.tag_DelayedLike, event_isDelayedLike)


        # Remove muons that happened greater than MUON_COUNT_TIME ago
        while len(recent_shower_muons[detector]) > 0 and (timestamp
                - recent_shower_muons[detector][0] > MUON_COUNT_TIME):
            recent_shower_muons[detector].popleft()

        # Remove prompt-likes that happened greater than PROMPT_COUNT_TIME ago
        while len(recent_promptlikes[detector]) > 0 and (timestamp -
                recent_promptlikes[detector][0] > PROMPT_COUNT_TIME):
            recent_promptlikes[detector].popleft()

        recent_showermuon_dts = [timestamp - muon_time for muon_time in
                recent_shower_muons[detector]]
        recent_promptlike_dts = [timestamp - prompt_time for prompt_time
                in recent_promptlikes[detector]]

        assign_value(buf.dt_previous_WSMuon, timestamp -
                last_WSMuon_time)
        assign_value(buf.dt_previous_ADMuon, timestamp -
                last_ADMuon_time[detector])
        assign_value(buf.dt_previous_ShowerMuon, timestamp -
                last_ShowerMuon_time[detector])
        assign_value(buf.dt_previous_PromptLike, timestamp -
                last_PromptLike_time[detector])
        assign_value(buf.num_ShowerMuons_5sec,
                len(recent_shower_muons[detector]))
        assign_value(buf.num_PromptLikes_400us,
                len(recent_promptlikes[detector]))
        for i, dt in enumerate(recent_showermuon_dts):
            assign_value(buf.dts_ShowerMuons_5sec, dt, i)
        for i, dt in enumerate(recent_promptlike_dts):
            assign_value(buf.dts_PromptLikes_400us, dt, i)

        # Compute muon vetoes
        assign_value(buf.tag_ADMuonVeto,
                muons.isVetoedByADMuon(buf.dt_previous_ADMuon[0]))
        assign_value(buf.tag_ShowerMuonVeto,
                muons.isVetoedByShowerMuon(buf.dt_previous_ShowerMuon[0]))

        # Update the dt_previous_* and dt_next_* values
        #
        # Ths comes after values are assigned to the buffer because we
        # don't want dt_previous_* to be 0.
        if event_isWSMuon:
            logging.debug("isWSMuon")
            last_WSMuon_time = timestamp
            for cached_event in event_cache[cache_assigned_WS:]:
                if cached_event.noTree_site[0] == site:
                    assign_value(cached_event.dt_next_WSMuon,
                            timestamp - cached_event.noTree_timestamp[0])
                    assign_value(cached_event.tag_WSMuonVeto,
                            muons.isVetoedByWSMuon(cached_event.dt_previous_WSMuon[0],
                                cached_event.dt_next_WSMuon[0]))
            cache_assigned_WS = len(event_cache)
        if event_isADMuon:
            last_ADMuon_time[detector] = timestamp
        if event_isShowerMuon:
            last_ShowerMuon_time[detector] = timestamp
            recent_shower_muons[detector].append(timestamp)
        if event_isPromptLike and not event_isFlasher:
            last_PromptLike_time[detector] = timestamp
            recent_promptlikes[detector].append(timestamp)
        if event_isDelayedLike and not event_isFlasher:
            start_point = cache_assigned_delayed[detector]
            for cached_event in event_cache[start_point:]:
                if cached_event.noTree_detector[0] == detector:
                    assign_value(cached_event.dt_next_DelayedLike,
                            timestamp
                            - cached_event.noTree_timestamp[0])
            cache_assigned_delayed[detector] = len(event_cache)

        # Determine which of the oldest events are ready to go into the
        # new TTree. It is possible (due to different ADs) that the
        # oldest event might not be ready but another event might be. In
        # that case, to preserve event order, we will wait for the
        # oldest event to be ready anyways.
        num_to_delete = 0
        all_done_with_cache = True
        cache_size = len(event_cache)
        while all_done_with_cache and num_to_delete < cache_size:
            cached_event = event_cache[num_to_delete]
            if done_with_cache(cached_event):
                logging.debug('event is done with cache')
                num_to_delete += 1
            else:
                logging.debug('event is not done with cache')
                all_done_with_cache = False
        logging.debug('first %d events have next_WS', cache_assigned_WS)
        logging.debug('first events that have next_delayed, by detector:')
        logging.debug(cache_assigned_delayed)
        logging.debug('deleting %d events', num_to_delete)
        cache_assigned_WS -= num_to_delete
        for key in cache_assigned_delayed:
            if cache_assigned_delayed[key] != 0:
                cache_assigned_delayed[key] -= num_to_delete
        # Remove the oldest events from the cache and fill them into the
        # new TTree
        for _ in range(num_to_delete):
            cached_event = event_cache.popleft()
            cached_event.copyTo(fill_buf)
            outdata.Fill()

        event_cache.append(buf)
        if debug:
            time.sleep(1)
        #outdata.Fill()
    for cached_event in event_cache:
        assign_value(cached_event.dt_next_WSMuon, -1)
        assign_value(cached_event.tag_WSMuonVeto, 2)
        assign_value(cached_event.dt_next_DelayedLike, -1)
        cached_event.copyTo(fill_buf)
        outdata.Fill()
    event_cache = []

    infile.Write()
    infile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', action='store_true')
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.debug)
