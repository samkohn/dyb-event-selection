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
from delayeds import isDelayedLike, isIBDDelayed
from translate import (TreeBuffer, float_value, assign_value,
        fetch_value, int_value, unsigned_int_value, long_value, git_describe)

def done_with_cache(buf):
    '''
    Test whether an event should be removed from the event cache, based
    on whether certain values have been assigned to it.

    '''
    # An event is done if it is a flasher, or if not, then if the next
    # WSMuon and the next delayed-like events have been processed. (If
    # the event is not an AD event then the delayed-like event does not
    # have to have been processed.)
    found_next_WSMuon = buf.dt_next_WSMuon[0] != 0
    detector = buf.noTree_detector[0]
    found_next_DelayedLike = (detector not in AD_DETECTORS
            or buf.dt_next_DelayedLike[0] != 0)
    return found_next_WSMuon and found_next_DelayedLike

def create_computed_TTree(host_file):
    host_file.cd()
    git_description = git_describe()
    outdata = TTree('computed', 'Computed quantities by Sam Kohn (git: %s)' %
            git_description)
    # Initialize the "buffer" used to fill values into the TTree
    fill_buf = TreeBuffer()
    fill_buf.noTree_loopIndex = int_value()
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
    fill_buf.tag_IBDDelayed = unsigned_int_value()
    fill_buf.dt_previous_WSMuon = long_value()
    fill_buf.dt_next_WSMuon = long_value()
    fill_buf.dt_previous_ADMuon = long_value()
    fill_buf.dt_previous_ShowerMuon = long_value()
    fill_buf.dt_previous_PromptLike = long_value()
    fill_buf.energy_previous_PromptLike = float_value()
    fill_buf.dt_next_DelayedLike = long_value()
    fill_buf.num_ShowerMuons_5sec = unsigned_int_value()
    fill_buf.dts_ShowerMuons_5sec = long_value(20)
    fill_buf.num_PromptLikes_400us = unsigned_int_value()
    fill_buf.dts_PromptLikes_400us = long_value(200)

    # Initialize the new TTree so that each TBranch reads from the
    # appropriate TreeBuffer attribute
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
    outdata.Branch('tag_IBDDelayed', fill_buf.tag_IBDDelayed,
            'tag_IBDDelayed/i')
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
    outdata.Branch('energy_previous_PromptLike',
            fill_buf.energy_previous_PromptLike,
            'energy_previous_PromptLike/F')
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
    return outdata, fill_buf

class ProcessHelper(object):
    MUON_COUNT_TIME = 5*10**9  # 5 seconds, in nanoseconds
    PROMPT_COUNT_TIME = int(400e3)  # 400 us, in nanoseconds
    def __init__(self):
        # Initialize the machinery for dt_next_* and dt_previous_*
        # attributes. The event_cache stores the TreeBuffer for each event
        # until the next WSMuon and next DelayedLike have been found, thus
        # preserving the event order between the input and output TTrees.
        self.event_cache = deque()
        self.last_WSMuon_time = 0
        self.last_ADMuon_time = {n:0 for n in range(9)}
        self.last_ShowerMuon_time = {n:0 for n in range(9)}
        self.recent_shower_muons = {n:deque() for n in range(9)}
        self.last_PromptLike_time = {n:-1 for n in range(9)}
        self.last_PromptLike_energy = {n:-1 for n in range(9)}
        self.recent_promptlikes = {n:deque() for n in range(9)}

def main(entries, debug):

    filename = 'out.root'
    infile = TFile(filename, 'UPDATE')
    indata = infile.Get('data')
    outdata, fill_buf = create_computed_TTree(infile)

    helper = ProcessHelper()

    if entries == -1:
        entries = xrange(indata.GetEntries())
    else:
        entries = xrange(entries)
    if debug:
        entries = xrange(10000)
    def callback(event_cache):
        print(event_cache.noTree_loopIndex[0])
    for event_number in entries:
        # Load the current event in the input TTree
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)
        indata_list = fetch_indata(indata)

        one_iteration(event_number, indata_list, outdata, fill_buf, helper,
                callback)
    finish_emptying_cache(outdata, fill_buf, helper.event_cache, callback)
    outdata.Write()
    infile.Close()

def finish_emptying_cache(outdata, fill_buf, cache, callback):
    # After the event loop is finished, fill the remaining events from
    # the event_cache into the output TTree
    for cached_event in cache:
        if cached_event.dt_next_WSMuon[0] == 0:
            assign_value(cached_event.dt_next_WSMuon, -1)
            assign_value(cached_event.tag_WSMuonVeto, 0)
        if cached_event.dt_next_DelayedLike[0] == 0:
            assign_value(cached_event.dt_next_DelayedLike, -1)
        e = cached_event
        ibd_delayed = isIBDDelayed(
                e.tag_DelayedLike[0],
                e.dt_previous_PromptLike[0],
                e.num_PromptLikes_400us[0],
                e.dt_next_DelayedLike[0],
                e.tag_WSMuonVeto[0],
                e.tag_ADMuonVeto[0],
                e.tag_ShowerMuonVeto[0],
                e.tag_flasher[0])
        logging.debug('cached event: %s', cached_event)
        logging.debug('isIBDDelayed: %s', ibd_delayed)
        assign_value(cached_event.tag_IBDDelayed, ibd_delayed)
        callback(cached_event)
        cached_event.copyTo(fill_buf)
        outdata.Fill()

def fetch_indata(indata):
    # Fetch the necessary values from the input TTree
    timestamp = fetch_value(indata, 'timeStamp', int)
    triggerType = fetch_value(indata, 'triggerType', int)
    detector = fetch_value(indata, 'detector', int)
    site = fetch_value(indata, 'site', int)
    nHit = fetch_value(indata, 'nHit', int)
    charge = fetch_value(indata, 'charge', float)
    fMax = fetch_value(indata, 'fMax', float)
    fQuad = fetch_value(indata, 'fQuad', float)
    fPSD_t1 = fetch_value(indata, 'fPSD_t1', float)
    fPSD_t2 = fetch_value(indata, 'fPSD_t2', float)
    f2inch_maxQ = fetch_value(indata, 'f2inch_maxQ', float)
    energy = fetch_value(indata, 'energy', float)
    return (
            timestamp,
            triggerType,
            detector,
            site,
            nHit,
            charge,
            fMax,
            fQuad,
            fPSD_t1,
            fPSD_t2,
            f2inch_maxQ,
            energy
            )

def one_iteration(event_number, relevant_indata, outdata, fill_buf, helper,
        callback=lambda e:None):
    (timestamp, triggerType, detector, site, nHit, charge, fMax, fQuad,
            fPSD_t1, fPSD_t2, f2inch_maxQ, energy) = relevant_indata

    logging.debug('event cache size: %d', len(helper.event_cache))

    buf = fill_buf.clone_type()
    assign_value(buf.noTree_loopIndex, event_number)
    assign_value(buf.noTree_timestamp, timestamp)
    assign_value(buf.noTree_detector, detector)
    assign_value(buf.noTree_site, site)

    # Compute simple tags and values (those that only require data
    # from the current event)
    event_fID = fID(fMax, fQuad)
    assign_value(buf.fID, event_fID)
    event_fPSD = fPSD(fPSD_t1, fPSD_t2)
    assign_value(buf.fPSD, event_fPSD)
    event_isFlasher = isFlasher(event_fID, event_fPSD, f2inch_maxQ,
            detector)
    assign_value(buf.tag_flasher, event_isFlasher)
    event_isWSMuon = muons.isWSMuon(detector, nHit, triggerType)
    assign_value(buf.tag_WSMuon, event_isWSMuon)
    event_isADMuon = muons.isADMuon(detector, charge)
    assign_value(buf.tag_ADMuon, event_isADMuon)
    event_isShowerMuon = muons.isShowerMuon(detector, charge)
    assign_value(buf.tag_ShowerMuon, event_isShowerMuon)
    event_isPromptLike = isPromptLike(detector, energy)
    assign_value(buf.tag_PromptLike, event_isPromptLike)
    event_isDelayedLike = isDelayedLike(detector, energy)
    assign_value(buf.tag_DelayedLike, event_isDelayedLike)


    # Compute tags and values that count the number of previous
    # events within a given time satisfying some criteria.

    # Ensure the list of previous shower muons only has events more
    # recent than MUON_COUNT_TIME ago
    while len(helper.recent_shower_muons[detector]) > 0 and (timestamp
            - helper.recent_shower_muons[detector][0] > ProcessHelper.MUON_COUNT_TIME):
        helper.recent_shower_muons[detector].popleft()

    # Ensure the list of previous prompt-like events only has events
    # more recent than PROMPT_COUNT_TIME ago
    while len(helper.recent_promptlikes[detector]) > 0 and (timestamp -
            helper.recent_promptlikes[detector][0] > ProcessHelper.PROMPT_COUNT_TIME):
        helper.recent_promptlikes[detector].popleft()

    # Compute the dts for both recent muons and recent prompt-likes
    recent_showermuon_dts = [timestamp - muon_time for muon_time in
            helper.recent_shower_muons[detector]]
    recent_promptlike_dts = [timestamp - prompt_time for prompt_time
            in helper.recent_promptlikes[detector]]

    # Assign the values for dt_previous_* and for the lists of
    # recent muons and prompt-likes.
    assign_value(buf.dt_previous_WSMuon, timestamp -
            helper.last_WSMuon_time)
    assign_value(buf.dt_previous_ADMuon, timestamp -
            helper.last_ADMuon_time[detector])
    assign_value(buf.dt_previous_ShowerMuon, timestamp -
            helper.last_ShowerMuon_time[detector])
    assign_value(buf.dt_previous_PromptLike, timestamp -
            helper.last_PromptLike_time[detector])
    assign_value(buf.energy_previous_PromptLike,
            helper.last_PromptLike_energy[detector])
    assign_value(buf.num_ShowerMuons_5sec,
            len(helper.recent_shower_muons[detector]))
    assign_value(buf.num_PromptLikes_400us,
            len(helper.recent_promptlikes[detector]))
    for i, dt in enumerate(recent_showermuon_dts):
        assign_value(buf.dts_ShowerMuons_5sec, dt, i)
    for i, dt in enumerate(recent_promptlike_dts):
        assign_value(buf.dts_PromptLikes_400us, dt, i)

    # Compute muon vetoes that only rely on previous events
    assign_value(buf.tag_ADMuonVeto,
            muons.isVetoedByADMuon(buf.dt_previous_ADMuon[0]))
    assign_value(buf.tag_ShowerMuonVeto,
            muons.isVetoedByShowerMuon(buf.dt_previous_ShowerMuon[0]))

    # Update the dt_previous_* and dt_next_* values

    # This comes after values are assigned to the buffer because we
    # don't want dt_previous_* to be 0.
    if event_isWSMuon:
        logging.debug("isWSMuon")
        helper.last_WSMuon_time = timestamp
        # Assign dt_next_WSMuon to the events in the event_cache
        for cached_event in helper.event_cache:
            # Some events might already have been assigned, so skip
            # those
            if cached_event.dt_next_WSMuon[0] != 0:
                continue
            # WSMuons are common across the entire EH/site
            if cached_event.noTree_site[0] == site:
                logging.debug('taggedNextWS%d', cached_event.noTree_timestamp[0])
                assign_value(cached_event.dt_next_WSMuon,
                        timestamp - cached_event.noTree_timestamp[0])
                assign_value(cached_event.tag_WSMuonVeto,
                        muons.isVetoedByWSMuon(cached_event.dt_previous_WSMuon[0],
                            cached_event.dt_next_WSMuon[0]))

    if event_isADMuon:
        helper.last_ADMuon_time[detector] = timestamp
    if event_isShowerMuon:
        helper.last_ShowerMuon_time[detector] = timestamp
        helper.recent_shower_muons[detector].append(timestamp)
    # Don't include prompt-like flasher events
    if event_isPromptLike and not event_isFlasher:
        helper.last_PromptLike_time[detector] = timestamp
        helper.last_PromptLike_energy[detector] = energy
        helper.recent_promptlikes[detector].append(timestamp)
    # Don't include prompt-like delayed events
    if event_isDelayedLike and not event_isFlasher:
        logging.debug("isDelayedLike")
        # Assign dt_next_PromptLike to events in the event_cache
        for cached_event in helper.event_cache:
            # Some events might already have been assigned, so skip
            # those
            if cached_event.dt_next_DelayedLike[0] != 0:
                continue
            # PromptLikes are restricted to a single AD
            if cached_event.noTree_detector[0] == detector:
                logging.debug('taggedNextDelayed%d',
                        cached_event.noTree_timestamp[0])
                assign_value(cached_event.dt_next_DelayedLike,
                        timestamp
                        - cached_event.noTree_timestamp[0])

    # Determine which of the oldest events are ready to go into the
    # new TTree. It is possible (due to different ADs) that the
    # oldest event might not be ready but another event might be. In
    # that case, to preserve event order, we will wait for the
    # oldest event to be ready anyways.
    num_to_delete = 0
    all_done_with_cache = True
    cache_size = len(helper.event_cache)
    while all_done_with_cache and num_to_delete < cache_size:
        cached_event = helper.event_cache[num_to_delete]
        if done_with_cache(cached_event):
            logging.debug('event is done with cache')
            num_to_delete += 1
        else:
            logging.debug('event is not done with cache')
            all_done_with_cache = False
    logging.debug('deleting %d events', num_to_delete)
    # Remove the oldest events from the cache and fill them into the
    # new TTree, computing the final IBD tags while we're at it.
    for _ in range(num_to_delete):
        cached_event = helper.event_cache.popleft()
        e = cached_event
        ibd_delayed = isIBDDelayed(
                e.tag_DelayedLike[0],
                e.dt_previous_PromptLike[0],
                e.num_PromptLikes_400us[0],
                e.dt_next_DelayedLike[0],
                e.tag_WSMuonVeto[0],
                e.tag_ADMuonVeto[0],
                e.tag_ShowerMuonVeto[0],
                e.tag_flasher[0])
        logging.debug('cached event: %s', cached_event)
        logging.debug('isIBDDelayed: %s', ibd_delayed)
        assign_value(cached_event.tag_IBDDelayed, ibd_delayed)
        callback(cached_event)
        cached_event.copyTo(fill_buf)
        outdata.Fill()

    helper.event_cache.append(buf)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--events', type=int, default=-1)
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.events, args.debug)
