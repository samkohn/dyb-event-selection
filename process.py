'''
Process the basic data by computing various cuts, quantities and tags.

'''
from __future__ import print_function

from collections import deque
import time
import argparse
import logging
import math

from ROOT import TFile, TTree
from common import *
from flashers import fID, fPSD
import flashers
import muons
import prompts
import delayeds
from adevent import isADEvent_THU
from translate import (TreeBuffer, float_value, assign_value,
        fetch_value, int_value, unsigned_int_value, long_value, git_describe)

def done_with_cache(buf, selection_name):
    '''
    Test whether an event should be removed from the event cache, based
    on whether certain values have been assigned to it.

    '''
    # An event is done if it is a flasher, or if not, then if the next
    # WSMuon and the next delayed-like events have been processed. (If
    # the event is not an AD event then the delayed-like event does not
    # have to have been processed.)
    found_next_WSMuon = buf.dt_next_WSMuon[0] != -1
    detector = buf.detector[0]
    found_next_DelayedLike = (detector not in AD_DETECTORS
            or buf.dt_next_DelayedLike[0] != -1)
    found_next_ADevent = (buf.dt_next_ADevent[0] != -1
            or detector not in AD_DETECTORS
            or selection_name != 'nh_THU')
    return found_next_WSMuon and found_next_DelayedLike and found_next_ADevent

def create_computed_TTree(name, host_file, selection_name, title=None):
    git_description = git_describe()
    if title is None:
        title = ('Computed quantities by Sam Kohn (git: %s)' %
            git_description)
    else:
        try:
            title = title % git_description
        except:
            pass
    host_file.cd()
    outdata = TTree(name, title)
    # Initialize the "buffer" used to fill values into the TTree
    fill_buf = TreeBuffer()
    fill_buf.noTree_loopIndex = int_value()
    fill_buf.timestamp = long_value()
    fill_buf.timestamp_seconds = int_value()
    fill_buf.timestamp_nanoseconds = int_value()
    fill_buf.detector = int_value()
    fill_buf.site = int_value()
    fill_buf.run = unsigned_int_value()
    fill_buf.fileno = unsigned_int_value()
    fill_buf.triggerNumber = int_value()
    fill_buf.detector = int_value()
    fill_buf.site = int_value()
    fill_buf.triggerType = unsigned_int_value()
    fill_buf.nHit = int_value()
    fill_buf.charge = float_value()
    fill_buf.fQuad = float_value()
    fill_buf.fMax = float_value()
    fill_buf.fPSD_t1 = float_value()
    fill_buf.fPSD_t2 = float_value()
    fill_buf.f2inch_maxQ = float_value()
    fill_buf.energy = float_value()
    fill_buf.x = float_value()
    fill_buf.y = float_value()
    fill_buf.z = float_value()
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
    fill_buf.nHit_previous_WSMuon = int_value()
    fill_buf.dt_next_WSMuon = long_value()
    fill_buf.nHit_next_WSMuon = int_value()
    fill_buf.dt_previous_ADMuon = long_value()
    fill_buf.charge_previous_ADMuon = float_value()
    fill_buf.dt_previous_ShowerMuon = long_value()
    fill_buf.charge_previous_ShowerMuon = float_value()
    fill_buf.dt_previous_PromptLike = long_value()
    fill_buf.energy_previous_PromptLike = float_value()
    fill_buf.muonVeto_previous_PromptLike = unsigned_int_value()
    fill_buf.x_previous_PromptLike = float_value()
    fill_buf.y_previous_PromptLike = float_value()
    fill_buf.z_previous_PromptLike = float_value()
    fill_buf.dr_previous_PromptLike = float_value()
    fill_buf.dt_next_ADevent = long_value()
    fill_buf.dt_next_DelayedLike = long_value()
    fill_buf.num_ShowerMuons_5sec = unsigned_int_value()
    fill_buf.dts_ShowerMuons_5sec = long_value(20)
    fill_buf.num_recent_PromptLikes = unsigned_int_value()
    fill_buf.dts_recent_PromptLikes = long_value(100)
    if selection_name == 'nh_THU':
        fill_buf.num_recent_ADevents_all = unsigned_int_value()
        fill_buf.dts_recent_ADevents_all = long_value(100)

    # Initialize the new TTree so that each TBranch reads from the
    # appropriate TreeBuffer attribute
    outdata.Branch('run', fill_buf.run, 'run/i')
    outdata.Branch('fileno', fill_buf.fileno, 'fileno/i')
    outdata.Branch('triggerNumber', fill_buf.triggerNumber, 'triggerNumber/I')
    outdata.Branch('timestamp_seconds', fill_buf.timestamp_seconds,
            'timestamp_seconds/I')
    outdata.Branch('timestamp_nanoseconds', fill_buf.timestamp_nanoseconds,
            'timestamp_nanoseconds/I')
    outdata.Branch('timestamp', fill_buf.timestamp, 'timestamp/L')
    outdata.Branch('detector', fill_buf.detector, 'detector/I')
    outdata.Branch('site', fill_buf.site, 'site/I')
    outdata.Branch('triggerType', fill_buf.triggerType, 'triggerType/i')
    outdata.Branch('nHit', fill_buf.nHit, 'nHit/I')
    outdata.Branch('charge', fill_buf.charge, 'charge/F')
    outdata.Branch('fQuad', fill_buf.fQuad, 'fQuad/F')
    outdata.Branch('fMax', fill_buf.fMax, 'fMax/F')
    outdata.Branch('fPSD_t1', fill_buf.fPSD_t1, 'fPSD_t1/F')
    outdata.Branch('fPSD_t2', fill_buf.fPSD_t2, 'fPSD_t2/F')
    outdata.Branch('f2inch_maxQ', fill_buf.f2inch_maxQ, 'f2inch_maxQ/F')
    outdata.Branch('energy', fill_buf.energy, 'energy/F')
    outdata.Branch('x', fill_buf.x, 'x/F')
    outdata.Branch('y', fill_buf.y, 'y/F')
    outdata.Branch('z', fill_buf.z, 'z/F')
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
    outdata.Branch('nHit_previous_WSMuon', fill_buf.nHit_previous_WSMuon,
            'nHit_previous_WSMuon/I')
    outdata.Branch('dt_next_WSMuon', fill_buf.dt_next_WSMuon,
            'dt_next_WSMuon/L')
    outdata.Branch('nHit_next_WSMuon', fill_buf.nHit_next_WSMuon,
            'nHit_next_WSMuon/I')
    outdata.Branch('dt_previous_ADMuon', fill_buf.dt_previous_ADMuon,
            'dt_previous_ADMuon/L')
    outdata.Branch('charge_previous_ADMuon', fill_buf.charge_previous_ADMuon,
            'charge_previous_ADMuon/F')
    outdata.Branch('dt_previous_ShowerMuon', fill_buf.dt_previous_ShowerMuon,
            'dt_previous_ShowerMuon/L')
    outdata.Branch('charge_previous_ShowerMuon',
            fill_buf.charge_previous_ShowerMuon,
            'charge_previous_ShowerMuon/F')
    outdata.Branch('dt_previous_PromptLike',
            fill_buf.dt_previous_PromptLike, 'dt_previous_PromptLike/L')
    outdata.Branch('energy_previous_PromptLike',
            fill_buf.energy_previous_PromptLike,
            'energy_previous_PromptLike/F')
    outdata.Branch('muonVeto_previous_PromptLike',
            fill_buf.muonVeto_previous_PromptLike,
            'muonVeto_previous_PromptLike/i')
    outdata.Branch('x_previous_PromptLike', fill_buf.x_previous_PromptLike,
            'x_previous_PromptLike/F')
    outdata.Branch('y_previous_PromptLike', fill_buf.y_previous_PromptLike,
            'y_previous_PromptLike/F')
    outdata.Branch('z_previous_PromptLike', fill_buf.z_previous_PromptLike,
            'z_previous_PromptLike/F')
    outdata.Branch('dr_previous_PromptLike',
            fill_buf.dr_previous_PromptLike, 'dr_previous_PromptLike/F')
    if selection_name == 'nh_THU':
        outdata.Branch('dt_next_ADevent', fill_buf.dt_next_ADevent,
                'dt_next_ADevent/L')
    outdata.Branch('dt_next_DelayedLike', fill_buf.dt_next_DelayedLike,
            'dt_next_DelayedLike/L')
    outdata.Branch('num_ShowerMuons_5sec', fill_buf.num_ShowerMuons_5sec,
            'num_ShowerMuons_5sec/i')
    outdata.Branch('dts_ShowerMuons_5sec', fill_buf.dts_ShowerMuons_5sec,
            'dts_ShowerMuons_5sec[num_ShowerMuons_5sec]/L')
    if selection_name == 'nh_THU':
        name_str = 'ADevents_800us_restricted'
    elif selection_name == 'nh_DMC':
        name_str = 'PromptLikes_800us'
    elif selection_name == 'ngd':
        name_str = 'PromptLikes_400us'
    else:
        raise ValueError('Invalid selection name: ' + selection_name)
    outdata.Branch('num_' + name_str,
            fill_buf.num_recent_PromptLikes, 'num_' + name_str + '/i')
    outdata.Branch('dts_' + name_str,
            fill_buf.dts_recent_PromptLikes, 'dts_' + name_str +
            '[num_' + name_str + ']/L')
    if selection_name == 'nh_THU':
        outdata.Branch('num_ADevents_800us_all',
                fill_buf.num_recent_ADevents_all,
                'num_ADevents_800us_all/i')
        outdata.Branch('dts_ADevents_800us_all',
                fill_buf.dts_recent_ADevents_all,
                'dts_ADevents_800us_all[num_ADevents_800us_all]/L')
    return outdata, fill_buf

class ProcessHelper(object):
    def __init__(self, selection_name):
        self.MUON_COUNT_TIME = 5*10**9  # 5 seconds, in nanoseconds
        if 'nh' in selection_name:
            self.PROMPT_COUNT_TIME = int(800e3)
        else:
            self.PROMPT_COUNT_TIME = int(400e3)  # 400 us, in nanoseconds
        self.run = 0
        self.fileno = 0
        # Initialize the machinery for dt_next_* and dt_previous_*
        # attributes. The event_cache stores the TreeBuffer for each event
        # until the next WSMuon and next DelayedLike have been found, thus
        # preserving the event order between the input and output TTrees.
        self.event_cache = deque()
        self.first_index_without_next_WSMuon = 0
        self.go_through_whole_cache_WSMuon = True
        self.go_through_whole_cache_ADevent = True
        self.last_WSMuon_time = 0
        self.last_WSMuon_nHit = 0
        self.last_ADMuon_time = {n:0 for n in range(9)}
        self.last_ADMuon_charge = {n:0 for n in range(9)}
        self.last_ShowerMuon_time = {n:0 for n in range(9)}
        self.last_ShowerMuon_charge = {n:0 for n in range(9)}
        self.recent_shower_muons = {n:deque() for n in range(9)}
        self.last_PromptLike_time = {n:-1 for n in range(9)}
        self.last_PromptLike_energy = {n:-1 for n in range(9)}
        self.last_PromptLike_muonVeto = {n: False for n in range(9)}
        self.last_PromptLike_x = {n:0 for n in range(9)}
        self.last_PromptLike_y = {n:0 for n in range(9)}
        self.last_PromptLike_z = {n:0 for n in range(9)}
        self.recent_promptlikes = {n:deque([], 100) for n in range(9)}
        self.coincidence_window_start = {n:0 for n in range(9)}
        if selection_name == 'nh_THU':
            # Keep track of all recent ADevents (within a fixed time delay) even though the nh_THU
            # selection requires ignoring any recent ADevents that occur within
            # a previous coincidence window (e.g. t0, t0 + 300us <--- ignored
            # by t0+500us)
            self.recent_ADevents_all = {n:deque([], 100) for n in range(9)}

def main(entries, infile, outfile, selection_name, debug):

    infile = TFile(infile, 'READ')
    indata = infile.Get('slimmed')
    outfile = TFile(outfile, 'RECREATE')
    outdata, fill_buf = create_computed_TTree('data', outfile, selection_name)
    out_IBDs, ibd_fill_buf = create_computed_TTree('ibds', outfile,
            selection_name, 'IBD candidates (git: %s)')

    helper = ProcessHelper(selection_name)
    if entries == -1:
        entries = xrange(indata.GetEntries())
    else:
        entries = xrange(entries)
    if debug:
        entries = xrange(10000)
    def callback(event_cache):
        logging.debug(event_cache.noTree_loopIndex[0])
    for event_number in entries:
        logging.debug('event %d', event_number)
        # Load the current event in the input TTree
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)
        fetch_indata(indata, fill_buf)

        one_iteration(event_number, outdata, fill_buf, out_IBDs,
                ibd_fill_buf, helper, selection_name, callback)
    finish_emptying_cache(outdata, fill_buf, out_IBDs, ibd_fill_buf,
            helper.event_cache, selection_name, callback)
    outfile.Write()
    outfile.Close()
    infile.Close()

def finish_emptying_cache(outdata, fill_buf, out_IBDs, ibd_fill_buf,
        cache, selection_name, callback):
    if selection_name == 'nh_THU':
        isIBDDelayed = delayeds.isIBDDelayed_nH_THU
    elif selection_name == 'nh_DMC':
        isIBDDelayed = delayeds.isIBDDelayed_nH_DMC
    elif selection_name == 'ngd':
        isIBDDelayed = delayeds.isIBDDelayed
    else:
        raise ValueError('Invalid selection name: ' + selection_name)
    # After the event loop is finished, fill the remaining events from
    # the event_cache into the output TTree
    for cached_event in cache:
        if cached_event.dt_next_WSMuon[0] == -1:
            assign_value(cached_event.dt_next_WSMuon, -1)
            assign_value(cached_event.nHit_next_WSMuon, -1)
            # Leave tag_WSMuonVeto as-is
            # assign_value(cached_event.tag_WSMuonVeto, 0)
        if cached_event.dt_next_DelayedLike[0] == -1:
            assign_value(cached_event.dt_next_DelayedLike, -1)
        e = cached_event
        ibd_delayed = isIBDDelayed(
                e.tag_DelayedLike[0],
                e.dt_previous_PromptLike[0],
                e.dr_previous_PromptLike[0],
                e.muonVeto_previous_PromptLike[0],
                e.num_recent_PromptLikes[0],
                e.dts_recent_PromptLikes,
                e.dt_next_ADevent[0],
                e.dt_next_WSMuon[0],
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
        if ibd_delayed:
            cached_event.copyTo(ibd_fill_buf)
            out_IBDs.Fill()

def fetch_indata(indata, buf):
    # Fetch the necessary values from the input TTree
    assign_value(buf.timestamp, fetch_value(indata, 'timestamp', int))
    assign_value(buf.timestamp_seconds, fetch_value(indata, 'timestamp_seconds',
        int))
    assign_value(buf.timestamp_nanoseconds, fetch_value(indata,
            'timestamp_nanoseconds', int))
    assign_value(buf.triggerType, fetch_value(indata, 'triggerType', int))
    assign_value(buf.detector, fetch_value(indata, 'detector', int))
    assign_value(buf.site, fetch_value(indata, 'site', int))
    assign_value(buf.nHit, fetch_value(indata, 'nHit', int))
    assign_value(buf.charge, fetch_value(indata, 'charge', float))
    assign_value(buf.fMax, fetch_value(indata, 'fMax', float))
    assign_value(buf.fQuad, fetch_value(indata, 'fQuad', float))
    assign_value(buf.fPSD_t1, fetch_value(indata, 'fPSD_t1', float))
    assign_value(buf.fPSD_t2, fetch_value(indata, 'fPSD_t2', float))
    assign_value(buf.f2inch_maxQ, fetch_value(indata, 'f2inch_maxQ', float))
    assign_value(buf.energy, fetch_value(indata, 'energy', float))
    assign_value(buf.x, fetch_value(indata, 'x', float))
    assign_value(buf.y, fetch_value(indata, 'y', float))
    assign_value(buf.z, fetch_value(indata, 'z', float))
    return buf

def one_iteration(event_number, outdata, fill_buf, out_IBDs,
        ibd_fill_buf, helper, selection_name, callback=lambda e:None):
    timestamp = fill_buf.timestamp[0]
    timestamp_seconds = fill_buf.timestamp_seconds[0]
    timestamp_nanoseconds = fill_buf.timestamp_nanoseconds[0]
    triggerType = fill_buf.triggerType[0]
    detector = fill_buf.detector[0]
    site = fill_buf.site[0]
    nHit = fill_buf.nHit[0]
    charge = fill_buf.charge[0]
    fMax = fill_buf.fMax[0]
    fQuad = fill_buf.fQuad[0]
    fPSD_t1 = fill_buf.fPSD_t1[0]
    fPSD_t2 = fill_buf.fPSD_t2[0]
    f2inch_maxQ = fill_buf.f2inch_maxQ[0]
    energy = fill_buf.energy[0]
    x = fill_buf.x[0]
    y = fill_buf.y[0]
    z = fill_buf.z[0]

    if selection_name == 'nh_THU' or selection_name == 'nh_DMC':
        isFlasher = flashers.isFlasher_nH
        muon_event_intensity = energy
        isWSMuon = muons.isWSMuon_nH
        isADMuon = muons.isADMuon_nH
        isShowerMuon = muons.isShowerMuon_nH
        isVetoedByWSMuon = muons.isVetoedByWSMuon_nH
        isVetoedByADMuon = muons.isVetoedByADMuon_nH
        isVetoedByShowerMuon = muons.isVetoedByShowerMuon_nH
        isPromptLike = prompts.isPromptLike_nH
        if selection_name == 'nh_THU':
            isIBDDelayed = delayeds.isIBDDelayed_nH_THU
            isDelayedLike = delayeds.isDelayedLike_nH_THU
        else:
            isIBDDelayed = delayeds.isIBDDelayed_nH_DMC
            isDelayedLike = delayeds.isDelayedLike_nH_DMC
    elif selection_name == 'ngd':
        isFlasher = flashers.isFlasher
        muon_event_intensity = charge
        isWSMuon = muons.isWSMuon
        isADMuon = muons.isADMuon
        isShowerMuon = muons.isShowerMuon
        isVetoedByWSMuon = muons.isVetoedByWSMuon
        isVetoedByADMuon = muons.isVetoedByADMuon
        isVetoedByShowerMuon = muons.isVetoedByShowerMuon
        isPromptLike = prompts.isPromptLike
        isDelayedLike = delayeds.isDelayedLike
        isIBDDelayed = delayeds.isIBDDelayed
    else:
        raise ValueError('Invalid selection name: ' + selection_name)

    logging.debug('event cache size: %d', len(helper.event_cache))

    buf = fill_buf.clone_type()

    # assign copied/input values
    assign_value(buf.noTree_loopIndex, event_number)
    assign_value(buf.run, helper.run)
    assign_value(buf.fileno, helper.fileno)
    assign_value(buf.timestamp, timestamp)
    assign_value(buf.timestamp_seconds, timestamp_seconds)
    assign_value(buf.timestamp_nanoseconds, timestamp_nanoseconds)
    assign_value(buf.triggerType, triggerType)
    assign_value(buf.detector, detector)
    assign_value(buf.site, site)
    assign_value(buf.nHit, nHit)
    assign_value(buf.charge, charge)
    assign_value(buf.fMax, fMax)
    assign_value(buf.fQuad, fQuad)
    assign_value(buf.fPSD_t1, fPSD_t1)
    assign_value(buf.fPSD_t2, fPSD_t2)
    assign_value(buf.f2inch_maxQ, f2inch_maxQ)
    assign_value(buf.energy, energy)
    assign_value(buf.x, x)
    assign_value(buf.y, y)
    assign_value(buf.z, z)

    # Initialize dt_next_WSMuon and dt_next_DelayedLike to -1
    assign_value(buf.dt_next_WSMuon, -1)
    assign_value(buf.dt_next_DelayedLike, -1)
    assign_value(buf.dt_next_ADevent, -1)

    # Compute simple tags and values (those that only require data
    # from the current event)
    event_fID = fID(fMax, fQuad)
    assign_value(buf.fID, event_fID)
    event_fPSD = fPSD(fPSD_t1, fPSD_t2)
    assign_value(buf.fPSD, event_fPSD)
    event_isFlasher = isFlasher(event_fID, event_fPSD, f2inch_maxQ,
            detector)
    assign_value(buf.tag_flasher, event_isFlasher)
    event_isWSMuon = isWSMuon(detector, nHit, triggerType)
    assign_value(buf.tag_WSMuon, event_isWSMuon)
    event_isADMuon = isADMuon(detector, muon_event_intensity)
    assign_value(buf.tag_ADMuon, event_isADMuon)
    event_isShowerMuon = isShowerMuon(detector, muon_event_intensity)
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
            - helper.recent_shower_muons[detector][0] > helper.MUON_COUNT_TIME):
        helper.recent_shower_muons[detector].popleft()

    # Ensure the list of previous prompt-like events only has events
    # more recent than PROMPT_COUNT_TIME ago
    while len(helper.recent_promptlikes[detector]) > 0 and (timestamp -
            helper.recent_promptlikes[detector][0] > helper.PROMPT_COUNT_TIME):
        helper.recent_promptlikes[detector].popleft()

    if selection_name == 'nh_THU':
        while len(helper.recent_ADevents_all[detector]) > 0 and (timestamp -
                helper.recent_ADevents_all[detector][0] > helper.PROMPT_COUNT_TIME):
            helper.recent_ADevents_all[detector].popleft()
        # for the THU analysis, ensure that all recent AD events are more
        # recent than the end of the previous coincidence window and muon
        # veto windows. This is accomplished by removing events after we
        # reach the end of one of these windows.
        if (timestamp - helper.coincidence_window_start[detector] >
                delayeds._NH_THU_MULT_PRE_MIN):
            while len(helper.recent_promptlikes[detector]) > 0 and (
                    helper.recent_promptlikes[detector][0] -
                    helper.coincidence_window_start[detector] <
                    delayeds._NH_THU_MULT_PRE_MIN):
                helper.recent_promptlikes[detector].popleft()
        if (timestamp - helper.last_WSMuon_time >
                muons._NH_WSMUON_VETO_LAST_NS):
            while len(helper.recent_promptlikes[detector]) > 0 and (
                    helper.recent_promptlikes[detector][0] -
                    helper.last_WSMuon_time <
                    muons._NH_WSMUON_VETO_LAST_NS):
                helper.recent_promptlikes[detector].popleft()
        if (timestamp - helper.last_ADMuon_time[detector] >
                muons._NH_ADMUON_VETO_LAST_NS):
            while len(helper.recent_promptlikes[detector]) > 0 and (
                    helper.recent_promptlikes[detector][0] -
                    helper.last_ADMuon_time[detector] <
                    muons._NH_ADMUON_VETO_LAST_NS):
                helper.recent_promptlikes[detector].popleft()
        if (timestamp - helper.last_ShowerMuon_time[detector] >
                muons._NH_SHOWER_MUON_VETO_LAST_NS):
            while len(helper.recent_promptlikes[detector]) > 0 and (
                    helper.recent_promptlikes[detector][0] -
                    helper.last_ShowerMuon_time[detector] <
                    muons._NH_SHOWER_MUON_VETO_LAST_NS):
                helper.recent_promptlikes[detector].popleft()

    # Compute the dts for both recent muons and recent prompt-likes
    recent_showermuon_dts = [timestamp - muon_time for muon_time in
            helper.recent_shower_muons[detector]]
    recent_promptlike_dts = [timestamp - prompt_time for prompt_time
            in helper.recent_promptlikes[detector]]
    if selection_name == 'nh_THU':
        recent_ADevents_all_dts = [timestamp - ADevent_time for ADevent_time in
                helper.recent_ADevents_all[detector]]

    # Assign the values for dt_previous_* and for the lists of
    # recent muons and prompt-likes.
    assign_value(buf.dt_previous_WSMuon, timestamp -
            helper.last_WSMuon_time)
    assign_value(buf.nHit_previous_WSMuon, helper.last_WSMuon_nHit)
    assign_value(buf.dt_previous_ADMuon, timestamp -
            helper.last_ADMuon_time[detector])
    assign_value(buf.charge_previous_ADMuon,
            helper.last_ADMuon_charge[detector])
    assign_value(buf.dt_previous_ShowerMuon, timestamp -
            helper.last_ShowerMuon_time[detector])
    assign_value(buf.charge_previous_ShowerMuon,
            helper.last_ShowerMuon_charge[detector])
    assign_value(buf.dt_previous_PromptLike, timestamp -
            helper.last_PromptLike_time[detector])
    assign_value(buf.energy_previous_PromptLike,
            helper.last_PromptLike_energy[detector])
    assign_value(buf.muonVeto_previous_PromptLike,
            helper.last_PromptLike_muonVeto[detector])
    assign_value(buf.x_previous_PromptLike, helper.last_PromptLike_x[detector])
    assign_value(buf.y_previous_PromptLike, helper.last_PromptLike_y[detector])
    assign_value(buf.z_previous_PromptLike, helper.last_PromptLike_z[detector])
    assign_value(buf.dr_previous_PromptLike, math.sqrt(
        (buf.x_previous_PromptLike[0]-x)**2
        + (buf.y_previous_PromptLike[0]-y)**2
        + (buf.z_previous_PromptLike[0]-z)**2
    ))
    assign_value(buf.num_ShowerMuons_5sec,
            len(helper.recent_shower_muons[detector]))
    assign_value(buf.num_recent_PromptLikes,
            len(helper.recent_promptlikes[detector]))
    for i, dt in enumerate(recent_showermuon_dts):
        assign_value(buf.dts_ShowerMuons_5sec, dt, i)
    for i, dt in enumerate(recent_promptlike_dts):
        assign_value(buf.dts_recent_PromptLikes, dt, i)
    if selection_name == 'nh_THU':
        assign_value(buf.num_recent_ADevents_all,
                len(helper.recent_ADevents_all[detector]))
        for i, dt in enumerate(recent_ADevents_all_dts):
            assign_value(buf.dts_recent_ADevents_all, dt, i)


    # Compute muon vetoes that only rely on previous events
    assign_value(buf.tag_ADMuonVeto,
            isVetoedByADMuon(buf.dt_previous_ADMuon[0]))
    assign_value(buf.tag_ShowerMuonVeto,
            isVetoedByShowerMuon(buf.dt_previous_ShowerMuon[0]))

    # Compute a temporary tag for WS Muon veto based only on the
    # dt_previous_WSMuon. It may be updated when the next WS Muon is processed.
    assign_value(buf.tag_WSMuonVeto, isVetoedByWSMuon(buf.dt_previous_WSMuon[0],
        float('inf')))

    # Update the dt_previous_* and dt_next_* values

    if selection_name == 'nh_THU':
        if isADEvent_THU(detector, energy) and not event_isFlasher:
            logging.debug("isADEvent")
            # Update the coincidence window start time, if this event is not
            # already within a coincicence window
            if (timestamp - helper.coincidence_window_start[detector] >
                    delayeds._NH_THU_MULT_PRE_MIN):
                helper.coincidence_window_start[detector] = timestamp
            helper.recent_promptlikes[detector].append(timestamp)
            helper.recent_ADevents_all[detector].append(timestamp)
            # Assign dt_next_ADevent to events in the event_cache
            for cached_event in helper.event_cache:
                # Some events might already have been assigned, so skip
                # those
                if cached_event.dt_next_ADevent[0] != -1:
                    continue
                # ADevents are restricted to a single AD
                if cached_event.detector[0] == detector:
                    logging.debug('taggedNextADevent%d',
                            cached_event.timestamp[0])
                    assign_value(cached_event.dt_next_ADevent,
                            timestamp
                            - cached_event.timestamp[0])

    # This comes after values are assigned to the buffer because we
    # don't want dt_previous_* to be 0.
    if event_isWSMuon:
        logging.debug("isWSMuon")
        helper.last_WSMuon_time = timestamp
        helper.last_WSMuon_nHit = nHit
        # Decide which events in the event cache to examine, based on whether
        # the event cache has been modified recently
        if helper.go_through_whole_cache_WSMuon:
            # Then the cache has changed substantially since the last
            # iteration, so we will just eat the cost of going through the
            # whole thing.
            pass
        else:
            # Save time by skipping the first chunk of events in the cache that
            # have already been assigned a next_WSMuon. This is achieved by
            # "rotating" the event_cache (a deque type) which means shifting
            # the first n entries to the end, so we start with the new
            # events that need to be labeled. When we reach an already-labeled
            # event, we will be able to exit the loop and un-rotate the cache.
            safe_start_index = helper.first_index_without_next_WSMuon
            helper.event_cache.rotate(-1 * safe_start_index)
        # Assign dt_next_WSMuon to the events in the event_cache
        for i, cached_event in enumerate(helper.event_cache):
            # Look for events that have already been assigned a next_WSMuon
            if cached_event.dt_next_WSMuon[0] != -1:
                if helper.go_through_whole_cache_WSMuon:
                    # We are going through every event in the cache, even if
                    # it's already been assigned
                    continue
                else:
                    # We have reached the end of the "new" events in the cache
                    # and now we can exit the loop.
                    break
            # WSMuons are common across the entire EH/site
            if cached_event.site[0] == site:
                logging.debug('taggedNextWS%d', cached_event.timestamp[0])
                assign_value(cached_event.dt_next_WSMuon,
                        timestamp - cached_event.timestamp[0])
                assign_value(cached_event.nHit_next_WSMuon, nHit)
                assign_value(cached_event.tag_WSMuonVeto,
                        isVetoedByWSMuon(cached_event.dt_previous_WSMuon[0],
                            cached_event.dt_next_WSMuon[0]))
        # Save the last index reached as the new starting location for next
        # time
        if helper.go_through_whole_cache_WSMuon:
            try:
                helper.first_index_without_next_WSMuon = i + 1
                helper.go_through_whole_cache_WSMuon = False
            except UnboundLocalError:
                helper.first_index_without_next_WSMuon = 0
                helper.go_through_whole_cache_WSMuon = True
        else:
            helper.event_cache.rotate(safe_start_index)
            try:
                helper.first_index_without_next_WSMuon = i + safe_start_index
            except UnboundLocalError:
                helper.first_index_without_next_WSMuon = 0
                helper.go_through_whole_cache_WSMuon = True

    if event_isADMuon:
        helper.last_ADMuon_time[detector] = timestamp
        helper.last_ADMuon_charge[detector] = charge
    if event_isShowerMuon:
        helper.last_ShowerMuon_time[detector] = timestamp
        helper.last_ShowerMuon_charge[detector] = charge
        helper.recent_shower_muons[detector].append(timestamp)
    # Don't include prompt-like flasher events
    if event_isPromptLike and not event_isFlasher:
        helper.last_PromptLike_time[detector] = timestamp
        helper.last_PromptLike_energy[detector] = energy
        # Figure out if this event is too soon after a recent muon. This
        # particular veto is used only in certain analyses. The WSMuon
        # term is weird because I don't care about dt_next_WSMuon, so I
        # can use a shortcut veto check function.
        helper.last_PromptLike_muonVeto[detector] = (buf.tag_ADMuonVeto[0]
                + buf.tag_ShowerMuonVeto[0]
                + isVetoedByWSMuon(buf.dt_previous_WSMuon[0], float('inf')) > 0)
        helper.last_PromptLike_x[detector] = x
        helper.last_PromptLike_y[detector] = y
        helper.last_PromptLike_z[detector] = z
        if selection_name != 'nh_THU':
            helper.recent_promptlikes[detector].append(timestamp)
    # Don't include prompt-like delayed events
    if event_isDelayedLike and not event_isFlasher:
        logging.debug("isDelayedLike")
        # Assign dt_next_PromptLike to events in the event_cache
        for cached_event in helper.event_cache:
            # Some events might already have been assigned, so skip
            # those
            if cached_event.dt_next_DelayedLike[0] != -1:
                continue
            # PromptLikes are restricted to a single AD
            if cached_event.detector[0] == detector:
                logging.debug('taggedNextDelayed%d',
                        cached_event.timestamp[0])
                assign_value(cached_event.dt_next_DelayedLike,
                        timestamp
                        - cached_event.timestamp[0])

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
        if done_with_cache(cached_event, selection_name):
            logging.debug('event is done with cache')
            num_to_delete += 1
        else:
            logging.debug('event is not done with cache')
            all_done_with_cache = False
    logging.debug('deleting %d events', num_to_delete)
    if num_to_delete > 0:
        helper.go_through_whole_cache_ADevent = True
        helper.go_through_whole_cache_WSMuon = True
    # Remove the oldest events from the cache and fill them into the
    # new TTree, computing the final IBD tags while we're at it.
    for _ in range(num_to_delete):
        cached_event = helper.event_cache.popleft()
        e = cached_event
        ibd_delayed = isIBDDelayed(
                e.tag_DelayedLike[0],
                e.dt_previous_PromptLike[0],
                e.dr_previous_PromptLike[0],
                e.muonVeto_previous_PromptLike[0],
                e.num_recent_PromptLikes[0],
                e.dts_recent_PromptLikes,
                e.dt_next_ADevent[0],
                e.dt_next_WSMuon[0],
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
        if ibd_delayed:
            cached_event.copyTo(ibd_fill_buf)
            out_IBDs.Fill()

    helper.event_cache.append(buf)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file name')
    parser.add_argument('-o', '--output', help='Output file name')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--events', type=int, default=-1)
    parser.add_argument('--selection', choices=['ngd', 'nh_THU', 'nh_DMC'])
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.events, args.input, args.output, args.selection, args.debug)
