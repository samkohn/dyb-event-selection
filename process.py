'''
Process the basic data by computing various cuts, quantities and tags.

'''
from __future__ import print_function

from collections import deque
import time
import argparse
import logging
import math
import os
import os.path
import json
import random
from functools import lru_cache

from common import *
from flashers import fID, fPSD
import flashers
import muons
import prompts
import delayeds
from adevent import isADEvent_THU
from translate import (TreeBuffer, float_value, assign_value,
        fetch_value, int_value, unsigned_int_value, long_value, git_describe,
        create_data_TTree, copy as raw_copy, initialize_indata_onefile as
        raw_initialize)

def create_computed_TTree(name, host_file, selection_name, title=None,
        flat=False):
    from ROOT import TTree
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
    buffer_depth = 10
    fill_buf = TreeBuffer()
    fill_buf.loopIndex = unsigned_int_value(buffer_depth)
    fill_buf.multiplicity = unsigned_int_value()
    fill_buf.timestamp = long_value(buffer_depth)
    fill_buf.timestamp_seconds = int_value(buffer_depth)
    fill_buf.timestamp_nanoseconds = int_value(buffer_depth)
    fill_buf.detector = int_value(buffer_depth)
    fill_buf.dt_to_prompt = long_value(buffer_depth)
    fill_buf.dr_to_prompt = float_value(buffer_depth)
    fill_buf.dt_cluster_to_prev_ADevent = long_value()
    fill_buf.site = int_value()
    fill_buf.run = unsigned_int_value()
    fill_buf.fileno = unsigned_int_value()
    fill_buf.triggerNumber = int_value(buffer_depth)
    fill_buf.triggerType = unsigned_int_value(buffer_depth)
    fill_buf.nHit = int_value(buffer_depth)
    fill_buf.charge = float_value(buffer_depth)
    fill_buf.fQuad = float_value(buffer_depth)
    fill_buf.fMax = float_value(buffer_depth)
    fill_buf.fPSD_t1 = float_value(buffer_depth)
    fill_buf.fPSD_t2 = float_value(buffer_depth)
    fill_buf.f2inch_maxQ = float_value(buffer_depth)
    fill_buf.energy = float_value(buffer_depth)
    fill_buf.x = float_value(buffer_depth)
    fill_buf.y = float_value(buffer_depth)
    fill_buf.z = float_value(buffer_depth)
    fill_buf.fID = float_value(buffer_depth)
    fill_buf.fPSD = float_value(buffer_depth)
    fill_buf.dt_previous_WSMuon = long_value(buffer_depth)
    fill_buf.dt_previous_ADMuon = long_value(buffer_depth)
    fill_buf.dt_previous_ShowerMuon = long_value(buffer_depth)

    def branch_multiple(name, typecode):
        outdata.Branch(name, getattr(fill_buf, name), '{}[multiplicity]/{}'.format(
            name, typecode))
        return
    def branch(name, typecode):
        outdata.Branch(name, getattr(fill_buf, name), '{}/{}'.format(name,
            typecode))
        return

    if flat:
        branch_multiple = branch

    # Initialize the new TTree so that each TBranch reads from the
    # appropriate TreeBuffer attribute
    branch('multiplicity', 'i')
    branch_multiple('loopIndex', 'i')
    branch_multiple('timestamp', 'L')
    branch_multiple('timestamp_seconds', 'I')
    branch_multiple('timestamp_nanoseconds', 'I')
    branch_multiple('detector', 'I')
    branch_multiple('dt_to_prompt', 'L')
    branch_multiple('dr_to_prompt', 'F')
    branch('dt_cluster_to_prev_ADevent', 'L')
    branch('site', 'I')
    branch('run', 'i')
    branch('fileno', 'i')
    branch_multiple('triggerNumber', 'I')
    branch_multiple('triggerType', 'I')
    branch_multiple('nHit', 'I')
    branch_multiple('charge', 'F')
    branch_multiple('fQuad', 'F')
    branch_multiple('fMax', 'F')
    branch_multiple('fPSD_t1', 'F')
    branch_multiple('fPSD_t2', 'F')
    branch_multiple('f2inch_maxQ', 'F')
    branch_multiple('energy', 'F')
    branch_multiple('x', 'F')
    branch_multiple('y', 'F')
    branch_multiple('z', 'F')
    branch_multiple('dt_previous_WSMuon', 'L')
    branch_multiple('dt_previous_ADMuon', 'L')
    branch_multiple('dt_previous_ShowerMuon', 'L')
    return outdata, fill_buf

def isValidTriggerType(triggerType):
    return (triggerType & 0x1100) > 0

def isValidEvent(indata, helpers):
        #last_ADMuon_timestamp, last_ShowerMuon_timestamp,
        #last_ADMuon_timestamp2, last_ShowerMuon_timestamp2):
    if not isValidTriggerType(indata.triggerType):
        return (False, ['not_adevent', 'trigger: {}'.format(indata.triggerType)])
    if indata.detector not in (1, 2, 5, 6):
        return (False, ['not_adevent', 'detector'])
    event_fID = fID(indata.fMax, indata.fQuad)
    event_fPSD = None
    isFlasher = flashers.isFlasher_nH(event_fID, event_fPSD,
            indata.f2inch_maxQ, indata.detector)
    reasons = []
    if isFlasher == 1:
        reasons.append('flasher')
        reasons.append('8-inch flasher')
    elif isFlasher == 4:
        reasons.append('flasher')
        reasons.append('2-inch flasher')
    if muons.isWSMuon_nH(indata.detector, indata.nHit, indata.triggerType):
        reasons.append('wsmuon')
        reasons.append('anymuon')
    if muons.isADMuon_nH(indata.detector, indata.energy):
        reasons.append('admuon')
        reasons.append('anymuon')
    if muons.isShowerMuon_nH(indata.detector, indata.energy):
        reasons.append('showermuon')
        reasons.append('anymuon')
    if indata.detector == 1:
        h = helpers[0]
        muon_veto_list = isMuonVeto(indata, h.shared_last_WSMuon_timestamp,
                h.last_ADMuon_timestamp, h.last_ShowerMuon_timestamp)
        reasons.extend(muon_veto_list)
    elif indata.detector == 2:
        h = helpers[1]
        muon_veto_list = isMuonVeto(indata, h.shared_last_WSMuon_timestamp,
                h.last_ADMuon_timestamp, h.last_ShowerMuon_timestamp)
        reasons.extend(muon_veto_list)
    detector = indata.detector
    energy = indata.energy
    if not isADEvent_THU(detector, energy):
        reasons.append('not_adevent')
    return (len(reasons) == 0, reasons)

class TimeTracker:
    def __init__(self, start_time, pre_veto):
        self.total_good_time = 0
        self.total_vetoed_time = 0
        self.end_of_veto_window = start_time
        self.last_event_tracked = start_time
        self.pre_veto = pre_veto
        self.start_time = start_time
        self.closed = False
        self.num_windows = 0

    def __repr__(self):
        total_time = self.last_event_tracked - self.start_time
        return 'Good: {}, Vetoed: {}, Now: {}, Window end: {}, Efficiency: {:.3}'.format(
                self.total_good_time, self.total_vetoed_time,
                self.last_event_tracked-self.start_time,
                self.end_of_veto_window-self.start_time,
                self.total_good_time/total_time)

    def export(self):
        total_time = self.total_good_time + self.total_vetoed_time
        return {
                'usable_livetime': self.total_good_time,
                'vetoed_livetime': self.total_vetoed_time,
                'daq_livetime': total_time,
                'good_fraction': self.total_good_time/total_time,
                'num_veto_windows': self.num_windows,
        }

    def close(self, end_time):
        if self.closed:
            raise ValueError('Already closed out. Cannot close twice!')
        if end_time < self.end_of_veto_window:
            self.total_vetoed_time += end_time - self.last_event_tracked
        else:
            dt_last_event_to_window_end = (self.end_of_veto_window
                    - self.last_event_tracked)
            dt_window_end_to_now = end_time - self.end_of_veto_window
            self.total_vetoed_time += dt_last_event_to_window_end
            self.total_good_time += dt_window_end_to_now
        self.closed = True
        return


    def log_new_veto(self, timestamp, veto_length):
        if self.closed:
            raise ValueError('Cannot log new veto after being closed out.')
        pre_veto = self.pre_veto
        potential_new_end = timestamp + veto_length
        # If already in a veto window
        if timestamp - pre_veto < self.end_of_veto_window:
            new_time = timestamp - self.last_event_tracked
            self.total_vetoed_time += new_time
            if potential_new_end < self.end_of_veto_window:
                pass
            else:
                self.end_of_veto_window = potential_new_end
        else:
            dt_last_event_to_window_end = (self.end_of_veto_window
                    - self.last_event_tracked)
            dt_window_end_to_pre_veto = (timestamp - self.end_of_veto_window -
                    pre_veto)
            dt_pre_veto_to_now = pre_veto
            self.total_vetoed_time += dt_last_event_to_window_end
            self.total_vetoed_time += dt_pre_veto_to_now
            self.total_good_time += dt_window_end_to_pre_veto
            self.end_of_veto_window = potential_new_end
            self.num_windows += 1
        self.last_event_tracked = timestamp
        return


def log_WSMuon(tracker, timestamp):
    tracker.log_new_veto(timestamp, muons._NH_WSMUON_VETO_LAST_NS)
    return

def log_ADMuon(tracker, timestamp):
    tracker.log_new_veto(timestamp, muons._NH_ADMUON_VETO_LAST_NS)
    return

def log_ShowerMuon(tracker, timestamp):
    tracker.log_new_veto(timestamp, muons._NH_SHOWER_MUON_VETO_LAST_NS)
    return

class RawFileAdapter():
    ATTR_LOOKUP = {
            'triggerNumber': ('triggerNumber', int),
            'timestamp_seconds': ('context.mTimeStamp.mSec', int),
            'timestamp_nanoseconds': ('context.mTimeStamp.mNanoSec',
                int),
            'detector': ('context.mDetId', int),
            'nHit': ('nHit', int),
            'charge': ('NominalCharge', float),
            'fQuad': ('Quadrant', float),
            'fMax': ('MaxQ', float),
            'fPSD_t1': ('time_PSD', float),
            'fPSD_t2': ('time_PSD1', float),
            'f2inch_maxQ': ('MaxQ_2inchPMT', float),
            'triggerType': ('triggerType', int),
            'energy': ('energy', float),
            'x': ('x', float),
            'y': ('y', float),
            'z': ('z', float),
    }
    def __init__(self, ttree_w_friend, run, fileno):
        from ROOT import TFile
        self.ttree = ttree_w_friend
        self.run = run
        self.fileno = fileno
        self.ttree.SetBranchStatus('execNumber', 1)
        # Hack to extract site number
        old_status = self.ttree.GetBranchStatus('context.mSite')
        if int(old_status) == 0:
            self.ttree.SetBranchStatus('context.mSite', 1)
        self.ttree.GetEntry(0)
        self.site = fetch_value(self.ttree, 'context.mSite', int)
        self.ttree.SetBranchStatus('context.mSite', int(old_status))
        # End hack

    def GetEntry(self, index):
        self.ttree.GetEntry(index)
        self.__getattr__.cache_clear()

    def GetEntries(self):
        return self.ttree.GetEntries()

    def GetListOfBranches(self):
        return self.ttree.GetListOfBranches()

    def _timestamp(self):
        sec = self.timestamp_seconds
        nano = self.timestamp_nanoseconds
        return sec*1000000000 + nano

    @lru_cache(maxsize=32)
    def __getattr__(self, name):
        if name == 'timestamp':
            result = self._timestamp()
            return result
        attr_name = self.ATTR_LOOKUP.get(name, None)
        if attr_name is None:
            raise AttributeError('No attribute "{}"'.format(name))
        return fetch_value(self.ttree, attr_name[0], attr_name[1])

    def SetBranchStatus(self, *args):
        self.ttree.SetBranchStatus(*args)

def update_muons(timestamp, reasons, time_tracker, last_WSMuon, last_ADMuon,
        last_ShowerMuon):
    if 'wsmuon' in reasons:
        last_WSMuon = timestamp
        log_WSMuon(time_tracker, timestamp)
    if 'admuon' in reasons:
        last_ADMuon = timestamp
        log_ADMuon(time_tracker, timestamp)
    if 'showermuon' in reasons:
        last_ShowerMuon = timestamp
        log_ShowerMuon(time_tracker, timestamp)
    return (last_WSMuon, last_ADMuon, last_ShowerMuon)

class CoincidenceHelper:
    shared_last_WSMuon_timestamp = 0
    def __init__(self):
        self.last_ADMuon_timestamp = 0
        self.last_ShowerMuon_timestamp = 0
        self.last_ADevent_timestamp = 0
        self.in_coincidence_window = False
        self.prompt_timestamp = 0
        self.prompt_x = 0
        self.prompt_y = 0
        self.prompt_z = 0
        self.multiplicity = 0

def fill_muon(ttree, in_ttree, fill_buf, helper):
    assign_value(fill_buf.multiplicity, 1)
    assign_value(fill_buf.timestamp, in_ttree.timestamp)
    assign_value(fill_buf.timestamp_seconds, in_ttree.timestamp_seconds)
    assign_value(fill_buf.timestamp_nanoseconds,
            in_ttree.timestamp_nanoseconds)
    assign_value(fill_buf.detector, in_ttree.detector)
    assign_value(fill_buf.dt_to_prompt, 0)
    assign_value(fill_buf.dr_to_prompt, 0)
    assign_value(fill_buf.dt_cluster_to_prev_ADevent, 0)
    assign_value(fill_buf.triggerNumber, in_ttree.triggerNumber)
    assign_value(fill_buf.triggerType, in_ttree.triggerType)
    assign_value(fill_buf.nHit, in_ttree.nHit)
    assign_value(fill_buf.fQuad, in_ttree.fQuad)
    assign_value(fill_buf.fMax, in_ttree.fMax)
    assign_value(fill_buf.fPSD_t1, in_ttree.fPSD_t1)
    assign_value(fill_buf.fPSD_t2, in_ttree.fPSD_t2)
    assign_value(fill_buf.f2inch_maxQ, in_ttree.f2inch_maxQ)
    assign_value(fill_buf.energy, in_ttree.energy)
    assign_value(fill_buf.x, in_ttree.x)
    assign_value(fill_buf.y, in_ttree.y)
    assign_value(fill_buf.z, in_ttree.z)
    assign_value(fill_buf.dt_previous_WSMuon, in_ttree.timestamp -
            helper.shared_last_WSMuon_timestamp)
    assign_value(fill_buf.dt_previous_ADMuon, in_ttree.timestamp -
            helper.last_ADMuon_timestamp)
    assign_value(fill_buf.dt_previous_ShowerMuon, in_ttree.timestamp -
            helper.last_ShowerMuon_timestamp)
    ttree.Fill()


def fill_flat(ttree, fill_buf):
    multiplicity = fill_buf.multiplicity[0]
    buf_clone = fill_buf.clone()
    assign_value(fill_buf.multiplicity, 1)
    for i in range(multiplicity):
        assign_value(fill_buf.loopIndex, buf_clone.loopIndex[i])
        assign_value(fill_buf.timestamp, buf_clone.timestamp[i])
        assign_value(fill_buf.timestamp_seconds, buf_clone.timestamp_seconds[i])
        assign_value(fill_buf.timestamp_nanoseconds, buf_clone.timestamp_nanoseconds[i])
        assign_value(fill_buf.detector, buf_clone.detector[i])
        assign_value(fill_buf.dt_to_prompt, 0)
        assign_value(fill_buf.dr_to_prompt, 0)
        assign_value(fill_buf.dt_cluster_to_prev_ADevent, 0)
        assign_value(fill_buf.triggerNumber, buf_clone.triggerNumber[i])
        assign_value(fill_buf.triggerType, buf_clone.triggerType[i])
        assign_value(fill_buf.nHit, buf_clone.nHit[i])
        assign_value(fill_buf.fQuad, buf_clone.fQuad[i])
        assign_value(fill_buf.fMax, buf_clone.fMax[i])
        assign_value(fill_buf.fPSD_t1, buf_clone.fPSD_t1[i])
        assign_value(fill_buf.fPSD_t2, buf_clone.fPSD_t2[i])
        assign_value(fill_buf.f2inch_maxQ, buf_clone.f2inch_maxQ[i])
        assign_value(fill_buf.energy, buf_clone.energy[i])
        assign_value(fill_buf.x, buf_clone.x[i])
        assign_value(fill_buf.y, buf_clone.y[i])
        assign_value(fill_buf.z, buf_clone.z[i])
        assign_value(fill_buf.dt_previous_WSMuon, buf_clone.dt_previous_WSMuon[i])
        assign_value(fill_buf.dt_previous_ADMuon, buf_clone.dt_previous_ADMuon[i])
        assign_value(fill_buf.dt_previous_ShowerMuon,
                buf_clone.dt_previous_ShowerMuon[i])
        ttree.Fill()

def fill_AD_muon(ttree, indata, fill_buf, helper):
    '''Fill the current muon event assuming it's in an AD (i.e. AD or Shower
       muon)
    '''
    fill_muon(ttree, indata, fill_buf, helper)

def fill_WS_muon(ttrees, indata, fill_bufs, helpers):
    '''Fill the current WS muon event into all output TTrees.
    '''
    for ttree, fill_buf, helper in zip(ttrees, fill_bufs, helpers):
        fill_muon(ttree, indata, fill_buf, helper)

def main_loop(indata, outdatas, fill_bufs, debug, limit, flatten):
    useStoredLoopIndex = 'loopIndex' in indata.GetListOfBranches()
    loopIndex = 0
    indata.GetEntry(loopIndex)
    trackers = [TimeTracker(indata.timestamp, delayeds._NH_THU_DT_MAX) for _ in
            range(2)]
    helpers = [CoincidenceHelper() for _ in range(2)]
    while loopIndex < limit:
        indata.GetEntry(loopIndex)
        isValid, reasons = isValidEvent(indata, helpers)
        if not isValid:
            #print('{}: det {} | in window {} | in window {} | {}'.format(
                #loopIndex, indata.detector, in_coincidence_window,
                #in_coincidence_window2, reasons))
            #input('  window1_dt {} | window2_dt {}'.format((indata.timestamp -
                #prompt_timestamp)/1e3, (indata.timestamp -
                    #prompt_timestamp2)/1e3))
            isMuon = 'anymuon' in reasons
            det = indata.detector
            if det in (1, 2):
                index = det - 1
                h = helpers[index]
                outdata = outdatas[index]
                fill_buf = fill_bufs[index]
                tracker = trackers[index]
                if h.in_coincidence_window and isMuon and (
                        indata.timestamp - h.prompt_timestamp <
                        delayeds._NH_THU_DT_MAX):
                    #input(('  leaving window {} because of muon veto').format(
                        #det))
                    h.in_coincidence_window = False
                    h.multiplicity = 0
                elif h.in_coincidence_window and (indata.timestamp -
                        h.prompt_timestamp >= delayeds._NH_THU_DT_MAX):
                    #input(('  leaving window {} because of coinc time'.format(
                        #det))
                    assign_value(fill_buf.multiplicity, h.multiplicity)
                    if flatten:
                        fill_flat(outdata, fill_buf)
                    else:
                        outdata.Fill()
                    h.in_coincidence_window = False
                    h.multiplicity = 0
                if isMuon:
                    if flatten:
                        assign_value(fill_buf.loopIndex, loopIndex)
                        fill_muon(outdata, indata, fill_buf, h)
                    new_timestamps = update_muons(indata.timestamp, reasons,
                            tracker, h.shared_last_WSMuon_timestamp,
                            h.last_ADMuon_timestamp,
                            h.last_ShowerMuon_timestamp)
                    h.last_ADMuon_timestamp = new_timestamps[1]
                    h.last_ShowerMuon_timestamp = new_timestamps[2]
                if ((isMuon or ('anymuonveto' in reasons))
                        and ('not_adevent' not in reasons)):
                    h.last_ADevent_timestamp = indata.timestamp
            elif det in (5, 6) and isMuon:
                for h, outdata, fill_buf, tracker in zip(helpers, outdatas,
                        fill_bufs, trackers):
                    if flatten:
                        assign_value(fill_buf.loopIndex, loopIndex)
                        fill_muon(outdata, indata, fill_buf, h)
                    if h.in_coincidence_window and (indata.timestamp
                            - h.prompt_timestamp < 400e3):
                        #input(('  leaving window {} because of WSmuon veto').format(
                            #det))
                        h.in_coincidence_window = False
                        h.multiplicity = 0
                    new_timestamps = update_muons(indata.timestamp, reasons,
                            tracker, h.shared_last_WSMuon_timestamp,
                            h.last_ADMuon_timestamp,
                            h.last_ShowerMuon_timestamp)
                    h.shared_last_WSMuon_timestamp = new_timestamps[0]
            loopIndex += 1
            continue
        det = indata.detector
        index = det - 1
        h = helpers[index]
        outdata = outdatas[index]
        fill_buf = fill_bufs[index]
        tracker = trackers[index]
        if not h.in_coincidence_window:
            #input('{}: entering window {}'.format(loopIndex, det))
            h.in_coincidence_window = True
            h.prompt_timestamp = indata.timestamp
            h.prompt_x = indata.x
            h.prompt_y = indata.y
            h.prompt_z = indata.z
            assign_value(fill_buf.dt_cluster_to_prev_ADevent, indata.timestamp -
                    h.last_ADevent_timestamp)
            assign_value(fill_buf.site, indata.site)
            assign_value(fill_buf.run, indata.run)
            assign_value(fill_buf.fileno, indata.fileno)
        if indata.timestamp - h.prompt_timestamp < delayeds._NH_THU_DT_MAX:
            h.multiplicity += 1
            if h.multiplicity == 11:
                raise RuntimeError('multiplicity overflow')
            #input('{}: in window {} | multiplicity {}'.format(loopIndex, det,
                #multiplicity))
            assign_event(indata, fill_buf, h.multiplicity - 1)
            assign_value(fill_buf.dt_previous_WSMuon, indata.timestamp -
                    h.shared_last_WSMuon_timestamp, h.multiplicity - 1)
            assign_value(fill_buf.dt_previous_ADMuon, indata.timestamp -
                    h.last_ADMuon_timestamp, h.multiplicity - 1)
            assign_value(fill_buf.dt_previous_ShowerMuon, indata.timestamp -
                    h.last_ShowerMuon_timestamp, h.multiplicity - 1)
            assign_value(fill_buf.dt_to_prompt,
                    indata.timestamp - h.prompt_timestamp,
                    h.multiplicity-1)
            assign_value(fill_buf.dr_to_prompt, math.sqrt(
                    (indata.x-h.prompt_x)**2 +
                    (indata.y-h.prompt_y)**2 +
                    (indata.z-h.prompt_z)**2),
                    h.multiplicity-1)
            if useStoredLoopIndex:
                index = indata.loopIndex
            else:
                index = loopIndex
            assign_value(fill_buf.loopIndex, index, h.multiplicity - 1)
            h.last_ADevent_timestamp = indata.timestamp
            loopIndex += 1
        else:
            #input('{}: leaving window {} because of coincidence time'.format(
                #loopIndex, det))
            assign_value(fill_buf.multiplicity, h.multiplicity)
            if flatten:
                fill_flat(outdata, fill_buf)
            else:
                outdata.Fill()
            h.in_coincidence_window = False
            h.multiplicity = 0
    for tracker in trackers:
        tracker.close(indata.timestamp)
    return trackers


def copy_to_buffer(ttree, buf, index, name):
    assign_value(getattr(buf, name), getattr(ttree, name), index)
    return

def assign_event(source, buf, index):
    copy_to_buffer(source, buf, index, 'timestamp')
    copy_to_buffer(source, buf, index, 'timestamp_seconds')
    copy_to_buffer(source, buf, index, 'timestamp_nanoseconds')
    copy_to_buffer(source, buf, index, 'detector')
    copy_to_buffer(source, buf, index, 'triggerNumber')
    copy_to_buffer(source, buf, index, 'triggerType')
    copy_to_buffer(source, buf, index, 'nHit')
    copy_to_buffer(source, buf, index, 'charge')
    copy_to_buffer(source, buf, index, 'fQuad')
    copy_to_buffer(source, buf, index, 'fMax')
    copy_to_buffer(source, buf, index, 'fPSD_t1')
    copy_to_buffer(source, buf, index, 'fPSD_t2')
    copy_to_buffer(source, buf, index, 'f2inch_maxQ')
    copy_to_buffer(source, buf, index, 'energy')
    copy_to_buffer(source, buf, index, 'x')
    copy_to_buffer(source, buf, index, 'y')
    copy_to_buffer(source, buf, index, 'z')

def isMuonVeto(source, time_last_WSMuon, time_last_ADMuon,
        time_last_ShowerMuon):
    reasons = []
    dt_previous_WSMuon = source.timestamp - time_last_WSMuon
    isWSMuonVeto = muons.isVetoedByWSMuon_nH(source.timestamp -
            time_last_WSMuon)
    if isWSMuonVeto:
        reasons.append('wsmuonveto')
    dt_previous_ADMuon = source.timestamp - time_last_ADMuon
    isADMuonVeto = muons.isVetoedByADMuon_nH(source.timestamp -
            time_last_ADMuon)
    if isADMuonVeto:
        reasons.append('admuonveto')
    dt_previous_ShowerMuon = source.timestamp - time_last_ShowerMuon
    isShowerMuonVeto = muons.isVetoedByShowerMuon_nH(source.timestamp -
            time_last_ShowerMuon)
    if isShowerMuonVeto:
        reasons.append('showermuonveto')
    isAnyMuonVeto = isWSMuonVeto or isADMuonVeto or isShowerMuonVeto
    if len(reasons) > 0:
        reasons.append('anymuonveto')
    return reasons

def prepare_indata_branches(indata):
    indata.SetBranchStatus('*', 0)
    indata.SetBranchStatus('timestamp', 1)
    indata.SetBranchStatus('timestamp_seconds', 1)
    indata.SetBranchStatus('timestamp_nanoseconds', 1)
    indata.SetBranchStatus('detector', 1)
    indata.SetBranchStatus('site', 1)
    indata.SetBranchStatus('run', 1)
    indata.SetBranchStatus('fileno', 1)
    indata.SetBranchStatus('triggerNumber', 1)
    indata.SetBranchStatus('triggerType', 1)
    indata.SetBranchStatus('nHit', 1)
    indata.SetBranchStatus('charge', 1)
    indata.SetBranchStatus('fQuad', 1)
    indata.SetBranchStatus('fMax', 1)
    indata.SetBranchStatus('fPSD_t1', 1)
    indata.SetBranchStatus('fPSD_t2', 1)
    indata.SetBranchStatus('f2inch_maxQ', 1)
    indata.SetBranchStatus('energy', 1)
    indata.SetBranchStatus('x', 1)
    indata.SetBranchStatus('y', 1)
    indata.SetBranchStatus('z', 1)
    #indata.SetBranchStatus('fID', 1)
    #indata.SetBranchStatus('fPSD', 1)
    #indata.SetBranchStatus('tag_flasher', 1)

def get_ads(run):
    return [1, 2]

def main(entries, infile, outfilename, runfile, is_partially_processed, debug,
        flatten):
    from ROOT import TFile

    run, fileno = runfile
    ads = get_ads(run)
    if is_partially_processed:
        infile = TFile(infile, 'READ')
        indata = infile.Get('ad_events')
        prepare_indata_branches(indata)
    else:
        infile = TFile(infile, 'READ')
        calibStats, adSimple = raw_initialize(infile)
        calibStats.AddFriend(adSimple)
        indata = RawFileAdapter(calibStats, run, fileno)
    outfilesplit = os.path.splitext(outfilename)
    outfilenames = [outfilesplit[0] + '_ad' + str(ad) + outfilesplit[1]
            for ad in ads]
    outfiles = [TFile(outfilename, 'RECREATE') for outfilename in outfilenames]
    ttree_name = 'ad_events'
    ttree_description = 'AD events (git: %s)'
    # Trick to transpose from [[out1, buf1], [out2, buf2], ...] to
    #[[out1, out2, ...], [buf1, buf2, ...]]
    outdatas, fill_bufs = zip(*(create_computed_TTree(ttree_name, outfile,
        ttree_description, flat=flatten) for outfile in outfiles))

    if entries == -1:
        entries = indata.GetEntries()
    trackers = main_loop(indata, outdatas, fill_bufs, debug, entries, flatten)
    for outfile in outfiles:
        outfile.Write()
        outfile.Close()
    if is_partially_processed:
        infile.Close()
    else:
        infile.Close()
    for outfilename, tracker in zip(outfilenames, trackers):
        with open(os.path.splitext(outfilename)[0] + '.json', 'w') as f:
            json.dump(tracker.export(), f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file name')
    parser.add_argument('-o', '--output', help='Output file name')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--events', type=int, default=-1)
    parser.add_argument('--new', action='store_true')
    parser.add_argument('-r', '--runfile', type=int, nargs=2,
        help='<run> <fileno>')
    parser.add_argument('--flat', action='store_true',
            help='create a flat TTree with one trigger per entry '
            '(rather than pairing up coincident events within one entry)')
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.events, args.input, args.output, args.runfile, args.new,
            args.debug, args.flat)
