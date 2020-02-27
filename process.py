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

def create_computed_TTree(name, host_file, selection_name, title=None):
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

def isValidEvent(indata, last_WSMuon_timestamp,
        last_ADMuon_timestamp, last_ShowerMuon_timestamp):
    if not isValidTriggerType(indata.triggerType):
        return (False, ['not_adevent', 'trigger: {}'.format(indata.triggerType)])
    if indata.detector not in (1, 5, 6):
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
    muon_veto_list = isMuonVeto(indata, last_WSMuon_timestamp,
            last_ADMuon_timestamp, last_ShowerMuon_timestamp)
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


def main_loop(indata, outdata, fill_buf, debug, limit):
    loopIndex = 0
    indata.GetEntry(loopIndex)
    time_tracker = TimeTracker(indata.timestamp, delayeds._NH_THU_DT_MAX)
    #time_tracker2 = TimeTracker(indata.timestamp, delayeds._NH_THU_DT_MAX)
    last_WSMuon_timestamp = 0
    last_ADMuon_timestamp = 0
    last_ShowerMuon_timestamp = 0
    last_ADevent_timestamp = 0
    in_coincidence_window = False
    prompt_timestamp = 0
    prompt_x = 0
    prompt_y = 0
    prompt_z = 0
    multiplicity = 0
    while loopIndex < limit:
        indata.GetEntry(loopIndex)
        isValid, reasons = isValidEvent(indata,
                last_WSMuon_timestamp, last_ADMuon_timestamp,
                last_ShowerMuon_timestamp)
        if not isValid:
            causeMuonVeto = 'anymuon' in reasons
            if in_coincidence_window and causeMuonVeto and (indata.timestamp -
                    prompt_timestamp < 400e3):
                in_coincidence_window = False
                multiplicity = 0
            elif in_coincidence_window and (indata.timestamp - prompt_timestamp
                    >= 400e3) and indata.detector == 1:
                assign_value(fill_buf.multiplicity, multiplicity)
                outdata.Fill()
                in_coincidence_window = False
                multiplicity = 0
            new_timestamps = update_muons(indata.timestamp, reasons,
                    time_tracker, last_WSMuon_timestamp, last_ADMuon_timestamp,
                    last_ShowerMuon_timestamp)
            last_WSMuon_timestamp = new_timestamps[0]
            last_ADMuon_timestamp = new_timestamps[1]
            last_ShowerMuon_timestamp = new_timestamps[2]
            if ((('anymuon' in reasons) or ('anymuonveto' in reasons))
                    and ('not_adevent' not in reasons)):
                last_ADevent_timestamp = indata.timestamp
            loopIndex += 1
            continue
        if not in_coincidence_window:
            in_coincidence_window = True
            prompt_timestamp = indata.timestamp
            prompt_x = indata.x
            prompt_y = indata.y
            prompt_z = indata.z
            assign_value(fill_buf.dt_cluster_to_prev_ADevent, indata.timestamp -
                    last_ADevent_timestamp)
            assign_value(fill_buf.site, indata.site)
            assign_value(fill_buf.run, indata.run)
            assign_value(fill_buf.fileno, indata.fileno)
            muonVeto_post_prompt = False
        if indata.timestamp - prompt_timestamp < 400e3:
            multiplicity += 1
            if multiplicity == 11:
                raise RuntimeError('multiplicity overflow')
            assign_event(indata, fill_buf, multiplicity - 1)
            assign_value(fill_buf.dt_previous_WSMuon, indata.timestamp -
                    last_WSMuon_timestamp, multiplicity - 1)
            assign_value(fill_buf.dt_previous_ADMuon, indata.timestamp -
                    last_ADMuon_timestamp, multiplicity - 1)
            assign_value(fill_buf.dt_previous_ShowerMuon, indata.timestamp -
                    last_ShowerMuon_timestamp, multiplicity - 1)
            assign_value(fill_buf.dt_to_prompt,
                    indata.timestamp - prompt_timestamp,
                    multiplicity-1)
            assign_value(fill_buf.dr_to_prompt, math.sqrt(
                    (indata.x-prompt_x)**2 +
                    (indata.y-prompt_y)**2 +
                    (indata.z-prompt_z)**2),
                    multiplicity-1)
            assign_value(fill_buf.loopIndex, loopIndex, multiplicity - 1)
            last_ADevent_timestamp = indata.timestamp
            loopIndex += 1
        else:
            assign_value(fill_buf.multiplicity, multiplicity)
            outdata.Fill()
            in_coincidence_window = False
            multiplicity = 0
    time_tracker.close(indata.timestamp)
    return time_tracker


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
    indata.SetBranchStatus('fID', 1)
    indata.SetBranchStatus('fPSD', 1)
    indata.SetBranchStatus('tag_flasher', 1)


def main(entries, infile, outfilename, runfile, is_partially_processed, debug):
    from ROOT import TFile

    run, fileno = runfile
    if is_partially_processed:
        infile = TFile(infile, 'READ')
        indata = infile.Get('computed')
        prepare_indata_branches(indata)
    else:
        infile = TFile(infile, 'READ')
        calibStats, adSimple = raw_initialize(infile)
        calibStats.AddFriend(adSimple)
        indata = RawFileAdapter(calibStats, run, fileno)
    outfilesplit = os.path.splitext(outfilename)
    outfilename = outfilesplit[0] + '_ad1' + outfilesplit[1]
    #outfilename2 = outfilesplit[0] + '_ad2' + outfilesplit[1]
    outfile = TFile(outfilename, 'RECREATE')
    #outfile2 = TFile(outfilename2, 'RECREATE')
    ttree_name = 'ad_events'
    ttree_description = 'AD events (git: %s)'
    outdata, fill_buf = create_computed_TTree(ttree_name, outfile,
            ttree_description)
    #outdata2, fill_buf2 = create_computed_TTree(ttree_name, outfile2,
            #ttree_description)

    if entries == -1:
        entries = indata.GetEntries()
    tracker = main_loop(indata, outdata, fill_buf, debug, entries)
    outfile.Write()
    outfile.Close()
    #outfile2.Write()
    #outfile2.Close()
    if is_partially_processed:
        infile.Close()
    else:
        infile.Close()
    with open(os.path.splitext(outfilename)[0] + '.json', 'w') as f:
        json.dump(tracker.export(), f)
    #with open(os.path.splitext(outfilename2)[0] + '.json', 'w') as f:
        #json.dump(tracker2.export(), f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file name')
    parser.add_argument('-o', '--output', help='Output file name')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--events', type=int, default=-1)
    parser.add_argument('--new', action='store_true')
    parser.add_argument('-r', '--runfile', type=int, nargs=2,
        help='<run> <fileno>')
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.events, args.input, args.output, args.runfile, args.new, args.debug)
