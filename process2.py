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
    found_next_DelayedLike = (selection_name == 'nh_THU'
            or detector not in AD_DETECTORS
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
    fill_buf.tag_flasher = unsigned_int_value(buffer_depth)
    fill_buf.tag_WSMuon = unsigned_int_value()
    fill_buf.tag_ADMuon = unsigned_int_value(buffer_depth)
    fill_buf.tag_ShowerMuon = unsigned_int_value(buffer_depth)
    fill_buf.tag_WSMuonVeto = unsigned_int_value(buffer_depth)
    fill_buf.tag_ADMuonVeto = unsigned_int_value(buffer_depth)
    fill_buf.tag_ShowerMuonVeto = unsigned_int_value(buffer_depth)
    fill_buf.tag_AnyMuonVeto = unsigned_int_value(buffer_depth)
    fill_buf.dt_previous_WSMuon = long_value(buffer_depth)
    fill_buf.nHit_previous_WSMuon = int_value(buffer_depth)
    fill_buf.dt_previous_ADMuon = long_value(buffer_depth)
    fill_buf.charge_previous_ADMuon = float_value(buffer_depth)
    fill_buf.dt_previous_ShowerMuon = long_value(buffer_depth)
    fill_buf.charge_previous_ShowerMuon = float_value(buffer_depth)

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
    branch_multiple('fID', 'F')
    branch_multiple('fPSD', 'F')
    branch_multiple('tag_flasher', 'i')
    branch_multiple('tag_WSMuon', 'i')
    branch_multiple('tag_ADMuon', 'i')
    branch_multiple('tag_ShowerMuon', 'i')
    branch_multiple('tag_WSMuonVeto', 'i')
    branch_multiple('tag_ADMuonVeto', 'i')
    branch_multiple('tag_ShowerMuonVeto', 'i')
    branch_multiple('tag_AnyMuonVeto', 'i')
    branch_multiple('dt_previous_WSMuon', 'L')
    branch_multiple('nHit_previous_WSMuon', 'L')
    branch_multiple('dt_previous_ADMuon', 'L')
    branch_multiple('charge_previous_ADMuon', 'F')
    branch_multiple('dt_previous_ShowerMuon', 'L')
    branch_multiple('charge_previous_ShowerMuon', 'F')
    return outdata, fill_buf

def isValidTriggerType(triggerType):
    return (triggerType & 0x1100) > 0

def main_loop(indata, outdata, fill_buf, debug, limit):
    loopIndex = 0
    last_WSMuon_timestamp = 0
    last_ADMuon_timestamp = 0
    last_ShowerMuon_timestamp = 0
    while loopIndex < limit:
        indata.GetEntry(loopIndex)
        if not isValidTriggerType(indata.triggerType):
            loopIndex += 1
            continue
        if indata.detector not in (1, 5, 6):
            loopIndex += 1
            continue
        if not(indata.tag_flasher == 0 or indata.tag_flasher == 2):
            loopIndex += 1
            continue
        isWSMuon = False
        isADMuon = False
        isShowerMuon = False
        if muons.isWSMuon_nH(indata.detector, indata.nHit, indata.triggerType):
            last_WSMuon_timestamp = indata.timestamp
            isWSMuon = True
        if muons.isADMuon_nH(indata.detector, indata.energy):
            last_ADMuon_timestamp = indata.timestamp
            isADMuon = True
        if muons.isShowerMuon_nH(indata.detector, indata.energy):
            last_ShowerMuon_timestamp = indata.timestamp
            isShowerMuon = True
        isAnyMuonVeto = assign_muons(indata, fill_buf, 0,
                last_WSMuon_timestamp, last_ADMuon_timestamp,
                last_ShowerMuon_timestamp,
                isWSMuon, isADMuon, isShowerMuon)
        if isAnyMuonVeto:
            loopIndex += 1
            continue
        detector = indata.detector
        energy = indata.energy
        if not isADEvent_THU(detector, energy):
            loopIndex += 1
            continue
        prompt_timestamp = indata.timestamp
        prompt_x = indata.x
        prompt_y = indata.y
        prompt_z = indata.z
        multiplicity = 1
        assign_value(fill_buf.site, indata.site)
        assign_value(fill_buf.run, indata.run)
        assign_value(fill_buf.fileno, indata.fileno)
        assign_event(indata, fill_buf, 0)
        assign_value(fill_buf.dt_to_prompt, 0, 0)
        assign_value(fill_buf.dr_to_prompt, 0, 0)
        assign_value(fill_buf.loopIndex, loopIndex, 0)
        loopIndex += 1
        if loopIndex == limit:
            break
        indata.GetEntry(loopIndex)
        while indata.timestamp - prompt_timestamp < 400e3:
            if not (isValidTriggerType(indata.triggerType)
                    and indata.detector == 1
                    and isADEvent_THU(indata.detector, indata.energy)):
                loopIndex += 1
                if loopIndex == limit:
                    break
                indata.GetEntry(loopIndex)
                continue
            logging.debug('  Inner loop: %d', loopIndex)
            multiplicity += 1
            assign_event(indata, fill_buf, multiplicity-1)
            assign_value(fill_buf.dt_to_prompt,
                    indata.timestamp - prompt_timestamp,
                    multiplicity-1)
            assign_value(fill_buf.dr_to_prompt,
                    math.sqrt(
                        (indata.x-prompt_x)**2 +
                        (indata.y-prompt_y)**2 +
                        (indata.z-prompt_z)**2),
                    multiplicity-1)
            assign_value(fill_buf.loopIndex, loopIndex, multiplicity-1)
            loopIndex += 1
            if loopIndex == limit:
                break
            indata.GetEntry(loopIndex)
        assign_value(fill_buf.multiplicity, multiplicity)
        logging.debug(fill_buf.multiplicity)
        outdata.Fill()



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

def assign_muons(source, buf, index, time_last_WSMuon, time_last_ADMuon,
        time_last_ShowerMuon, isWSMuon, isADMuon, isShowerMuon):
    assign_value(buf.tag_WSMuon, isWSMuon)
    assign_value(buf.tag_ADMuon, isADMuon)
    assign_value(buf.tag_ShowerMuon, isShowerMuon)
    isWSMuonVeto = muons.isVetoedByWSMuon_nH(source.timestamp -
            time_last_WSMuon)
    isADMuonVeto = muons.isVetoedByADMuon_nH(source.timestamp -
            time_last_ADMuon)
    isShowerMuonVeto = muons.isVetoedByShowerMuon_nH(source.timestamp -
            time_last_ShowerMuon)
    isAnyMuonVeto = isWSMuonVeto or isADMuonVeto or isShowerMuonVeto
    assign_value(buf.tag_WSMuonVeto, isWSMuonVeto)
    assign_value(buf.tag_ADMuonVeto, isADMuonVeto)
    assign_value(buf.tag_ShowerMuonVeto, isShowerMuonVeto)
    assign_value(buf.tag_AnyMuonVeto, isAnyMuonVeto)
    return isAnyMuonVeto

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


def main(entries, infile, outfile, debug):

    infile = TFile(infile, 'READ')
    indata = infile.Get('computed')
    prepare_indata_branches(indata)
    outfile = TFile(outfile, 'RECREATE')
    ttree_name = 'ad_events'
    ttree_description = 'AD events (git: %s)'
    outdata, fill_buf = create_computed_TTree(ttree_name, outfile,
            ttree_description)

    if entries == -1:
        entries = indata.GetEntries()
    main_loop(indata, outdata, fill_buf, debug, entries)
    outfile.Write()
    outfile.Close()
    infile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file name')
    parser.add_argument('-o', '--output', help='Output file name')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--events', type=int, default=-1)
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.events, args.input, args.output, args.debug)
