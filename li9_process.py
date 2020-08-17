'''
Process the basic data by computing various cuts, quantities and tags.

'''
import argparse
from collections import deque
from dataclasses import dataclass
import itertools
import logging
import math
import os
import pprint
import json

from common import *
import muons
from muons import (
        _NH_WSMUON_VETO_LAST_NS as WSMUON_VETO,
        _NH_ADMUON_VETO_LAST_NS as ADMUON_VETO,
        _NH_SHOWER_MUON_VETO_LAST_NS as SHOWER_MUON_VETO
)
import delayeds
from flashers import isFlasher
from adevent import isADEvent_THU, isADEvent_THU_lowenergy
from root_util import (TreeBuffer, float_value, assign_value,
        int_value, unsigned_int_value, long_value)
from translate import git_describe

COINCIDENCE_WINDOW = delayeds._NH_THU_MAX_TIME
LI9_LOWENERGY_MIN = 20
LI9_LOWENERGY_MAX = 1000
LI9_HIGHENERGY_MIN = 2500
LI9_MIDENERGY_MIN = LI9_LOWENERGY_MAX
LI9_MIDENERGY_MAX = LI9_HIGHENERGY_MIN

def li9_lowenergy_muon(detector, energy):
    return (detector in AD_DETECTORS
            and energy > LI9_LOWENERGY_MIN
            and energy < LI9_LOWENERGY_MAX)

def li9_midenergy_muon(detector, energy):
    return (detector in AD_DETECTORS
            and energy > LI9_MIDENERGY_MIN
            and energy < LI9_MIDENERGY_MAX)

def li9_highenergy_muon(detector, energy):
    return (detector in AD_DETECTORS
            and energy > LI9_HIGHENERGY_MIN)

@dataclass
class MuonTimingNTag:
    timestamp: int
    charge: float
    label: str
    ntag: bool
    entry: int

def muon_TTree(host_file):
    from ROOT import TTree
    host_file.cd()
    muons = TTree('muon_tags', 'Muon tags')
    fill_buf = TreeBuffer()
    fill_buf.neutron_tag = unsigned_int_value()
    fill_buf.low_energy = unsigned_int_value()
    fill_buf.mid_energy = unsigned_int_value()
    fill_buf.high_energy = unsigned_int_value()
    fill_buf.entry = unsigned_int_value()
    def branch(name, typecode):
        muons.Branch(name, getattr(fill_buf, name), '{}/{}'.format(name,
            typecode))
        return
    branch('neutron_tag', 'i')
    branch('low_energy', 'i')
    branch('mid_energy', 'i')
    branch('high_energy', 'i')
    branch('entry', 'i')
    return muons, fill_buf


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
    buffer_depth = 20
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
    fill_buf.dt_lowenergy_muon_ntag = long_value()
    fill_buf.dt_midenergy_muon_ntag = long_value()
    fill_buf.dt_highenergy_muon_ntag = long_value()
    fill_buf.dt_lowenergy_muon = long_value()
    fill_buf.dt_midenergy_muon = long_value()
    fill_buf.dt_highenergy_muon = long_value()

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
    branch_multiple('fID', 'F')
    branch_multiple('fPSD', 'F')
    branch('dt_lowenergy_muon', 'L')
    branch('dt_midenergy_muon', 'L')
    branch('dt_highenergy_muon', 'L')
    branch('dt_lowenergy_muon_ntag', 'L')
    branch('dt_midenergy_muon_ntag', 'L')
    branch('dt_highenergy_muon_ntag', 'L')
    return outdata, fill_buf

def isValidTriggerType(triggerType):
    return (triggerType & 0x1100) > 0

def isValidEvent(indata, helpers):
        #last_ADMuon_timestamp, last_ShowerMuon_timestamp,
        #last_ADMuon_timestamp2, last_ShowerMuon_timestamp2):
    if not isValidTriggerType(indata.triggerType):
        return (False, ['not_adevent', 'trigger: {}'.format(indata.triggerType)])
    reasons = []
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
        self.muon_counts = 'uninitialized'

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
                'muon_counts': self.muon_counts,
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
        # If veto window extends from before to after start_time
        # (crosses boundary), then count it as a real window. Otherwise
        # if the window does not extend to after start_time (all before
        # start_time), ignore it.
        #
        # The comparison is to self.end_of_veto_window because this
        # attribute is initialized to start_time but there may be
        # multiple muons with varying veto window sizes overlapping just
        # before start_time. This ensures the window doesn't shrink
        # because e.g. a shower muon is followed by a WP muon.
        if timestamp < self.start_time and potential_new_end > self.end_of_veto_window:
            # This is the start of the run so if we're in this code
            # block, then there's a window and it is the first and only
            # window so far. This direct assignment (rather than increment)
            # neatly handles the case where there are multiple overlapping
            # muon veto windows caused by muons just before self.start_time.
            self.num_windows = 1
            self.end_of_veto_window = potential_new_end
            # Return without messing with self.last_event_tracked or any other
            # logic. Better than wrapping the below in an else block.
            return
        elif timestamp > self.start_time:
            pass  # Execute the rest of this function
        else:  # before start_time and veto window ends before start_time
            return

        # If already in a veto window
        if timestamp - pre_veto < self.end_of_veto_window:
            logging.debug('Already in veto window')
            new_time = timestamp - self.last_event_tracked
            self.total_vetoed_time += new_time
            if potential_new_end < self.end_of_veto_window:
                pass
            else:
                self.end_of_veto_window = potential_new_end
            # This code block also handles the case where a muon occurs just
            # after start_time, and the pre-veto extends to before start_time,
            # and there was no muon just before start_time to create a window
            # that this new muon's pre-veto would overlap with. Certainly if we
            # are in this block, there should be at least 1 window. Without the
            # following line, the above-described situation would result in 0
            # windows.
            self.num_windows = max(1, self.num_windows)
        else:
            logging.debug('new veto window')
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
    def __init__(self):
        self.last_ADevent_timestamp = 0
        self.in_coincidence_window = False
        self.prompt_timestamp = 0
        self.prompt_x = 0
        self.prompt_y = 0
        self.prompt_z = 0
        self.multiplicity = 0

class MuonHelper:
    def __init__(self, muon_ttree, time_tracker, ad_num, out_muons, muon_fill_buf):
        from ROOT import TH1F
        self.ttree = muon_ttree
        self.ad = ad_num
        self.out_muons = out_muons
        self.muon_fill_buf = muon_fill_buf
        self.time_WSMuon = 0
        INIT_TIME = int(1.8e19)  # larger than all Dyb timestamps but fits in C long
        self.time_lowenergy_muon = INIT_TIME
        self.time_midenergy_muon = INIT_TIME
        self.time_highenergy_muon = INIT_TIME
        self.time_lowenergy_muon_ntag = 0
        self.time_midenergy_muon_ntag = 0
        self.time_highenergy_muon_ntag = 0
        self.time_next_muon = None
        self._event_timestamp = None
        self._muon_entry = 0
        self._last_good_muon_entry = -1
        self.ttree.GetEntry(self._muon_entry)
        self.time_tracker = time_tracker
        self.muon_counts = {
                'low': 0,
                'mid': 0,
                'high': 0,
                'low_ntag': 0,
                'mid_ntag': 0,
                'high_ntag': 0,
                }
        self.recent_muons = {
                'low': deque([], 20),
                'mid': deque([], 20),
                'high': deque([], 20),
                'all': deque([]),  # just to keep track of the order
        }
        self.muon_spectrum_all = TH1F(f'muonspec_all', 'muonspec_all', 100, 0, 5e5)
        self.muon_spectrum_ntag = TH1F(f'muonspec_ntag', 'muonspec_ntag', 100, 0, 5e5)
        return

    #@staticmethod
    #def _search_most_recent(iterable, search, key=lambda x:x):
        #last = None
        #for item in iterable:
            #if key(item) > search:
                #return last
            #else:
                #last = item
        #return last

    #def _search_recent_muon(self, timestamp):
        #most_recent = {}
        #logging.debug(pprint.pformat(self.recent_muons))
        #for label, timestamps_charges in self.recent_muons.items():
            #timestamp_charge_pair = self._search_most_recent(timestamps_charges,
                    #timestamp, key=lambda x:x[0])
            #if timestamp_charge_pair is not None:
                #most_recent[label] = timestamp_charge_pair
        #logging.debug(pprint.pformat(most_recent))
        #return max(most_recent.items(), key=lambda pair:pair[1])

    def attempt_neutron_tag_muon(self, ntag_timestamp):
        # Find whether any muons are within (20, 200)us before ntag_timestamp
        for label, muons in self.recent_muons.items():
            no_tags_yet = True
            # Most recent first --> reversed
            for muon in reversed(muons):
                if muon.ntag:
                    no_tags_yet = False
                    continue
                dt = ntag_timestamp - muon.timestamp
                if dt > 20000 and dt < 200000:
                    # TAG THIS MUON
                    muon.ntag = True
                    self.muon_counts[f'{label}_ntag'] += 1
                    if muon.charge >= 5000:
                        # Trying to reproduce Chris's plot which cuts off at 5000 PE
                        self.muon_spectrum_ntag.Fill(muon.charge)
                    if no_tags_yet:
                        setattr(self, f'time_{label}energy_muon_ntag', muon.timestamp)
                        no_tags_yet = False
                    logging.debug('tagging recent muon!')
                    logging.debug('neutron tag event timestamp: %d', ntag_timestamp)
                    logging.debug('muon event: %s', muon)
                    logging.debug('dt: %d', ntag_timestamp - muon.timestamp)
                elif dt > 200000:
                    # Then this event will not tag any muons of this type
                    break  # out of inner loop

    #def neutron_tag_muon(self, muon_timestamp, label, only_increment=False):
        #if label == 'low':
            #self.muon_counts['low_ntag'] += 1
            #if only_increment:
                #logging.debug('Flagged for only increment')
            #else:
                #self.time_lowenergy_muon_ntag = muon_timestamp
            #logging.debug('tagging low energy muon at timestamp %d', muon_timestamp)
        #elif label == 'mid':
            #if self.time_midenergy_muon_ntag != muon_timestamp:
                ## only increment counter for new tags (don't want to
                ## double-count if there are 2 neutrons)
                #self.muon_counts['mid_ntag'] += 1
            #self.time_midenergy_muon_ntag = muon_timestamp
            #logging.debug('tagging last mid energy muon at timestamp %d',
                    #self.time_midenergy_muon_ntag)
        #elif label == 'high':
            #if self.time_highenergy_muon_ntag != muon_timestamp:
                ## only increment counter for new tags (don't want to
                ## double-count if there are 2 neutrons)
                #self.muon_counts['high_ntag'] += 1
            #self.time_highenergy_muon_ntag = muon_timestamp
            #logging.debug('tagging last high energy muon at timestamp %d',
                    #self.time_highenergy_muon_ntag)
        #else:
            #raise RuntimeError("Can't find that muon")
        #self.muon_spectrum_ntag.Fill(charge)
        #logging.debug('neutron tagged muon has charge %d', charge)
        #logging.debug('  (energy ~ %.01f MeV)', charge / 170)
        #return

    def dt_previous_WSMuon(self):
        return self._event_timestamp - self.time_WSMuon

    def dt_previous_showermuons(self, tags=None):
        if tags is None:
            tags = ['', '_ntag']
        now = self._event_timestamp
        labels = itertools.product(['low', 'mid', 'high'], tags)
        result = {}
        for energy, tag in labels:
            attribute = f'time_{energy}energy_muon{tag}'
            diff = now - getattr(self, attribute)
            if diff < 0:  # this happens when there hasn't yet been a muon
                diff = now  # this is a sentinel value for "no muon yet"
            result[energy, tag] = diff
        return result

    def dt_next_muon(self):
        return self.time_next_muon - self._event_timestamp

    def isVetoed_strict(self):
        logging.debug('In isVetoed_strict')
        logging.debug('  dt next WSMuon: %d', self.dt_next_muon())
        return self.isVetoed() or self.dt_next_muon() < COINCIDENCE_WINDOW

    def isVetoed(self):
        logging.debug('In isVetoed')
        logging.debug('  dt to WSMuon: %d', self.dt_previous_WSMuon())
        logging.debug('  dt to AD/Shower: %d',
                min(self.dt_previous_showermuons().values()))
        return self.dt_previous_WSMuon() < muons._NH_WSMUON_VETO_LAST_NS

    def close(self, last_entry=-1, output=None):
    #def close(self, last_entry=-1):
        from ROOT import TFile
        if last_entry == -1:
            total_entries = self.ttree.GetEntries()
        else:
            total_entries = min(last_entry, self.ttree.GetEntries())
        mu_data = self.ttree
        self._muon_entry = self._last_good_muon_entry + 1
        while self._muon_entry < total_entries:
            mu_data.GetEntry(self._muon_entry)
            is_WS = muons.isWSMuon_nH(mu_data.detector, mu_data.nHit,
                    mu_data.triggerType)
            if is_WS:
                log_WSMuon(self.time_tracker, mu_data.timestamp)
            muon_object = MuonTimingNTag(0, 0, 'not relevant', False, self._muon_entry)
            self.recent_muons['all'].append(muon_object)
            self._muon_entry += 1
        for remaining_muon in self.recent_muons['all']:
            label = remaining_muon.label
            assign_value(self.muon_fill_buf.low_energy, int(label == 'low'))
            assign_value(self.muon_fill_buf.mid_energy, int(label == 'mid'))
            assign_value(self.muon_fill_buf.high_energy, int(label == 'high'))
            assign_value(self.muon_fill_buf.neutron_tag, int(remaining_muon.ntag))
            assign_value(self.muon_fill_buf.entry, remaining_muon.entry)
            self.out_muons.Fill()
            logging.debug('[close] Filling entry %d', remaining_muon.entry)

        if output is not None:
            outfile = TFile(output, 'RECREATE')
            self.muon_spectrum_all.SetDirectory(outfile)
            self.muon_spectrum_ntag.SetDirectory(outfile)
            outfile.Write()
            outfile.Close()


    def load(self, timestamp):
        self._event_timestamp = timestamp
        mu_data = self.ttree
        if timestamp < mu_data.timestamp:
            if self.time_WSMuon == 0:
                # This is during the first few events, when there may be an AD
                # event before the first muon event.
                pass
            else:
                raise RuntimeError('cluster events are out of order')
        mu_data.GetEntry(self._muon_entry)
        total_entries = mu_data.GetEntries()
        # Find the dts to previous muons
        while timestamp > mu_data.timestamp and self._muon_entry < total_entries:
            # This if-statement avoids double-counting the last muon
            # found the previous time load() was called, since we
            # backtrack and repeat the last index used to ensure that mu_data.timestamp <
            # timestamp is always true.
            if self._muon_entry <= self._last_good_muon_entry:
                logging.debug("I believe we've already counted entry %d",
                        self._muon_entry)
                self._muon_entry += 1
                mu_data.GetEntry(self._muon_entry)
                continue
            #if mu_data.timestamp in (self.time_lowenergy_muon,
                    #self.time_midenergy_muon, self.time_highenergy_muon):
                #self._muon_entry += 1
                #mu_data.GetEntry(self._muon_entry)
                #continue
            is_WS = muons.isWSMuon_nH(mu_data.detector, mu_data.nHit,
                    mu_data.triggerType)
            is_lowenergy = (mu_data.detector == self.ad
                    and li9_lowenergy_muon(mu_data.detector,
                        mu_data.energy))
            is_midenergy = (mu_data.detector == self.ad
                    and li9_midenergy_muon(mu_data.detector,
                        mu_data.energy))
            is_highenergy = (mu_data.detector == self.ad
                    and li9_highenergy_muon(mu_data.detector,
                        mu_data.energy))
            is_muon = False
            if is_WS:
                self.time_WSMuon = mu_data.timestamp
                logging.debug('Logging WSMuon, entry %d', mu_data.loopIndex)
                log_WSMuon(self.time_tracker, mu_data.timestamp)
            elif is_lowenergy:
                is_muon = True
                label = 'low'
            elif is_midenergy:
                is_muon = True
                label = 'mid'
            elif is_highenergy:
                is_muon = True
                label = 'high'
            if isFlasher(mu_data.fID, 1, mu_data.f2inch_maxQ, mu_data.detector) & 0b101 > 0:
                # Reject flasher candidates (but ignore PSD cut, hence the bit logic)
                is_muon = False
            if is_muon:
                muon_object = MuonTimingNTag(
                        mu_data.timestamp,
                        mu_data.charge,
                        label,
                        False,
                        self._muon_entry
                )
                if mu_data.charge > 5000:
                    self.muon_spectrum_all.Fill(mu_data.charge)
                self.recent_muons['all'].append(muon_object)
                self.muon_counts[label] += 1
                self.recent_muons[label].append(muon_object)
                setattr(self, f'time_{label}energy_muon', mu_data.timestamp)
                logging.debug('new %s energy muon', label)
            else:
                muon_object = MuonTimingNTag(0, 0, 'not relevant', False, self._muon_entry)
                self.recent_muons['all'].append(muon_object)
            self._last_good_muon_entry = self._muon_entry
            self._muon_entry += 1
            mu_data.GetEntry(self._muon_entry)
        # Test to see if we have simply run out of muons, in which
        # case we do not need to search for a "next" muon. The latest
        # "next" muon is the last muon, meaning there is an automatic
        # veto window at the end of each file due to the pre-veto. Note
        # that the "next muon" timestamp will be earlier than the event
        # timestamp, resulting in all future events being vetoed.
        if timestamp > mu_data.timestamp:  # i.e. *still* greater than
            return
        # Now find the dt to the next WS muon
        not_found = True
        while not_found and self._muon_entry < total_entries:
            mu_data.GetEntry(self._muon_entry)
            is_WS = muons.isWSMuon_nH(mu_data.detector, mu_data.nHit,
                    mu_data.triggerType)
            if is_WS:
                self.time_next_muon = mu_data.timestamp
                not_found = False
            self._muon_entry += 1
        self._muon_entry = self._last_good_muon_entry
        mu_data.GetEntry(self._muon_entry)

        # Process some of the cache by filling it into the muon TTree
        far_enough_past = True
        any_muons_in_cache = len(self.recent_muons['all']) > 0
        while far_enough_past and any_muons_in_cache:
            early_muon = self.recent_muons['all'][0]
            if timestamp - early_muon.timestamp > 200000:
                self.recent_muons['all'].popleft()
                label = early_muon.label
                assign_value(self.muon_fill_buf.low_energy, int(label == 'low'))
                assign_value(self.muon_fill_buf.mid_energy, int(label == 'mid'))
                assign_value(self.muon_fill_buf.high_energy, int(label == 'high'))
                assign_value(self.muon_fill_buf.neutron_tag, int(early_muon.ntag))
                assign_value(self.muon_fill_buf.entry, early_muon.entry)
                self.out_muons.Fill()
                logging.debug('Filling entry %d', early_muon.entry)
            else:
                far_enough_past = False
            any_muons_in_cache = len(self.recent_muons['all']) > 0


def is_neutron_tag(energy):
    return energy > 1.8 and energy < 12


def main_loop(clusters, muons, outdata, fill_buf, debug, limit, out_muons, muon_fill_buf):
#def main_loop(clusters, muons, outdata, fill_buf, debug, limit:
    clusters.GetEntry(0)
    muons.GetEntry(0)
    # Start tracking DAQ and veto time after the "phantom muon" veto window
    max_veto_window = max(WSMUON_VETO, ADMUON_VETO, SHOWER_MUON_VETO)
    t0 = min(clusters.timestamp, muons.timestamp) + max_veto_window
    tracker = TimeTracker(t0, COINCIDENCE_WINDOW)
    helper = CoincidenceHelper()
    muon_helper = MuonHelper(muons, tracker, clusters.detector, out_muons, muon_fill_buf)
    #muon_helper = MuonHelper(muons, tracker, clusters.detector)
    clusters_index = 0
    # Process the initial "phantom muon" veto window without tracking DAQ time
    # or muon veto count, and without processing coincidences.
    while clusters.timestamp < t0 and clusters_index < limit:
        clusters_index += 1
        clusters.GetEntry(clusters_index)
    logging.debug(limit)
    # Set up the neutron tag isolation time window buffer
    ntag_buf = fill_buf.clone_type()
    while clusters_index < limit:
        logging.debug(clusters_index)
        clusters.GetEntry(clusters_index)
        logging.debug(clusters.loopIndex)
        isValid, reasons = isValidEvent(clusters, helper)
        timestamp = clusters.timestamp
        logging.debug(timestamp)
        if isValid:
            muon_helper.load(timestamp)
            # Neutron tag checks:
            if is_neutron_tag(clusters.energy):
                muon_helper.attempt_neutron_tag_muon(clusters.timestamp)
            # Perform the search for double coincidences
            if not helper.in_coincidence_window:
                isVetoed = muon_helper.isVetoed_strict()
                if isVetoed:
                    logging.debug('  -- is vetoed')
                    logging.debug('  WS: %d', muon_helper.dt_previous_WSMuon())
                    logging.debug('  Sh: %s', muon_helper.dt_previous_showermuons())
                    logging.debug('  Ne: %d', muon_helper.dt_next_muon())
                    clusters_index += 1
                    continue
                logging.debug('survives muon veto')
                helper.in_coincidence_window = True
                helper.prompt_timestamp = timestamp
                helper.prompt_x = clusters.x
                helper.prompt_y = clusters.y
                helper.prompt_z = clusters.z
                assign_value(fill_buf.site, clusters.site)
                assign_value(fill_buf.run, clusters.run)
                assign_value(fill_buf.fileno, clusters.fileno)
                assign_value(fill_buf.dt_cluster_to_prev_ADevent,
                        clusters.timestamp - helper.last_ADevent_timestamp)
                assign_muons(fill_buf, muon_helper)
            else:
                isVetoed = muon_helper.isVetoed()
                if isVetoed:
                    logging.debug('  -- is vetoed')
                    logging.debug('  WS: %d', muon_helper.dt_previous_WSMuon())
                    logging.debug('  Sh: %s', muon_helper.dt_previous_showermuons())
                    logging.debug('  Ne: %d', muon_helper.dt_next_muon())
                    clusters_index += 1
                    helper.last_ADevent_timestamp = clusters.timestamp
                    continue
                logging.debug('survives muon veto')
            if timestamp - helper.prompt_timestamp < COINCIDENCE_WINDOW:
                helper.multiplicity += 1
                logging.debug('  multiplicity = %d', helper.multiplicity)
                if helper.multiplicity == 21:
                    helper.multiplicity -= 1
                assign_event(clusters, fill_buf, helper.multiplicity - 1)
                assign_value(fill_buf.dt_to_prompt,
                        timestamp - helper.prompt_timestamp,
                        helper.multiplicity-1)
                assign_value(fill_buf.dr_to_prompt, math.sqrt(
                        (clusters.x-helper.prompt_x)**2 +
                        (clusters.y-helper.prompt_y)**2 +
                        (clusters.z-helper.prompt_z)**2),
                        helper.multiplicity-1)
                assign_value(fill_buf.loopIndex, clusters.loopIndex,
                        helper.multiplicity - 1)
                clusters_index += 1
                helper.last_ADevent_timestamp = clusters.timestamp
            else:  # The coincidence window has ended
                # Fill the TTree and reset, without incrementing the loop index
                if helper.multiplicity == 2:
                    assign_value(fill_buf.multiplicity, helper.multiplicity)
                    outdata.Fill()
                helper.in_coincidence_window = False
                helper.multiplicity = 0
        else:  # i.e. if not isValid
            if helper.in_coincidence_window:
                if timestamp - helper.prompt_timestamp >= COINCIDENCE_WINDOW:
                    # We've left the window
                    if helper.multiplicity == 2:
                        assign_value(fill_buf.multiplicity,
                                helper.multiplicity)
                        outdata.Fill()
                    helper.in_coincidence_window = False
                    helper.multiplicity = 0
            clusters_index += 1
    logging.debug(clusters_index)
    tracker.muon_counts = muon_helper.muon_counts
    #muon_helper.close(10*limit)
    return tracker, muon_helper
    #return tracker

def copy_to_buffer(ttree, buf, index, name):
    assign_value(getattr(buf, name), getattr(ttree, name), index)
    return

def assign_muons(buf, muon_helper):
    dts = muon_helper.dt_previous_showermuons()
    logging.debug('in assign_muons')
    logging.debug(dts)
    try:
        for (energy, tag), value in dts.items():
            name = f'dt_{energy}energy_muon{tag}'
            assign_value(getattr(buf, name), value)
    except:
        logging.debug(value)
        raise

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
    copy_to_buffer(source, buf, index, 'fID')
    copy_to_buffer(source, buf, index, 'fPSD')

def main(entries, events_filename, muon_filename, out_filename, runfile,
        detector, debug):
    from ROOT import TFile

    run, fileno = runfile
    events_file = TFile(events_filename, 'READ')
    muon_file = TFile(muon_filename, 'READ')
    in_events = events_file.Get('events')
    muons = muon_file.Get('muons')
    outfile = TFile(out_filename, 'RECREATE')
    ttree_name = 'ad_events'
    ttree_description = 'AD events by Sam Kohn (git: %s)'
    outdata, fill_buf = create_computed_TTree(ttree_name, outfile,
        ttree_description)
    out_muons, muon_fill_buf = muon_TTree(outfile)

    if entries == -1:
        entries = in_events.GetEntries()
    tracker, muon_helper = main_loop(in_events, muons, outdata, fill_buf, debug, entries,
            out_muons, muon_fill_buf)
    muon_helper.close(100*entries, output=out_filename+'_muonrates.root')
    outfile.Write()
    outfile.Close()
    events_file.Close()
    muon_file.Close()
    with open(os.path.splitext(out_filename)[0] + '.json', 'w') as f:
        json.dump(tracker.export(), f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--adevents', help='Input AD events file name')
    parser.add_argument('--muons', help='Input muons file name')
    parser.add_argument('-o', '--output', help='Output file name')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--events', type=int, default=-1)
    parser.add_argument('-r', '--runfile', type=int, nargs=2,
        help='<run> <fileno>')
    parser.add_argument('--det', type=int, help='AD number (1-4)')
    parser.add_argument('--lowenergy', action='store_true',
            help='Lower energy threshold')
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    if args.lowenergy:
        isADEvent_THU = isADEvent_THU_lowenergy
    main(args.events, args.adevents, args.muons, args.output, args.runfile,
            args.det, args.debug)
