'''
Process the basic data by computing various cuts, quantities and tags.

'''
import argparse
import itertools
import logging
import math
import os
import json

from common import *
import muons
from muons import (
        _NH_WSMUON_VETO_LAST_NS as WSMUON_VETO,
        _NH_ADMUON_VETO_LAST_NS as ADMUON_VETO,
        _NH_SHOWER_MUON_VETO_LAST_NS as SHOWER_MUON_VETO
)
import delayeds
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
    def __init__(self, muon_ttree, time_tracker, ad_num):
        self.ttree = muon_ttree
        self.ad = ad_num
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
        return

    def neutron_tag_last_muon(self):
        logging.debug('tagging last muon!')
        (label, _), dt = min(self.dt_previous_showermuons(tags=['']).items(),
                key=lambda pair:pair[1])
        if label == 'low':
            if self.time_lowenergy_muon_ntag != self.time_lowenergy_muon:
                # only increment counter for new tags (don't want to
                # double-count if there are 2 neutrons)
                self.muon_counts['low_ntag'] += 1
            self.time_lowenergy_muon_ntag = self.time_lowenergy_muon
            logging.debug('tagging last low energy muon!')
        elif label == 'mid':
            if self.time_midenergy_muon_ntag != self.time_midenergy_muon:
                # only increment counter for new tags (don't want to
                # double-count if there are 2 neutrons)
                self.muon_counts['mid_ntag'] += 1
            self.time_midenergy_muon_ntag = self.time_midenergy_muon
            logging.debug('tagging last mid energy muon!')
        elif label == 'high':
            if self.time_highenergy_muon_ntag != self.time_highenergy_muon:
                # only increment counter for new tags (don't want to
                # double-count if there are 2 neutrons)
                self.muon_counts['high_ntag'] += 1
            self.time_highenergy_muon_ntag = self.time_highenergy_muon
            logging.debug('tagging last high energy muon!')
        else:
            raise RuntimeError("Can't find that muon")
        return

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
        return self.isVetoed() or self.dt_next_muon() < COINCIDENCE_WINDOW

    def isVetoed(self):
        # Only veto on WS Muons that aren't accompanied by a shower muon
        return (self.dt_previous_WSMuon() < muons._NH_WSMUON_VETO_LAST_NS
                    and min(self.dt_previous_showermuons().values()) >
                        muons._NH_WSMUON_VETO_LAST_NS)

    def close(self, last_entry=-1):
        if last_entry == -1:
            total_entries = self.ttree.GetEntries()
        else:
            total_entries = min(last_entry, self.ttree.GetEntries())
        mu_data = self.ttree
        while self._muon_entry < total_entries:
            mu_data.GetEntry(self._muon_entry)
            is_WS = muons.isWSMuon_nH(mu_data.detector, mu_data.nHit,
                    mu_data.triggerType)
            if is_WS:
                log_WSMuon(self.time_tracker, mu_data.timestamp)
            self._muon_entry += 1

    def load(self, timestamp):
        self._event_timestamp = timestamp
        mu_data = self.ttree
        if timestamp < mu_data.timestamp:
            if self.time_previous_WSMuon == 0:
                # This is during the first few events, when there may be an AD
                # event before the first muon event.
                pass
            else:
                raise RuntimeError('cluster events are out of order')
        mu_data.GetEntry(self._muon_entry)
        total_entries = mu_data.GetEntries()
        # Find the dts to previous muons
        while timestamp > mu_data.timestamp and self._muon_entry < total_entries:
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
            if is_WS:
                self.time_WSMuon = mu_data.timestamp
                log_WSMuon(self.time_tracker, mu_data.timestamp)
            elif is_lowenergy:
                self.muon_counts['low'] += 1
                logging.debug('new low energy muon')
                self.time_lowenergy_muon = mu_data.timestamp
            elif is_midenergy:
                self.muon_counts['mid'] += 1
                logging.debug('new mid energy muon')
                self.time_midenergy_muon = mu_data.timestamp
            elif is_highenergy:
                self.muon_counts['high'] += 1
                logging.debug('new high energy muon')
                self.time_highenergy_muon = mu_data.timestamp
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
        # Save the muon entry number to revert back to after the search forward
        saved_entry = max(self._muon_entry - 1, 0)
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
        self._muon_entry = saved_entry
        mu_data.GetEntry(self._muon_entry)


def is_neutron_tag(buf, multiplicity):
    if multiplicity != 1:
        # Since the coincidence time is so large,
        # only return false if the second event
        # is closer than 400us to the potential neutron tag.
        if buf.dt_to_prompt[1] < 400e3:
            return False
    energy = buf.energy[0]
    if energy < 1.8 or energy > 12:
        return False
    SOONEST = 20000
    LATEST = 200000
    actual = min(buf.dt_lowenergy_muon[0],
            buf.dt_midenergy_muon[0], buf.dt_highenergy_muon[0])
    logging.debug('neutron tag? actual time since most recent: %d', actual)
    return (actual > SOONEST and actual < LATEST)


def main_loop(clusters, muons, outdata, fill_buf, debug, limit):
    clusters.GetEntry(0)
    muons.GetEntry(0)
    # Start tracking DAQ and veto time after the "phantom muon" veto window
    max_veto_window = max(WSMUON_VETO, ADMUON_VETO, SHOWER_MUON_VETO)
    t0 = min(clusters.timestamp, muons.timestamp) + max_veto_window
    tracker = TimeTracker(t0, COINCIDENCE_WINDOW)
    helper = CoincidenceHelper()
    muon_helper = MuonHelper(muons, tracker, clusters.detector)
    clusters_index = 0
    # Process the initial "phantom muon" veto window without tracking DAQ time
    # or muon veto count, and without processing coincidences.
    while clusters.timestamp < t0 and clusters_index < limit:
        clusters_index += 1
        clusters.GetEntry(clusters_index)
    logging.debug(limit)
    while clusters_index < limit:
        logging.debug(clusters_index)
        clusters.GetEntry(clusters_index)
        logging.debug(clusters.loopIndex)
        isValid, reasons = isValidEvent(clusters, helper)
        timestamp = clusters.timestamp
        if isValid:
            muon_helper.load(timestamp)
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
                # Also determine if this event is a neutron tag
                if is_neutron_tag(fill_buf, helper.multiplicity):
                    muon_helper.neutron_tag_last_muon()
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
                    # Also determine if this event is a neutron tag
                    if is_neutron_tag(fill_buf, helper.multiplicity):
                        muon_helper.neutron_tag_last_muon()
                    if helper.multiplicity == 2:
                        assign_value(fill_buf.multiplicity,
                                helper.multiplicity)
                        outdata.Fill()
                    helper.in_coincidence_window = False
                    helper.multiplicity = 0
            clusters_index += 1
    logging.debug(clusters_index)
    muon_helper.close(limit * 100)
    tracker.muon_counts = muon_helper.muon_counts
    return tracker

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

    if entries == -1:
        entries = in_events.GetEntries()
    tracker = main_loop(in_events, muons, outdata, fill_buf, debug, entries)
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
