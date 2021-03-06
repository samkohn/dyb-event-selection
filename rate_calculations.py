'''
Compute the livetime inefficiency due to muons.

'''
from __future__ import print_function
import logging
import argparse
import json
from math import exp
from pprint import pprint

from ROOT import TFile, TTree

from translate import fetch_value
import muons
import delayeds

def AD_dict(nADs, default=0):
    return {n: default for n in range(1, nADs+1)}

class RateHelper(object):
    def __init__(self, run, fileno, site):
        nADs = {1: 2, 2: 2, 3: 4}[site]
        self.site = None
        self.run = run
        self.fileno = fileno
        self.total_nonvetoed_livetime = AD_dict(nADs)
        self.start_time = 0
        self.end_time = 0
        self.current_nomuon_livetime_start = AD_dict(nADs)  # timestamp of end of last muon veto
        self.number_IBD_candidates = AD_dict(nADs)
        self.number_prompts = AD_dict(nADs)
        self.number_delayeds = AD_dict(nADs)
        self.number_prompt_singles = AD_dict(nADs)
        self.number_delayed_singles = AD_dict(nADs)
        self.last_promptlike_event_timestamp = AD_dict(nADs)
        self.last_event_maybe_single = AD_dict(nADs, False)
        self.last_event_isDelayedLike = AD_dict(nADs, False)
        self.current_singles_livetime_start = AD_dict(nADs)
        self.singles_livetime = AD_dict(nADs)
        self.events_since_last_valid_time = 0

    def compute_results(self):
        NS_PER_DAY = 1e9*60*60*24
        S_PER_DAY = 60*60*24
        daq_livetime = self.end_time - self.start_time
        mu_eff = {n: float(nonvetoed)/daq_livetime for n, nonvetoed in
                self.total_nonvetoed_livetime.items()}
        if all(livetime > 0 for livetime in
                self.total_nonvetoed_livetime.values()):
            prompt_rate_Hz = {n: (1e9*num)/self.total_nonvetoed_livetime[n] for n,
                    num in self.number_prompts.items()}
            delayed_rate_Hz = {n: (1e9*num)/self.total_nonvetoed_livetime[n] for n,
                    num in self.number_delayeds.items()}
        else:
            prompt_rate_Hz = {n: 0 for n in self.total_nonvetoed_livetime}
            delayed_rate_Hz = {n: 0 for n in self.total_nonvetoed_livetime}
        if all(livetime > 0 for livetime in self.singles_livetime.values()):
            prompt_singles_rate_Hz = {n: ((1e9*num)/self.singles_livetime[n] if num >
                    0 else 0) for n,
                    num in self.number_prompt_singles.items()}
            delayed_singles_rate_Hz = {n: ((1e9*num)/self.singles_livetime[n] if num >
                    0 else 0) for n,
                    num in self.number_delayed_singles.items()}
        else:
            prompt_singles_rate_Hz = {n: 0 for n in self.singles_livetime}
            delayed_singles_rate_Hz = {n: 0 for n in self.singles_livetime}
        livetime_days = {n: t/NS_PER_DAY for n, t in
                self.total_nonvetoed_livetime.items()}
        dt_mult_cut_sec = 200e-6
        mult_eff = {n: self.multiplicity_eff(
                prompt_singles_rate_Hz[n],
                delayed_singles_rate_Hz[n],
                dt_mult_cut_sec)
            for n in prompt_rate_Hz}
        if all(eff > 0 for eff in mu_eff.values()):
            ibd_rate_perday = {n:
                    num/(mu_eff[n]*mult_eff[n]*(daq_livetime/NS_PER_DAY)) for n, num in
                    self.number_IBD_candidates.items()}
        else:
            ibd_rate_perday = {n: 0 for n in self.number_IBD_candidates}
        acc_rate_perday = {n: S_PER_DAY*self.accidental_rate(
                prompt_singles_rate_Hz[n],
                delayed_singles_rate_Hz[n],
                dt_mult_cut_sec)
            for n in prompt_rate_Hz}
        return {
                'run': self.run,
                'fileno': self.fileno,
                'site': self.site,
                'daq_livetime': daq_livetime,
                'usable_livetime': self.total_nonvetoed_livetime,
                'singles_livetime': self.singles_livetime,
                'start_time': self.start_time,
                'end_time': self.end_time,
                'muon_efficiency': mu_eff,
                'prompt_rate_Hz': prompt_rate_Hz,
                'delayed_rate_Hz': delayed_rate_Hz,
                'prompt_singles_rate_Hz': prompt_singles_rate_Hz,
                'delayed_singles_rate_Hz': delayed_singles_rate_Hz,
                'ibd_rate_perday': ibd_rate_perday,
                'multiplicity_efficiency': mult_eff,
                'accidental_rate_perday': acc_rate_perday,
                'number_prompts': self.number_prompts,
                'number_delayeds': self.number_delayeds,
                'number_prompt_singles': self.number_prompt_singles,
                'number_delayed_singles': self.number_delayed_singles,
                'number_IBDs': self.number_IBD_candidates,
                }

    @staticmethod
    def multiplicity_eff(rp, rd, dt):
        return exp(-2*rp*dt - rd*dt)

    @staticmethod
    def accidental_rate(rp, rd, dt):
        return rd*rp*dt*exp(-2*rp*dt - rd*dt)



def main(filename, site, run, fileno, num_events, start_event, debug,
        selection_name):
    infile = TFile(filename, 'READ')
    indata = infile.Get('data')

    helper = RateHelper(run, fileno, site)

    total_entries = indata.GetEntries()
    entries = min(num_events+start_event, total_entries) if num_events > 0 else total_entries
    for event_number in xrange(start_event, entries):
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)
        logging.debug('Event %d', event_number)
        data_list = fetch_data(indata)
        if selection_name == 'nh_THU':
            one_iteration_nh_THU(event_number, data_list, helper, start_event,
                    entries)
        else:
            one_iteration(event_number, data_list, helper, start_event, entries)
    print_results(helper)
    outname = filename.split('.')[0] + '.json'
    with open(outname, 'w') as f:
        json.dump(helper.compute_results(), f)


def print_results(helper):
    pprint(helper.compute_results())

def callback_adapter(helper, start_event, entries, selection_name):
    def callback(event_cache):
        is_crosstrigger = (event_cache.triggerType[0] & 0x2) != 0
        data_list = []
        data_list.append(event_cache.timestamp[0])
        data_list.append(event_cache.detector[0])
        data_list.append(event_cache.tag_WSMuon[0] if not is_crosstrigger else
                0)
        data_list.append(event_cache.tag_ADMuon[0] if not is_crosstrigger else
                0)
        data_list.append(event_cache.tag_ShowerMuon[0] if not is_crosstrigger
                else 0)
        data_list.append(event_cache.tag_PromptLike[0])
        if selection_name != 'nh_THU':
            data_list.append(event_cache.tag_IBDDelayed[0])
            data_list.append(event_cache.tag_DelayedLike[0])
        data_list.append(event_cache.tag_WSMuonVeto[0])
        data_list.append(event_cache.tag_ADMuonVeto[0])
        data_list.append(event_cache.tag_ShowerMuonVeto[0])
        data_list.append(event_cache.tag_flasher[0])
        data_list.append(event_cache.dt_previous_WSMuon[0])
        data_list.append(event_cache.dt_next_WSMuon[0])
        data_list.append(event_cache.dt_previous_ADMuon[0])
        data_list.append(event_cache.dt_previous_ShowerMuon[0])
        event_number = event_cache.noTree_loopIndex[0]
        if selection_name == 'nh_THU':
            one_iteration_nh_THU(event_number, data_list, helper, start_event,
                    entries)
        else:
            one_iteration(event_number, data_list, helper, start_event, entries)
    return callback

def one_iteration_nh_THU(event_number, data_list, helper, start_event,
        entries):
        (
            timestamp,
            detector,
            isWSMuon,
            isADMuon,
            isShowerMuon,
            isPromptLike,
            isWSMuonVetoed,
            isADMuonVetoed,
            isShowerMuonVetoed,
            isFlasher,
            dt_last_WSMuon,
            dt_next_WSMuon,
            dt_last_ADMuon,
            dt_last_ShowerMuon
            ) = data_list

        muon_type = 'None'
        if isWSMuon:
            muon_type = 'WS'
        if isADMuon:
            muon_type = 'AD'
        if isShowerMuon:
            muon_type = 'Shower'
        logging.debug('muon type: %s', muon_type)

        logging.debug('end of most recent muon vetoes, before this event: %s',
                helper.current_nomuon_livetime_start)
        logging.debug('timestamp: %d', timestamp)

        if event_number == start_event:
            helper.start_time = timestamp
            for key in helper.current_nomuon_livetime_start:
                helper.current_nomuon_livetime_start[key] = timestamp
        if event_number == entries - 1:
            helper.end_time = timestamp
            if not (isWSMuon or isADMuon or isShowerMuon):
                for key, window_start_time in (
                        helper.current_nomuon_livetime_start.items()):
                    new_time = max(0, timestamp - window_start_time)
                    helper.total_nonvetoed_livetime[key] += new_time

        if isWSMuon:
            for det in helper.current_nomuon_livetime_start:
                current_start = helper.current_nomuon_livetime_start[det]
                if current_start < timestamp - delayeds._NH_THU_MULT_POST_MIN:
                    helper.total_nonvetoed_livetime[det] += (timestamp
                            - current_start - delayeds._NH_THU_MULT_POST_MIN)
                    helper.current_nomuon_livetime_start[det] = (timestamp
                            + muons._NH_WSMUON_VETO_LAST_NS)
                else:
                    helper.current_nomuon_livetime_start[det] = max(
                            current_start, timestamp + muons._NH_WSMUON_VETO_LAST_NS)
        if isADMuon:
            current_start = helper.current_nomuon_livetime_start[detector]
            if current_start < timestamp - delayeds._NH_THU_MULT_POST_MIN:
                helper.total_nonvetoed_livetime[detector] += (timestamp
                        - current_start - delayeds._NH_THU_MULT_POST_MIN)
                helper.current_nomuon_livetime_start[detector] = (timestamp
                        + muons._NH_ADMUON_VETO_LAST_NS)
            else:
                helper.current_nomuon_livetime_start[detector] = max(
                        current_start, timestamp + muons._NH_ADMUON_VETO_LAST_NS)
        if isShowerMuon:
            current_start = helper.current_nomuon_livetime_start[detector]
            if current_start < timestamp - delayeds._NH_THU_MULT_POST_MIN:
                helper.total_nonvetoed_livetime[detector] += (timestamp
                        - current_start - delayeds._NH_THU_MULT_POST_MIN)
                helper.current_nomuon_livetime_start[detector] = (timestamp
                        + muons._NH_SHOWER_MUON_VETO_LAST_NS)
            else:
                helper.current_nomuon_livetime_start[detector] = max(
                        current_start, timestamp
                        + muons._NH_SHOWER_MUON_VETO_LAST_NS)


def one_iteration(event_number, data_list, helper, start_event, entries):
        (
            timestamp,
            detector,
            isWSMuon,
            isADMuon,
            isShowerMuon,
            isIBDDelayed,
            isPromptLike,
            isDelayedLike,
            isWSMuonVetoed,
            isADMuonVetoed,
            isShowerMuonVetoed,
            isFlasher,
            dt_last_WSMuon,
            dt_next_WSMuon,
            dt_last_ADMuon,
            dt_last_ShowerMuon
            ) = data_list

        muon_type = 'None'
        if isWSMuon:
            muon_type = 'WS'
        if isADMuon:
            muon_type = 'AD'
        if isShowerMuon:
            muon_type = 'Shower'
        logging.debug('muon type: %s', muon_type)

        logging.debug('start of singles livetime windows, before this event: %s',
                helper.current_singles_livetime_start)
        logging.debug('end of most recent muon vetoes, before this event: %s',
                helper.current_nomuon_livetime_start)
        logging.debug('timestamp: %d', timestamp)

        if event_number == start_event:
            helper.start_time = timestamp
            for key in helper.current_nomuon_livetime_start:
                helper.current_nomuon_livetime_start[key] = timestamp
            for key in helper.current_singles_livetime_start:
                helper.current_singles_livetime_start[key] = timestamp
        if event_number == entries - 1:
            helper.end_time = timestamp
            if not (isWSMuon or isADMuon or isShowerMuon):
                for key, window_start_time in (
                        helper.current_nomuon_livetime_start.items()):
                    new_time = max(0, timestamp - window_start_time)
                    helper.total_nonvetoed_livetime[key] += new_time


        if isIBDDelayed:
            helper.number_IBD_candidates[detector] += 1
        if (isDelayedLike
                and not isFlasher
                and not isWSMuonVetoed
                and not isADMuonVetoed
                and not isShowerMuonVetoed):
            helper.number_delayeds[detector] += 1
        if (isPromptLike
                and not isFlasher
                and not isWSMuonVetoed
                and not isADMuonVetoed
                and not isShowerMuonVetoed):
            helper.number_prompts[detector] += 1
            # Count uncorrelated singles (uncorrelated at the timescale of
            # 200us)
            singles_dt_ns = int(200e3)
            logging.debug('promptlike. Timestamp: %d', timestamp)
            # Check for muon veto
            if timestamp > helper.current_nomuon_livetime_start[detector]:
                logging.debug('not muon-vetoed')
                # Check if this event is far enough from last event for either
                # to potentially be single events
                if (timestamp - singles_dt_ns >
                        helper.last_promptlike_event_timestamp[detector]):
                    logging.debug('far enough from last event to maybe be a single')
                    # Check if last event was also isolated from
                    # next-to-last event (if so, it was a single)
                    #
                    # Note: this ignores the possibility that a muon event
                    # occurs between this event and the last single-like event,
                    # but I believe that effect is negligible.
                    if helper.last_event_maybe_single[detector]:
                        logging.debug('last single-like event was a single')
                        helper.number_prompt_singles[detector] += 1
                        if helper.last_event_isDelayedLike[detector]:
                            helper.number_delayed_singles[detector] += 1
                    # Flag current event as potentially a single event
                    helper.last_event_maybe_single[detector] = True
                    helper.last_event_isDelayedLike[detector] = isDelayedLike
                else:
                    logging.debug('not far enough from last event to be a single')
                    # Add to the total singles livetime all the time from the
                    # start of the current window up until singles_dt before
                    # the last event.
                    new_time = (helper.last_promptlike_event_timestamp[detector]
                            - singles_dt_ns
                            - helper.current_singles_livetime_start[detector])
                    logging.debug('new singles livetime = max(%d, 0)',
                            new_time)
                    helper.singles_livetime[detector] += max(new_time, 0)
                    # Reset the current singles livetime start to singles_dt
                    # past the current event
                    helper.current_singles_livetime_start[detector] = (timestamp
                            + singles_dt_ns)
                    helper.last_event_maybe_single[detector] = False
                # Save this event's timestamp whether it might be a single or
                # not
                helper.last_promptlike_event_timestamp[detector] = timestamp

        # Deal with intersection between singles and muons
        if isWSMuon or isADMuon or isShowerMuon:
            window_end = timestamp
            if isWSMuon:
                window_end -= muons._WSMUON_VETO_NEXT_NS
                detectors = helper.current_nomuon_livetime_start.keys()
                veto_length = muons._WSMUON_VETO_LAST_NS
            else:
                detectors = (detector,)
                if isADMuon:
                    veto_length = muons._ADMUON_VETO_LAST_NS
                elif isShowerMuon:
                    veto_length = muons._SHOWER_MUON_VETO_LAST_NS
            for det in detectors:
                # Add to singles livetime the time between the last veto and the
                # current muon.
                if window_end > helper.current_nomuon_livetime_start[det]:
                    window_start = helper.current_singles_livetime_start[det]
                    logging.debug('detector %d, window_start: %d', det,
                            window_start)
                    dt = window_end - window_start
                    if dt > 0:
                        logging.debug('new time: %d', dt)
                        helper.singles_livetime[det] += dt
                # Update singles livetime window start to be the end of the muon
                # veto window
                veto_end = timestamp + veto_length
                logging.debug('end of veto window: %d', veto_end)
                if (veto_end > helper.current_singles_livetime_start[det]):
                    helper.current_singles_livetime_start[det] = veto_end

        if isWSMuon:
            for det, window_start_time in helper.current_nomuon_livetime_start.items():
                new_livetime = (timestamp - muons._WSMUON_VETO_NEXT_NS
                        - window_start_time)
                if new_livetime > 0:
                    helper.total_nonvetoed_livetime[det] += new_livetime
                new_start = timestamp + muons._WSMUON_VETO_LAST_NS
                if new_start > window_start_time:
                    helper.current_nomuon_livetime_start[det] = new_start
        if isADMuon:
            window_start_time = helper.current_nomuon_livetime_start[detector]
            new_livetime = timestamp - window_start_time
            if new_livetime > 0:
                helper.total_nonvetoed_livetime[detector] += new_livetime
            new_start = timestamp + muons._ADMUON_VETO_LAST_NS
            if new_start > window_start_time:
                helper.current_nomuon_livetime_start[detector] = new_start
        if isShowerMuon:
            window_start_time = helper.current_nomuon_livetime_start[detector]
            new_livetime = timestamp - window_start_time
            if new_livetime > 0:
                helper.total_nonvetoed_livetime[detector] += new_livetime
            new_start = timestamp + muons._SHOWER_MUON_VETO_LAST_NS
            if new_start > window_start_time:
                helper.current_nomuon_livetime_start[detector] = new_start
        logging.debug('start of singles livetime windows, after this event: %s',
                helper.current_singles_livetime_start)
        logging.debug('end of most recent muon vetoes, after this event: %s',
                helper.current_nomuon_livetime_start)
        logging.debug('new non-vetoed livetime total: %s', helper.total_nonvetoed_livetime)

def fetch_data(indata):
    timestamp = fetch_value(indata, 'timestamp', int)
    detector = fetch_value(indata, 'detector', int)
    isWSMuon = fetch_value(indata, 'tag_WSMuon', bool)
    isADMuon = fetch_value(indata, 'tag_ADMuon', bool)
    isShowerMuon = fetch_value(indata, 'tag_ShowerMuon', bool)
    isIBDDelayed = fetch_value(indata, 'tag_IBDDelayed', bool)
    isPromptLike = fetch_value(indata, 'tag_PromptLike', bool)
    isDelayedLike = fetch_value(indata, 'tag_DelayedLike', bool)
    isWSMuonVetoed = fetch_value(indata, 'tag_WSMuonVeto', bool)
    isADMuonVetoed = fetch_value(indata, 'tag_ADMuonVeto', bool)
    isShowerMuonVetoed = fetch_value(indata,
            'tag_ShowerMuonVeto', bool)
    isFlasher = fetch_value(indata, 'tag_flasher', bool)
    dt_last_WSMuon = fetch_value(indata, 'dt_previous_WSMuon', int)
    dt_next_WSMuon = fetch_value(indata, 'dt_next_WSMuon', int)
    dt_last_ADMuon = fetch_value(indata, 'dt_previous_ADMuon', int)
    dt_last_ShowerMuon = fetch_value(indata,
            'dt_previous_ShowerMuon', int)
    return (
            timestamp,
            detector,
            isWSMuon,
            isADMuon,
            isShowerMuon,
            isIBDDelayed,
            isPromptLike,
            isDelayedLike,
            isWSMuonVetoed,
            isADMuonVetoed,
            isShowerMuonVetoed,
            isFlasher,
            dt_last_WSMuon,
            dt_next_WSMuon,
            dt_last_ADMuon,
            dt_last_ShowerMuon
            )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--num-events', type=int, default=-1)
    parser.add_argument('-s', '--start-event', type=int, default=0)
    parser.add_argument('-f', '--filename', default='out.root')
    parser.add_argument('--run', default=0, type=int)
    parser.add_argument('--fileno', default=0, type=int)
    parser.add_argument('--site', type=int)
    parser.add_argument('--selection', default='ngd')
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.filename, args.site, args.run, args.fileno, args.num_events,
            args.start_event, args.debug, args.selection)

