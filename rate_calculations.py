'''
Compute the livetime inefficiency due to muons.

'''
from __future__ import print_function
import logging
import argparse
from math import exp

from ROOT import TFile, TTree

from translate import fetch_value
import muons

class RateHelper(object):
    def __init__(self, run, fileno):
        self.run = run
        self.fileno = fileno
        self.total_nonvetoed_livetime = {1:0, 2:0}  # by AD
        self.start_time = 0
        self.end_time = 0
        self.next_livetime_start = {1:0, 2:0}
        self.number_IBD_candidates = {1:0, 2:0}
        self.number_prompts = {1:0, 2:0}
        self.number_delayeds = {1:0, 2:0}
        self.events_since_last_valid_time = 0

    def compute_results(self):
        NS_PER_DAY = 1e9*60*60*24
        S_PER_DAY = 60*60*24
        daq_livetime = self.end_time - self.start_time
        mu_eff = {n: float(nonvetoed)/daq_livetime for n, nonvetoed in
                self.total_nonvetoed_livetime.items()}
        prompt_rate_Hz = {n: (1e9*num)/self.total_nonvetoed_livetime[n] for n,
                num in self.number_prompts.items()}
        delayed_rate_Hz = {n: (1e9*num)/self.total_nonvetoed_livetime[n] for n,
                num in self.number_delayeds.items()}
        livetime_days = {n: t/NS_PER_DAY for n, t in
                self.total_nonvetoed_livetime.items()}
        ibd_rate_perday = {n: num/livetime_days[n] for n, num in
                self.number_IBD_candidates.items()}
        dt_mult_cut_sec = 200e-6
        mult_eff = {n: self.multiplicity_eff(
            prompt_rate_Hz[n], delayed_rate_Hz[n], dt_mult_cut_sec)
            for n in prompt_rate_Hz}
        acc_rate_perday = {n: S_PER_DAY*self.accidental_rate(
            prompt_rate_Hz[n], delayed_rate_Hz[n], dt_mult_cut_sec)
            for n in prompt_rate_Hz}
        return {
                'run': self.run,
                'fileno': self.fileno,
                'daq_livetime': daq_livetime,
                'usable_livetime': self.total_nonvetoed_livetime,
                'start_time': self.start_time,
                'end_time': self.end_time,
                'muon_efficiency': mu_eff,
                'prompt_rate_Hz': prompt_rate_Hz,
                'delayed_rate_Hz': delayed_rate_Hz,
                'ibd_rate_perday': ibd_rate_perday,
                'multiplicity_efficiency': mult_eff,
                'accidental_rate_perday': acc_rate_perday,
                'number_prompts': self.number_prompts,
                'number_delayeds': self.number_delayeds,
                'number_IBDs': self.number_IBD_candidates,
                }

    @staticmethod
    def multiplicity_eff(rp, rd, dt):
        return exp(-2*rp*dt - rd*dt)

    @staticmethod
    def accidental_rate(rp, rd, dt):
        return rd*rp*dt*exp(-2*rp*dt - rd*dt)



def main(num_events, start_event, debug):
    filename = 'out.root'
    infile = TFile(filename, 'READ')
    indata = infile.Get('data')
    incomputed = infile.Get('computed')

    helper = RateHelper()

    total_entries = min(indata.GetEntries(), incomputed.GetEntries())
    entries = min(num_events+start_event, total_entries) if num_events > 0 else total_entries
    for event_number in xrange(start_event, entries):
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)
        incomputed.LoadTree(event_number)
        incomputed.GetEntry(event_number)
        logging.debug('Event %d', event_number)
        data_list = fetch_data(indata, incomputed)
        one_iteration(event_number, data_list, helper, start_event, entries)
    print_results(helper)

def print_results(helper):
    print('total DAQ livetime:')
    daq_livetime = helper.end_time - helper.start_time
    print(daq_livetime)
    print('total nonvetoed livetime:')
    print(helper.total_nonvetoed_livetime)
    efficiency = {n: float(nonvetoed)/daq_livetime for n, nonvetoed in
            helper.total_nonvetoed_livetime.items()}
    print('efficiency:')
    print(efficiency)
    nonvetoed_livetime_days = {n: t/(1e9 * 60 * 60 * 24) for n, t in
            helper.total_nonvetoed_livetime.items()}
    print('IBD candidates:')
    print(helper.number_IBD_candidates)
    print('delayed-like:')
    print(helper.number_delayeds)
    print('prompt-like:')
    print(helper.number_prompts)
    print('IBD rate per day:')
    print({n: num/nonvetoed_livetime_days[n] for n, num in
        helper.number_IBD_candidates.items()})
    print('delayed-like rate (Hz):')
    print({n: float(num)/helper.total_nonvetoed_livetime[n]*1e9 for n, num in
        helper.number_delayeds.items()})
    print('prompt-like rate (Hz):')
    print({n: float(num)/helper.total_nonvetoed_livetime[n]*1e9 for n, num in
        helper.number_prompts.items()})

def callback_adapter(helper, start_event, entries):
    def callback(event_cache):
        data_list = []
        data_list.append(event_cache.noTree_timestamp[0])
        data_list.append(event_cache.noTree_detector[0])
        data_list.append(event_cache.tag_WSMuon[0])
        data_list.append(event_cache.tag_ADMuon[0])
        data_list.append(event_cache.tag_ShowerMuon[0])
        data_list.append(event_cache.tag_IBDDelayed[0])
        data_list.append(event_cache.tag_PromptLike[0])
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
        one_iteration(event_number, data_list, helper, start_event, entries)
    return callback

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


        if event_number == start_event:
            helper.start_time = timestamp
            for key in helper.next_livetime_start:
                helper.next_livetime_start[key] = (timestamp +
                        muons._SHOWER_MUON_VETO_LAST_NS)
        if event_number == entries - 1:
            helper.end_time = timestamp

        if isIBDDelayed:
            helper.number_IBD_candidates[detector] += 1
        if (isDelayedLike
                and not isFlasher
                and not isWSMuonVetoed
                and not isADMuonVetoed
                and not isShowerMuonVetoed):
            helper.number_delayeds[detector] += 1
        EXTRA_PROMPT_DT = 400e3
        if (isPromptLike
                and not isFlasher
                and dt_next_WSMuon > muons._WSMUON_VETO_NEXT_NS
                and dt_last_WSMuon > (muons._WSMUON_VETO_LAST_NS
                    - EXTRA_PROMPT_DT)
                and dt_last_ADMuon > (muons._ADMUON_VETO_LAST_NS
                    - EXTRA_PROMPT_DT)
                and dt_last_ShowerMuon >
                    (muons._SHOWER_MUON_VETO_LAST_NS
                        - EXTRA_PROMPT_DT)):
            helper.number_prompts[detector] += 1

        logging.debug('next livetime start: %s', helper.next_livetime_start)
        logging.debug('timestamp: %d', timestamp)
        if isWSMuon:
            for det, next_start in helper.next_livetime_start.items():
                new_livetime = (timestamp - muons._WSMUON_VETO_NEXT_NS
                        - next_start)
                if new_livetime > 0:
                    helper.total_nonvetoed_livetime[det] += new_livetime
                new_start = timestamp + muons._WSMUON_VETO_LAST_NS
                if new_start > next_start:
                    helper.next_livetime_start[det] = new_start
        if isADMuon:
            next_start = helper.next_livetime_start[detector]
            new_livetime = timestamp - next_start
            if new_livetime > 0:
                helper.total_nonvetoed_livetime[detector] += new_livetime
            new_start = timestamp + muons._ADMUON_VETO_LAST_NS
            if new_start > next_start:
                helper.next_livetime_start[detector] = new_start
        if isShowerMuon:
            next_start = helper.next_livetime_start[detector]
            new_livetime = timestamp - next_start
            if new_livetime > 0:
                helper.total_nonvetoed_livetime[detector] += new_livetime
            new_start = timestamp + muons._SHOWER_MUON_VETO_LAST_NS
            if new_start > next_start:
                helper.next_livetime_start[detector] = new_start
        logging.debug('new next livetime start: %s', helper.next_livetime_start)
        logging.debug('new total: %s', helper.total_nonvetoed_livetime)

def fetch_data(indata, incomputed):
    timestamp = fetch_value(indata, 'timeStamp', int)
    detector = fetch_value(indata, 'detector', int)
    isWSMuon = fetch_value(incomputed, 'tag_WSMuon', bool)
    isADMuon = fetch_value(incomputed, 'tag_ADMuon', bool)
    isShowerMuon = fetch_value(incomputed, 'tag_ShowerMuon', bool)
    isIBDDelayed = fetch_value(incomputed, 'tag_IBDDelayed', bool)
    isPromptLike = fetch_value(incomputed, 'tag_PromptLike', bool)
    isDelayedLike = fetch_value(incomputed, 'tag_DelayedLike', bool)
    isWSMuonVetoed = fetch_value(incomputed, 'tag_WSMuonVeto', bool)
    isADMuonVetoed = fetch_value(incomputed, 'tag_ADMuonVeto', bool)
    isShowerMuonVetoed = fetch_value(incomputed,
            'tag_ShowerMuonVeto', bool)
    isFlasher = fetch_value(incomputed, 'tag_flasher', bool)
    dt_last_WSMuon = fetch_value(incomputed, 'dt_previous_WSMuon', int)
    dt_next_WSMuon = fetch_value(incomputed, 'dt_next_WSMuon', int)
    dt_last_ADMuon = fetch_value(incomputed, 'dt_previous_ADMuon', int)
    dt_last_ShowerMuon = fetch_value(incomputed,
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
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.num_events, args.start_event, args.debug)

