'''
Compute the livetime inefficiency due to muons.

'''
from __future__ import print_function
import logging
import argparse
import json
from math import exp

from ROOT import TFile, TTree

from translate import fetch_value
import muons

def AD_dict(nADs, default=0):
    return {n: default for n in range(1, nADs+1)}

class RateHelper(object):
    def __init__(self, run, fileno):
        nADs = 4
        self.site = None
        self.run = run
        self.fileno = fileno
        self.total_nonvetoed_livetime = AD_dict(nADs)
        self.start_time = 0
        self.end_time = 0
        self.next_livetime_start = AD_dict(nADs)
        self.number_IBD_candidates = AD_dict(nADs)
        self.number_prompts = AD_dict(nADs)
        self.number_delayeds = AD_dict(nADs)
        self.number_prompt_singles = AD_dict(nADs)
        self.number_delayed_singles = AD_dict(nADs)
        self.next_singles_livetime_start = AD_dict(nADs)
        self.last_singles_buffer = AD_dict(nADs)
        self.last_single_vetoed = AD_dict(nADs, False)
        self.singles_livetime = AD_dict(nADs)
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
        prompt_singles_rate_Hz = {n: (1e9*num)/self.singles_livetime[n] for n,
                num in self.number_prompt_singles.items()}
        delayed_singles_rate_Hz = {n: (1e9*num)/self.singles_livetime[n] for n,
                num in self.number_delayed_singles.items()}
        livetime_days = {n: t/NS_PER_DAY for n, t in
                self.total_nonvetoed_livetime.items()}
        ibd_rate_perday = {n: num/livetime_days[n] for n, num in
                self.number_IBD_candidates.items()}
        dt_mult_cut_sec = 200e-6
        mult_eff = {n: self.multiplicity_eff(
                prompt_singles_rate_Hz[n],
                delayed_singles_rate_Hz[n],
                dt_mult_cut_sec)
            for n in prompt_rate_Hz}
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



def main(filename, run, fileno, num_events, start_event, debug):
    infile = TFile(filename, 'READ')
    indata = infile.Get('data')

    helper = RateHelper(run, fileno)

    total_entries = indata.GetEntries()
    entries = min(num_events+start_event, total_entries) if num_events > 0 else total_entries
    for event_number in xrange(start_event, entries):
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)
        logging.debug('Event %d', event_number)
        data_list = fetch_data(indata)
        one_iteration(event_number, data_list, helper, start_event, entries)
    print_results(helper)
    outname = filename.split('.')[0] + '.json'
    with open(outname, 'w') as f:
        json.dump(helper.compute_results(), f)


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
    print('delayed-like singles:')
    print(helper.number_delayed_singles)
    print('prompt-like singles:')
    print(helper.number_prompt_singles)
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
        data_list.append(event_cache.timestamp[0])
        data_list.append(event_cache.detector[0])
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
                helper.next_livetime_start[key] = timestamp
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
            buffer_start_time = max(helper.next_livetime_start[detector],
                    helper.last_singles_buffer[detector])
            logging.debug('buffer start time: %d', buffer_start_time)
            if timestamp > helper.next_singles_livetime_start[detector]:
                if helper.last_single_vetoed[detector]:
                    logging.debug('last single was vetoed')
                    helper.last_singles_buffer[detector] = (
                            helper.next_singles_livetime_start[detector])
                    logging.debug('next buffer start time: %d',
                            helper.last_singles_buffer[detector])
                    helper.last_single_vetoed[detector] = False
                else:
                    logging.debug('last single was not vetoed')
                    new_time = timestamp-singles_dt_ns - buffer_start_time
                    helper.singles_livetime[detector] += new_time
                    logging.debug('new time: %d', new_time)
                    helper.number_prompt_singles[detector] += 1
                    if isDelayedLike:
                        helper.number_delayed_singles[detector] += 1
                    helper.last_singles_buffer[detector] = (timestamp
                            - singles_dt_ns)
                    logging.debug('next buffer start time: %d',
                            helper.last_singles_buffer[detector])
            else:
                logging.debug('current event is vetoed')
                helper.last_single_vetoed[detector] = True
                helper.last_singles_buffer[detector] = timestamp+singles_dt_ns
            helper.next_singles_livetime_start[detector] = (timestamp
                    + singles_dt_ns)
            logging.debug('next singles start: %d', timestamp+singles_dt_ns)
        # Deal with intersection between singles and muons
        if isWSMuon or isADMuon or isShowerMuon:
            # Add to singles livetime the time between the last veto and the
            # current muon.
            window_end = timestamp
            if isWSMuon:
                window_end -= muons._WSMUON_VETO_NEXT_NS
                detectors = helper.next_livetime_start.keys()
            else:
                detectors = (detector,)
            for det in detectors:
                window_start = max(helper.last_singles_buffer[det],
                        helper.next_livetime_start[det])
                logging.debug('detector %d, window_start: %d', det,
                        window_start)
                dt = window_end - window_start
                if dt > 0:
                    logging.debug('new time: %d', dt)
                    helper.singles_livetime[det] += dt



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
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.filename, args.run, args.fileno, args.num_events, args.start_event, args.debug)

