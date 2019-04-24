'''
Compute the livetime inefficiency due to muons.

'''
from __future__ import print_function
import logging
import argparse

from ROOT import TFile, TTree

from translate import fetch_value
import muons

def main(num_events, debug):
    filename = 'out.root'
    infile = TFile(filename, 'READ')
    indata = infile.Get('data')
    incomputed = infile.Get('computed')

    total_nonvetoed_livetime = {1:0, 2:0}  # by AD
    start_time = 0
    end_time = 0
    number_IBD_candidates = {1:0, 2:0}
    number_prompts = {1:0, 2:0}
    number_delayeds = {1:0, 2:0}

    entries = num_events if num_events > 0 else indata.GetEntries()
    for event_number in xrange(entries):
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)
        incomputed.LoadTree(event_number)
        incomputed.GetEntry(event_number)

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


        if event_number == 0:
            start_time = timestamp
        if event_number == entries - 1:
            end_time = timestamp

        if isIBDDelayed:
            number_IBD_candidates[detector] += 1
        if (isDelayedLike
                and not isFlasher
                and not isWSMuonVetoed
                and not isADMuonVetoed
                and not isShowerMuonVetoed):
            number_delayeds[detector] += 1
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
            number_prompts[detector] += 1

        if isWSMuon or isADMuon or isShowerMuon:
            dt_past_WSVeto = max(0, dt_last_WSMuon -
                    muons._WSMUON_VETO_LAST_NS)
            dt_past_ADVeto = max(0, dt_last_ADMuon -
                    muons._ADMUON_VETO_LAST_NS)
            dt_past_ShowerVeto = max(0, dt_last_ShowerMuon -
                    muons._SHOWER_MUON_VETO_LAST_NS)
            logging.debug('last WSMuon: %d', dt_last_WSMuon)
            logging.debug('dt WS Veto: %d', dt_past_WSVeto)
            logging.debug('last_ADMuon: %d', dt_last_ADMuon)
            logging.debug('dt AD Veto: %d', dt_past_ADVeto)
            logging.debug('last shower: %d', dt_last_ShowerMuon)
            logging.debug('dt shower veto: %d', dt_past_ShowerVeto)

            new_livetime = min(dt_past_WSVeto, dt_past_ADVeto,
                    dt_past_ShowerVeto)

            if isWSMuon:
                new_livetime -= muons._WSMUON_VETO_NEXT_NS
                new_livetime = max(0, new_livetime)
                for det in total_nonvetoed_livetime:
                    total_nonvetoed_livetime[det] += new_livetime
            else:
                total_nonvetoed_livetime[detector] += new_livetime
            logging.debug('new livetime: %d', new_livetime)
    print('total DAQ livetime:')
    daq_livetime = end_time - start_time
    print(daq_livetime)
    print('total nonvetoed livetime:')
    print(total_nonvetoed_livetime)
    efficiency = {n: float(nonvetoed)/daq_livetime for n, nonvetoed in
            total_nonvetoed_livetime.items()}
    print('efficiency:')
    print(efficiency)
    nonvetoed_livetime_days = {n: t/(1e9 * 60 * 60 * 24) for n, t in
            total_nonvetoed_livetime.items()}
    print('IBD candidates:')
    print(number_IBD_candidates)
    print('delayed-like:')
    print(number_delayeds)
    print('prompt-like:')
    print(number_prompts)
    print('IBD rate per day:')
    print({n: num/nonvetoed_livetime_days[n] for n, num in
        number_IBD_candidates.items()})
    print('delayed-like rate (Hz):')
    print({n: float(num)/total_nonvetoed_livetime[n]*1e9 for n, num in
        number_delayeds.items()})
    print('prompt-like rate (Hz):')
    print({n: float(num)/total_nonvetoed_livetime[n]*1e9 for n, num in
        number_prompts.items()})

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--num-events', type=int, default=-1)
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.num_events, args.debug)

