'''
Aggregate the output JSON files to get total livetime, accidentals rate, etc.

'''

from __future__ import print_function
import argparse
import os
import json

from rate_calculations import RateHelper

def main(run, files, site):
    helper = RateHelper(run, 'all', site)
    daq_livetime = 0
    for filename in files:
        with open(filename, 'r') as f:
            results = json.load(f)
        daq_livetime += results['daq_livetime']
        for detector_str, livetime in results['usable_livetime'].items():
            detector = int(detector_str)
            helper.total_nonvetoed_livetime[detector] += livetime
        for detector_str, livetime in results['singles_livetime'].items():
            detector = int(detector_str)
            helper.singles_livetime[detector] += livetime
        for detector_str, number_prompts in results['number_prompts'].items():
            detector = int(detector_str)
            helper.number_prompts[detector] += number_prompts
        for detector_str, number_delayeds in (
                results['number_delayeds'].items()):
            detector = int(detector_str)
            helper.number_delayeds[detector] += number_delayeds
        for detector_str, number_prompts in (
                results['number_prompt_singles'].items()):
            detector = int(detector_str)
            helper.number_prompt_singles[detector] += number_prompts
        for detector_str, number_delayeds in (
                results['number_delayed_singles'].items()):
            detector = int(detector_str)
            helper.number_delayed_singles[detector] += number_delayeds
        for detector_str, number_IBDs in results['number_IBDs'].items():
            detector = int(detector_str)
            helper.number_IBD_candidates[detector] += number_IBDs
    helper.start_time = 0
    helper.end_time = daq_livetime
    aggregated = helper.compute_results()
    print(aggregated)
    with open('out_%d_all.json' % run, 'w') as f:
        json.dump(aggregated, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    parser.add_argument('-r', '--run', type=int)
    parser.add_argument('--site', type=int)
    args = parser.parse_args()
    main(args.run, args.files, args.site)

