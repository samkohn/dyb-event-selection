'''
Aggregate the output JSON files to get total livetime, accidentals rate, etc.

'''

from __future__ import print_function
import argparse
import os
import json
import sqlite3

def main2(files, site, outfile):
    daq_livetime = 0
    usable_livetime = 0
    num_veto_windows = 0
    muon_counts = None
    for filename in files:
        with open(filename, 'r') as f:
            results = json.load(f)
        daq_livetime += results['daq_livetime']
        usable_livetime += results['usable_livetime']
        num_veto_windows += results['num_veto_windows']
        if muon_counts is None:
            muon_counts = results['muon_counts']
        else:
            try:
                for key, value in results['muon_counts'].items():
                    muon_counts[key] += value
            except KeyError:
                print(filename)
                raise
    efficiency = usable_livetime / daq_livetime
    rate = num_veto_windows*1e9/usable_livetime
    with open(outfile, 'w') as f:
        json.dump({
            'site': site,
            'daq_livetime': daq_livetime,
            'usable_livetime': usable_livetime,
            'usable_fraction': efficiency,
            'num_veto_windows': num_veto_windows,
            'muon_counts': muon_counts,
            }, f)
    #if db is not None:
        #with sqlite3.Connection(db) as conn:
            #cursor = conn.cursor()
            #cursor.execute('INSERT OR REPLACE INTO muon_rates '
                    #'VALUES (?, ?, ?, ?, ?, ?)',
                    #(run, ad, num_veto_windows, usable_livetime,
                        #rate, efficiency))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    parser.add_argument('--site', type=int)
    parser.add_argument('-o', '--output')
    #parser.add_argument('--update-db',
            #help='Optional, register the run info and muon rate with the given '
            #'database')
    args = parser.parse_args()
    main2(args.files, args.site, args.output)

