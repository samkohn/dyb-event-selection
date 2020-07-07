'''
Aggregate the output JSON files to get total livetime, accidentals rate, etc.

'''

from __future__ import print_function
import argparse
import os
import json
import sqlite3

def main2(run, files, site, ad, outfile, db):
    daq_livetime = 0
    usable_livetime = 0
    num_veto_windows = 0
    for filename in files:
        with open(filename, 'r') as f:
            results = json.load(f)
        daq_livetime += results['daq_livetime']
        usable_livetime += results['usable_livetime']
        num_veto_windows += results['num_veto_windows']
    efficiency = usable_livetime / daq_livetime
    rate = num_veto_windows*1e9/usable_livetime
    with open(outfile, 'w') as f:
        json.dump({
            'run': run,
            'site': site,
            'ad': ad,
            'daq_livetime': daq_livetime,
            'usable_livetime': usable_livetime,
            'usable_fraction': efficiency,
            'num_veto_windows': num_veto_windows,
            }, f)
    if db is not None:
        with sqlite3.Connection(db) as conn:
            cursor = conn.cursor()
            cursor.execute('INSERT OR REPLACE INTO muon_rates '
                    'VALUES (?, ?, ?, ?, ?, ?)',
                    (run, ad, num_veto_windows, usable_livetime,
                        rate, efficiency))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    parser.add_argument('-r', '--run', type=int)
    parser.add_argument('--site', type=int)
    parser.add_argument('--ad', type=int)
    parser.add_argument('-o', '--output')
    parser.add_argument('--update-db',
            help='Optional, register the run info and muon rate with the given '
            'database')
    args = parser.parse_args()
    main2(args.run, args.files, args.site, args.ad, args.output, args.update_db)

