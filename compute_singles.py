from __future__ import print_function
import argparse
import sqlite3
import os.path
import json

def main(filename, update_db):
    if update_db:
        conn = sqlite3.Connection('parameters.db')
        cursor = conn.cursor()
    with open(os.path.splitext(filename)[0] + '.json', 'r') as f:
        stats = json.load(f)
        livetime = stats['usable_livetime']
    import ROOT
    ad = 1
    #ad_map = {
            #1: [1, 2],
            #2: [1, 2],
            #4: [1, 2, 3, 4],
            #}
    ch = ROOT.TChain('ad_events')
    ch.Add(filename)
    ch.GetEntry(0)
    runNo = ch.run
    site = ch.site
    start_time = ch.timestamp[0]
    for AD in [ad]:
        upper = ch.Draw('energy', 'detector == {}'.format(AD), 'goff')
        lower = ch.Draw('energy', ('detector == {} && multiplicity =='
                '1').format(AD), 'goff')
        lower += ch.Draw('energy', ('detector == {} && multiplicity '
                '== 2 && dr_to_prompt[1] > 500').format(AD), 'goff')
        counts = (upper + lower) / 2.
        error = (upper - lower) / 2.
        rate = counts/(livetime/1e9)
        rate_error = error/(livetime/1e9)
        print(AD)
        print('Upper:', upper)
        print('Lower:', lower)
        print('Counts:', counts)
        print('Error:', error)
        print('Rate:', rate, 'Hz')
        print('Rate error:', rate_error, 'Hz')
        if counts > 0:
            print('Pct:', error/counts)
        if update_db:
            cursor.execute('SELECT COUNT(*) FROM runs WHERE RunNo = ?',
                    (runNo,))
            if cursor.fetchone()[0] == 0:
                cursor.execute('INSERT INTO runs VALUES (?, ?, ?)',
                        (runNo, site, start_time))
            cursor.execute('INSERT OR REPLACE INTO singles_rates '
                    'VALUES (?, ?, ?, ?, ?, ?, ?)',
                    (runNo, AD, int(counts), livetime, rate, error,
                        rate_error))
            conn.commit()
    if update_db:
        conn.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(args.infile, args.update_db)
