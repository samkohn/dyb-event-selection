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
        livetimes = stats['usable_livetime']
    import ROOT
    ad_map = {
            1: [1, 2],
            2: [1, 2],
            4: [1, 2, 3, 4],
            }
    ch = ROOT.TChain('computed')
    ch.Add(filename)
    ch.GetEntry(0)
    runNo = ch.run
    site = ch.site
    ch.Draw('>>basic_eventlist', 'energy < 12 && energy > 1.5 && '
        '!tag_AnyMuonVeto && '
        '(tag_flasher == 0 || tag_flasher == 2) && '
        '(dt_coinc_window_start + dt_next_WSMuon > 400000)',
        'goff')
    basic_eventlist = ROOT.gDirectory.Get('basic_eventlist')
    ch.SetEventList(basic_eventlist)
    for AD in ad_map[site]:
        upper = ch.Draw('energy', 'detector == {}'.format(AD), 'goff')
        lower = ch.Draw('energy', ('detector == {} && coincidence_number =='
                '1').format(AD), 'goff')
        lower += 2 * ch.Draw('energy', ('detector == {} && coincidence_number '
                '== 2 && dr_previous_PromptLike > 500').format(AD), 'goff')
        counts = (upper + lower) / 2.
        error = (upper - lower) / 2.
        rate = counts/(livetimes[str(AD)]/1e9)
        rate_error = error/(livetimes[str(AD)]/1e9)
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
            cursor.execute('INSERT OR REPLACE INTO singles_rates '
                    'VALUES (?, ?, ?, ?, ?, ?, ?)',
                    (runNo, AD, int(counts), livetimes[str(AD)], rate, error,
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
