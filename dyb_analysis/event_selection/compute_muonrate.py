from __future__ import print_function
import argparse
import sqlite3
import os.path
import json

def main(filename, update_db):
    if update_db:
        conn = sqlite3.Connection('/global/u2/s/skohn/dyb-event-selection-production'
            '/parameters.db')
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
    ch.Draw('>>basic_eventlist',
        '(tag_WSMuon || tag_ADMuon || tag_ShowerMuon) && '
        '(tag_flasher == 0 || tag_flasher == 2) && '
        'dt_previous_WSMuon > 400000 && dt_previous_ADMuon > 800000 &&'
        'dt_previous_ShowerMuon > 1000000000',
        'goff')
    basic_eventlist = ROOT.gDirectory.Get('basic_eventlist')
    ch.SetEventList(basic_eventlist)
    for AD in ad_map[site]:
        counts = ch.Draw('energy', ('detector == 5 || detector == 6 ||'
            'detector == {}').format(AD), 'goff')
        rate = counts/(livetimes[str(AD)]/1e9)
        print(AD)
        print('Counts:', counts)
        print('Rate:', rate, 'Hz')
        if update_db:
            cursor.execute('INSERT OR REPLACE INTO muon_rates '
                    'VALUES (?, ?, ?, ?, ?)',
                    (runNo, AD, int(counts), livetimes[str(AD)], rate))
            conn.commit()
    conn.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(args.infile, args.update_db)
