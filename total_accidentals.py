"""Loop through all runs and tally the accidentals counts and uncertainty.

Refer to DocDB-12201 for the derivation of the equations used here.
The rates are *not* corrected for muon-livetime or multiplicity efficiencies.
"""
import argparse
import json
import os.path
import pprint
import sqlite3

import common
import delayeds

def main(database, datafile_base, hall_constraint, det_constraint, update_db):
    import ROOT
    bins = [1500, 12000]
    source = "Nominal rate-only 9/22/2020"
    with sqlite3.Connection(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('''
        SELECT RunNo, DetNo, Hall, DistanceTime_DT_Eff, singles.Rate_Hz,
            singles.IsolatedEventRate_Hz, muons.Livetime_ns
        FROM (runs INNER JOIN singles_rates as singles USING (RunNo))
            INNER JOIN muon_rates as muons USING (RunNo, DetNo)
            INNER JOIN distance_time_eff_study USING (RunNo, DetNo)
        WHERE PairingType = "random_many_expo_time"
        ORDER BY RunNo
        ''')
        runs = cursor.fetchall()
        cursor.execute(''' SELECT Hall, DetNo, Peak, Resolution
            FROM delayed_energy_fits''')
        delayed_fits = cursor.fetchall()
    energy_bounds = {}
    for row in delayed_fits:
        mean = row['Peak']
        width = row['Resolution']
        upper = mean + 3 * width
        lower = mean - 3 * width
        energy_bounds[row['Hall'], row['DetNo']] = (lower, upper)
    accs_by_det = {halldet: 0 for halldet in common.all_ads}
    COINCIDENCE_WINDOW_S = (delayeds._NH_THU_MAX_TIME - delayeds._NH_THU_MIN_TIME)/1e9
    for run, det, hall, acc_DT_eff, r_s, r_1fold, livetime_ns in runs:
        if hall_constraint is not None and hall != hall_constraint:
            continue
        if det_constraint is not None and det != det_constraint:
            continue
        if acc_DT_eff is None:
            continue
        r_ss = r_s * r_1fold * COINCIDENCE_WINDOW_S
        livetime_s = livetime_ns/1e9
        N_acc_all = r_ss * livetime_s * acc_DT_eff
        filename = os.path.join(
                datafile_base,
                f'EH{hall}',
                f'acc_ad{det}',
                f'acc_ad{det}_{run}.root'
        )
        delayed_bounds = energy_bounds[hall, det]
        datafile = ROOT.TFile(filename, 'READ')
        events = datafile.Get('accidentals')
        num_acc_entries = events.GetEntries()
        num_candidates = events.Draw(
                'energy[0]',
                'multiplicity == 2 && '
                f'energy[1] > {delayed_bounds[0]} && '
                f'energy[1] < {delayed_bounds[1]}',
                'goff'
        )
        # scale N_acc_all by fraction of entries with the right delayed energy
        N_acc_candidates = N_acc_all * num_candidates / num_acc_entries
        accs_by_det[hall, det] += N_acc_candidates
        if update_db:
            with sqlite3.Connection(database) as conn:
                cursor = conn.cursor()
                cursor.execute('''INSERT OR REPLACE INTO num_accs_by_run
                    VALUES (?, ?, ?, ?, ?)''',
                    (run, det, json.dumps([N_acc_candidates]), json.dumps(bins), source))
        print(f'Finished Run {run} EH{hall}-AD{det}')
    if update_db:
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            for (hall, det), num in accs_by_det.items():
                if num > 0:
                    cursor.execute('''INSERT OR REPLACE INTO num_accs
                        VALUES (?, ?, ?, ?, ?)''',
                        (hall, det, json.dumps([num]), json.dumps(bins), source))
    else:
        pprint.pprint(accs_by_det)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('datafile_base',
            help='base directory of data files with subdirs EH1, EH2, EH3')
    parser.add_argument('--hall', type=int)
    parser.add_argument('--det', type=int)
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()

    main(args.database, args.datafile_base, args.hall, args.det, args.update_db)
