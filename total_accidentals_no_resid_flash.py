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
    source = "No resid. flashers rate-only 9/17/2020"
    with sqlite3.Connection(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        #cursor.execute('''
        #SELECT RunNo, DetNo, Hall, DistanceTime_DT_Eff, singles.Rate_Hz,
            #singles.IsolatedEventRate_Hz, muons.Livetime_ns
        #FROM (runs INNER JOIN singles_rates as singles USING (RunNo))
            #INNER JOIN muon_rates as muons USING (RunNo, DetNo)
            #INNER JOIN distance_time_eff_study USING (RunNo, DetNo)
        #WHERE PairingType = "random_many_resid_flasher"
        #ORDER BY RunNo
        #''')
        cursor.execute('''
        SELECT RunNo, DetNo, Hall, DistanceTime_DT_Eff, muons.Livetime_ns
        FROM runs INNER JOIN muon_rates as muons USING (RunNo)
            INNER JOIN distance_time_eff_study USING (RunNo, DetNo)
        WHERE PairingType = "random_many_resid_flasher"
        ORDER BY RunNo, DetNo
        ''')
        runs = cursor.fetchall()
        cursor.execute(''' SELECT Hall, DetNo, Peak, Resolution
            FROM delayed_energy_fits''')
        delayed_fits = cursor.fetchall()
    with sqlite3.Connection('/global/cscratch1/sd/skohn/dyb30/parameters_resid_flashers.db') as conn:
        cursor = conn.cursor()
        cursor.execute('''
        SELECT RunNo, DetNo, Hall, singles.Rate_Hz, singles.IsolatedEventRate_Hz
        FROM runs INNER JOIN singles_rates as singles USING (RunNo)
        ORDER BY RunNo, DetNo
        ''')
        singles = cursor.fetchall()
    energy_bounds = {}
    for row in delayed_fits:
        mean = row['Peak']
        width = row['Resolution']
        upper = mean + 3 * width
        lower = mean - 3 * width
        energy_bounds[row['Hall'], row['DetNo']] = (lower, upper)
    accs_by_det = {halldet: 0 for halldet in common.all_ads}
    COINCIDENCE_WINDOW_S = (delayeds._NH_THU_MAX_TIME - delayeds._NH_THU_MIN_TIME)/1e9
    #for run, det, hall, acc_DT_eff, r_s, r_1fold, livetime_ns in runs:
    index1 = 0
    index2 = 0
    n_skipped = 0
    while index1 < len(runs) and index2 < len(singles):
        run, det, hall, acc_DT_eff, livetime_ns = runs[index1]
        run2, det2, hall2, r_s, r_1fold = singles[index2]
        #for (run, det, hall, acc_DT_eff, livetime_ns), (run2,det2, hall2, r_s, r_1fold) in zip(runs,
                #singles):
        if run != run2 or det!= det2 or hall != hall2:
            index2 += 1
            if update_db:
                with sqlite3.Connection(database) as conn:
                    cursor = conn.cursor()
                    cursor.execute('''DELETE FROM num_accs_by_run
                        WHERE RunNo = ? AND DetNo = ? AND BinEdges = ? AND Source = ?''',
                    (run, det, json.dumps(bins), source))
            continue
            #raise ValueError()
        if hall_constraint is not None and hall != hall_constraint:
            index1 += 1
            index2 += 1
            continue
        if det_constraint is not None and det != det_constraint:
            index1 += 1
            index2 += 1
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
        num_acc_entries = events.Draw(
                'energy[0]',
                'multiplicity == 2 && '
                '!(z[0] > 2200 && x[0]*x[0] + y[0]*y[0] > 500) && '
                '!(z[1] > 2200 && x[1]*x[1] + y[1]*y[1] > 500)',
                'goff'
        )
        num_candidates = events.Draw(
                'energy[0]',
                'multiplicity == 2 && '
                f'energy[1] > {delayed_bounds[0]} && '
                f'energy[1] < {delayed_bounds[1]} &&'
                '!(z[0] > 2200 && x[0]*x[0] + y[0]*y[0] > 500) && '
                '!(z[1] > 2200 && x[1]*x[1] + y[1]*y[1] > 500)',
                'goff'
        )
        # scale N_acc_all by fraction of entries with the right energy
        N_acc_candidates = N_acc_all * num_candidates / num_acc_entries
        accs_by_det[hall, det] += N_acc_candidates
        if update_db:
            with sqlite3.Connection(database) as conn:
                cursor = conn.cursor()
                cursor.execute('''INSERT OR REPLACE INTO num_accs_by_run
                    VALUES (?, ?, ?, ?, ?)''',
                    (run, det, json.dumps([N_acc_candidates]), json.dumps(bins), source))
        print(f'Finished Run {run} EH{hall}-AD{det}')
        index1 += 1
        index2 += 1
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
