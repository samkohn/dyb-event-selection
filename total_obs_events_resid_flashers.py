"""Loop through all runs and tally the observed prompt spectra.

Backgrounds including accidentals are still included.
This is the raw number of counts,
even ignoring the muon-veto and multiplicity efficiencies,
which can be divided out later using a weighted average.
Therefore the statistical uncertainty for each bin
is the standard Poisson uncertainty.
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
    with sqlite3.Connection(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        #cursor.execute('''SELECT RunNo, Hall FROM runs ORDER BY Hall, RunNo''')
        #cursor.execute('''SELECT runs.RunNo, Hall FROM runs LEFT JOIN
                #(SELECT * FROM num_coincidences_by_run AS t1
                #WHERE t1.Source = "No resid. flasher rate-only 9/17/2020") AS existing
            #ON runs.RunNo = existing.RunNo
            #WHERE existing.RunNo IS NULL;''')
        cursor.execute('''SELECT runs.RunNo, Hall, runs.DetNo
            FROM (runs NATURAL JOIN muon_rates) AS runs
                LEFT JOIN (SELECT * FROM num_coincidences_by_run AS t1
                    WHERE t1.Source = "No resid. flasher rate-only 9/17/2020") AS existing
                ON runs.RunNo = existing.RunNo AND runs.DetNo = existing.DetNo
            WHERE existing.DetNo IS NULL''')
        runs = cursor.fetchall()
        cursor.execute('''SELECT Hall, DetNo, Peak, Resolution
            FROM delayed_energy_fits''')
        delayed_fits = cursor.fetchall()
    energy_bounds = {}
    for row in delayed_fits:
        mean = row['Peak']
        width = row['Resolution']
        upper = mean + 3 * width
        lower = mean - 3 * width
        energy_bounds[row['Hall'], row['DetNo']] = (lower, upper)
    coincidences_by_det = {halldet: 0 for halldet in common.all_ads}
    for run, hall, det in runs:
        if hall_constraint is not None and hall != hall_constraint:
            continue
        #dets = common.dets_for(hall, run)  # Fetch AD numbers given EH and 6/8/7AD period
        #for det in dets:
        for det in [det]:  # hack to avoid rewriting code
            if det_constraint is not None and det != det_constraint:
                continue
            filename = os.path.join(
                    datafile_base,
                    f'EH{hall}',
                    f'hadded_ad{det}',
                    f'out_ad{det}_{run}.root'
            )
            delayed_bounds = energy_bounds[hall, det]
            datafile = ROOT.TFile(filename, 'READ')
            events = datafile.Get('ad_events')
            num_candidates = events.Draw(
                    'energy[0]',
                    'multiplicity == 2 && '
                    f'energy[1] > {delayed_bounds[0]} && '
                    f'energy[1] < {delayed_bounds[1]} && '
                    '!(z[0] > 2200 && x[0]*x[0] + y[0]*y[0] > 500) && '
                    '!(z[1] > 2200 && x[1]*x[1] + y[1]*y[1] > 500) && '
                    f'({delayeds._NH_THU_DIST_TIME_CUT_STR})',
                    'goff'
            )
            if update_db:
                with sqlite3.Connection(database, timeout=60) as conn:
                    cursor = conn.cursor()
                    cursor.execute('''INSERT OR REPLACE INTO num_coincidences_by_run
                        VALUES (?, ?, ?, ?, "No resid. flasher rate-only 9/17/2020")''',
                        (run, det, json.dumps([num_candidates]), json.dumps([1500, 12000])))
            else:
                print(f'Finished Run {run} EH{hall}-AD{det}: {num_candidates}')
            coincidences_by_det[hall, det] += num_candidates
    if update_db:
        with sqlite3.Connection(database, timeout=60) as conn:
            cursor = conn.cursor()
            for (hall, det), num in coincidences_by_det.items():
                if num > 0:
                    cursor.execute('''INSERT OR REPLACE INTO num_coincidences
                        VALUES (?, ?, ?, ?, "No resid. flasher rate-only 9/17/2020")''',
                        (hall, det, json.dumps([num]), json.dumps([1500, 12000])))
    else:
        pprint.pprint(coincidences_by_det)


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
