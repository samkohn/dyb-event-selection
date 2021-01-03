"""Loop through all "acc" data files and save the total accidental sample efficiency.

This is the fraction of all "ss" pairs that pass the IBD cut (energy && DT).
"""

import argparse
import itertools as it
import math
import multiprocessing

import common
import delayeds

def one_file(run_key, data_file_path, energy_lookup):
    import ROOT
    run, site, ad = run_key
    data_file_name = data_file_path.format(run=run, site=site, ad=ad)
    data_file = ROOT.TFile(data_file_name, 'READ')
    all_pairs = data_file.Get('all_pairs')
    delayed_min, delayed_max = energy_lookup[site, ad]
    num_passing_cuts = all_pairs.Draw('energy[0]',
        f'energy[1] > {delayed_min} && energy[1] < {delayed_max} && '
        f'{delayeds._NH_THU_DIST_TIME_CUT_STR}',
        'goff'
    )
    total = all_pairs.GetEntries()
    efficiency = num_passing_cuts / total
    # Binomial error / total = absolute error on efficiency
    error = math.sqrt(total * efficiency * (1 - efficiency)) / total
    return efficiency, error, run, ad

def main(main_database, energy_cuts_database, data_file_path, label, update_db):
    import ROOT
    # Fetch all triplets of RunNo, Hall, DetNo to use to find files
    with common.get_db(main_database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT RunNo, Hall, DetNo
            FROM runs NATURAL JOIN accidental_subtraction
            WHERE Label = ?
            ORDER BY RunNo, Hall, DetNo''', (label, ))
        run_keys = cursor.fetchall()

    # Look up delayed energy cuts
    with common.get_db(energy_cuts_database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT Hall, DetNo, Peak - 3 * Resolution,
            Peak + 3 * Resolution
            FROM delayed_energy_fits''')
        energy_bounds = cursor.fetchall()
    energy_lookup = {}
    for site, ad, low_bound, up_bound in energy_bounds:
        energy_lookup[site, ad] = (low_bound, up_bound)

    with multiprocessing.Pool(63) as pool:
        results = pool.starmap(one_file, zip(run_keys, it.repeat(data_file_path),
            it.repeat(energy_lookup)))
    if update_db:
        new_results = []
        for x in results:
            new_results.append(x + (label,))
        results = new_results

        with common.get_db(main_database) as conn:
            cursor = conn.cursor()
            cursor.executemany('''UPDATE accidental_subtraction
                SET Total_Acc_Eff = ?,
                    Total_Acc_Eff_err = ?
                WHERE RunNo = ? AND DetNo = ? AND Label = ?''',
            results)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('main_database', help='Database for runs & for results')
    parser.add_argument('energy_cuts_database',
        help='Database for delayed energy cut lookup'
    )
    parser.add_argument('--acc-file-path',
        help='File path template for run-by-run accidentals files, e.g. '
        '/home/skohn/EH{site}/acc_ad{ad}/acc_ad{ad}_{run}.root'
    )
    parser.add_argument('--label')
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(args.main_database, args.energy_cuts_database, args.acc_file_path,
        args.label, args.update_db)
