"""Loop through all "hadded" data files and save the total number of coincidences."""

import argparse
import itertools as it
import multiprocessing
import sqlite3

import delayeds

def one_file(run_key, data_file_path, energy_lookup):
    import ROOT
    run, site, ad = run_key
    data_file_name = data_file_path.format(run=run, site=site, ad=ad)
    data_file = ROOT.TFile(data_file_name, 'READ')
    ad_events = data_file.Get('ad_events')
    delayed_min, delayed_max = energy_lookup['nominal', site, ad]
    num_coincidences_nominal = ad_events.Draw('energy[0]',
        f'energy[1] > {delayed_min} && energy[1] < {delayed_max} && '
        f'{delayeds._NH_THU_DIST_TIME_CUT_STR}',
        'goff'
    )
    delayed_min, delayed_max = energy_lookup['adtime', site, ad]
    num_coincidences_adtime = ad_events.Draw('energy[0]',
        f'energy[1] > {delayed_min} && energy[1] < {delayed_max} && '
        f'dr_to_prompt_AdTime[1] + {delayeds._NH_THU_DIST_TIME_CONST} '
        f' * dt_to_prompt[1] < {delayeds._NH_THU_DIST_TIME_MAX}',
        'goff'
    )
    return run, ad, num_coincidences_nominal, num_coincidences_adtime

def main(main_database, energy_cuts_database, data_file_path, update_db):
    import ROOT
    # Fetch all triplets of RunNo, Hall, DetNo to use to find files
    with sqlite3.Connection(main_database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT RunNo, Hall, DetNo
            FROM runs NATURAL JOIN muon_rates
            ORDER BY RunNo, Hall, DetNo''')
        run_keys = cursor.fetchall()

    # Look up delayed energy cuts
    with sqlite3.Connection(energy_cuts_database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT Hall, DetNo, Peak - 3 * Resolution,
            Peak + 3 * Resolution
            FROM delayed_energy_fits
            WHERE Source = "nominal"''')
        nominal_energy_bounds = cursor.fetchall()
        cursor.execute('''SELECT Hall, DetNo, Peak - 3 * Resolution,
            Peak + 3 * Resolution
            FROM delayed_energy_fits
            WHERE Source = "adtime"''')
        adtime_energy_bounds = cursor.fetchall()

    energy_lookup = {}
    for site, ad, low_bound, up_bound in nominal_energy_bounds:
        energy_lookup['nominal', site, ad] = (low_bound, up_bound)
    for site, ad, low_bound, up_bound in adtime_energy_bounds:
        energy_lookup['adtime', site, ad] = (low_bound, up_bound)

    with multiprocessing.Pool(63) as pool:
        results = pool.starmap(one_file, zip(run_keys, it.repeat(data_file_path),
            it.repeat(energy_lookup)))
    if update_db:
        results_nominal = [(run, ad, num_nominal, 'new nominal') for run, ad, num_nominal, _ in results]
        results_adtime = [(run, ad, num_adtime, 'new adtime') for run, ad, _, num_adtime in results]
        with sqlite3.Connection(main_database) as conn:
            cursor = conn.cursor()
            cursor.executemany('''INSERT OR REPLACE INTO num_coincidences_by_run
                VALUES (?, ?, ?, ?)''', results_nominal + results_adtime)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('main_database', help='Database for runs & for results')
    parser.add_argument('energy_cuts_database',
        help='Database for delayed energy cut lookup'
    )
    parser.add_argument('--data-file-path',
        help='File path template for run-by-run data files, e.g. '
        '/home/skohn/EH{site}/hadded_ad{ad}/out_ad{ad}_{run}.root'
    )
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(args.main_database, args.energy_cuts_database, args.data_file_path,
        args.update_db)
