"""Loop through all "hadded" data files and save the total number of coincidences."""

import argparse
from collections import defaultdict
import itertools as it
import json
import multiprocessing

import numpy as np

import common
import delayeds
import adevent

def one_file(run_key, data_file_path, energy_lookup, bin_edges):
    import ROOT
    run, site, ad = run_key
    if run % 10 == 0:
        print(run_key)
    data_file_name = data_file_path.format(run=run, site=site, ad=ad)
    data_file = ROOT.TFile(data_file_name, 'READ')
    ad_events = data_file.Get('ad_events')

    bin_edges = np.array(bin_edges)
    num_bins = len(bin_edges) - 1
    low_edge = bin_edges[0]
    up_edge = bin_edges[-1]
    nominal_hist = ROOT.TH1F("nominal_hist", "nominal_hist", num_bins, low_edge, up_edge)
    nominal_hist.GetXaxis().Set(num_bins, bin_edges)
    adtime_hist = ROOT.TH1F("adtime_hist", "adtime_hist", num_bins, low_edge, up_edge)
    adtime_hist.GetXaxis().Set(num_bins, bin_edges)

    delayed_min, delayed_max = energy_lookup['nominal', site, ad]
    num_coincidences_nominal = ad_events.Draw('energy[0] >> nominal_hist',
        'multiplicity == 2 && '
        f'energy[0] < {adevent._EMAX_THU} && '
        f'energy[1] > {delayed_min} && energy[1] < {delayed_max} && '
        f'{delayeds._NH_THU_DIST_TIME_CUT_STR}',
        'goff'
    )
    delayed_min, delayed_max = energy_lookup['adtime', site, ad]
    num_coincidences_adtime = ad_events.Draw('energy[0] >> adtime_hist',
        'multiplicity == 2 && '
        f'energy[0] < {adevent._EMAX_THU} && '
        f'energy[1] > {delayed_min} && energy[1] < {delayed_max} && '
        f'dr_to_prompt_AdTime[1] + {delayeds._NH_THU_DIST_TIME_CONST} '
        f' * dt_to_prompt[1] < {delayeds._NH_THU_DIST_TIME_MAX}',
        'goff'
    )

    nominal_spec_list = []
    adtime_spec_list = []
    for i in range(num_bins):
        nominal_spec_list.append(nominal_hist.GetBinContent(i + 1))
        adtime_spec_list.append(adtime_hist.GetBinContent(i + 1))
    return (
        run,
        ad,
        num_coincidences_nominal,
        num_coincidences_adtime,
        nominal_spec_list,
        adtime_spec_list,
    )

def main(
    main_database,
    data_file_path,
    binning_id,
    save_spectrum,
    save_total,
    labels,
    update_db
):
    import ROOT
    # Fetch all triplets of RunNo, Hall, DetNo to use to find files
    with common.get_db(main_database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT RunNo, Hall, DetNo
            FROM runs NATURAL JOIN muon_rates
            ORDER BY RunNo, Hall, DetNo''')
        run_keys = cursor.fetchall()
        cursor.execute('''
            SELECT
                BinEdgeEnergy_keV
            FROM
                reco_binnings
            WHERE
                Id = ?
            ORDER BY
                BinEdgeIndex
            ''',
            (binning_id,)
        )
        # returned as list of 1-tuples in keV
        bin_edges = [x[0]/1000 for x in cursor.fetchall()]

    # Look up delayed energy cuts
    with common.get_db(main_database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT Hall, DetNo, Peak - 3 * Resolution,
            Peak + 3 * Resolution
            FROM delayed_energy_fits
            WHERE Source = ?''', (labels['nominal'],))
        nominal_energy_bounds = cursor.fetchall()
        cursor.execute('''SELECT Hall, DetNo, Peak - 3 * Resolution,
            Peak + 3 * Resolution
            FROM delayed_energy_fits
            WHERE Source = ?''', (labels['adtime'],))
        adtime_energy_bounds = cursor.fetchall()

    energy_lookup = {}
    for site, ad, low_bound, up_bound in nominal_energy_bounds:
        energy_lookup['nominal', site, ad] = (low_bound, up_bound)
    for site, ad, low_bound, up_bound in adtime_energy_bounds:
        energy_lookup['adtime', site, ad] = (low_bound, up_bound)

    def arg_generator():
        for run_key in run_keys:
            yield (run_key, data_file_path, energy_lookup, bin_edges)

    with multiprocessing.Pool() as pool:
        results = pool.starmap(one_file, arg_generator())

    # Set up dicts to hold combined histograms
    # The first time a (site, ad) pair is encountered it will produce 0 (int() == 0)
    # and numpy will add the array to 0 to produce the first stored array.
    adsimple_hist = defaultdict(int)
    adtime_hist = defaultdict(int)
    # Combine result histograms
    for (run, site, ad), result in zip(run_keys, results):
        adsimple_hist[site, ad] += np.array(result[4], dtype=int)
        adtime_hist[site, ad] += np.array(result[5], dtype=int)

    if update_db:
        spectrum_rows = []
        for (site, ad), spectrum in adsimple_hist.items():
            spectrum_rows.append((
                site,
                ad,
                binning_id,
                json.dumps(spectrum.tolist()),
                "adsimple 1/27/2021",
            ))
        for (site, ad), spectrum in adtime_hist.items():
            spectrum_rows.append((
                site,
                ad,
                binning_id,
                json.dumps(spectrum.tolist()),
                "adtime 1/27/2021",
            ))
        results_nominal = [(*result[:3], labels['nominal']) for result in results]
        results_adtime = [(*result[:2], result[3], labels['adtime']) for result in results]
        with common.get_db(main_database) as conn:
            cursor = conn.cursor()
            if save_total:
                cursor.executemany('''INSERT OR REPLACE INTO num_coincidences_by_run
                    VALUES (?, ?, ?, ?)''', results_nominal + results_adtime)
            if save_spectrum:
                cursor.executemany('''
                    INSERT OR REPLACE INTO
                        num_coincidences
                    VALUES
                        (?, ?, ?, ?, ?)
                    ''',
                    spectrum_rows
                )
    else:
        print(results[:10])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('main_database', help='Database for runs & for results')
    parser.add_argument('--data-file-path',
        help='File path template for run-by-run data files, e.g. '
        '/home/skohn/EH{site}/hadded_ad{ad}/out_ad{ad}_{run}.root'
    )
    parser.add_argument('--binning-id', type=int)
    parser.add_argument('--spectrum', action='store_true')
    parser.add_argument('--total', action='store_true')
    parser.add_argument('--update-db', action='store_true')
    parser.add_argument('--label-nominal', required=True)
    parser.add_argument('--label-adtime', required=True)
    args = parser.parse_args()
    if not (args.total or args.spectrum):
        raise ValueError("Must specify at least one of --total or --spectrum")
    main(args.main_database, args.data_file_path,
        args.binning_id, args.spectrum, args.total,
        {'nominal': args.label_nominal, 'adtime': args.label_adtime},
        args.update_db
    )
