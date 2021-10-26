"""Sum up the accidentals spectrum for each AD."""

import argparse
from collections import defaultdict
import multiprocessing

import numpy as np

import adevent
import common
import delayeds


def one_file(run_key, data_file_path, energy_lookup, bin_edges):
    import ROOT
    run, site, ad = run_key
    if run % 10 == 0:
        print(run_key)
    data_file_name = data_file_path.format(run=run, site=site, ad=ad)
    data_file = ROOT.TFile(data_file_name, 'READ')
    acc_pairs = data_file.Get('all_pairs')

    bin_edges = np.array(bin_edges)
    num_bins = len(bin_edges) - 1
    low_edge = bin_edges[0]
    up_edge = bin_edges[-1]
    acc_hist = ROOT.TH1F("acc_hist", "acc_hist", num_bins, low_edge, up_edge)
    acc_hist.GetXaxis().Set(num_bins, bin_edges)

    delayed_min, delayed_max = energy_lookup[site, ad]
    acc_pairs.Draw('energy[0] >> acc_hist',
        'multiplicity == 2 && '
        f'energy[0] < {adevent._EMAX_THU} && '
        f'energy[1] > {delayed_min} && energy[1] < {delayed_max} && '
        f'{delayeds._NH_THU_DIST_TIME_CUT_STR}',
        'goff'
    )

    spec_list = []
    for i in range(num_bins):
        spec_list.append(acc_hist.GetBinContent(i + 1))
    return (
        run,
        ad,
        spec_list,
    )


def main(
    bgs_database,
    main_database,
    data_file_path,
    binning_id,
    input_label,
    spectrum_label,
    update_db
):
    import ROOT
    # Fetch all triplets of RunNo, Hall, DetNo to use to find files
    with common.get_db(main_database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT RunNo, Hall, DetNo
            FROM runs NATURAL JOIN hall_dets
            WHERE (RunNo, DetNo) IN (SELECT RunNo, DetNo FROM muon_rates)
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
            WHERE Source = ?''', (input_label,))
        energy_bounds = cursor.fetchall()

    energy_lookup = {}
    for site, ad, low_bound, up_bound in energy_bounds:
        energy_lookup[site, ad] = (low_bound, up_bound)

    def arg_generator():
        for run_key in run_keys:
            yield (run_key, data_file_path, energy_lookup, bin_edges)

    with multiprocessing.Pool() as pool:
        results = pool.starmap(one_file, arg_generator())

    # Set up dicts to hold combined histograms
    # The first time a (site, ad) pair is encountered it will produce 0 (int() == 0)
    # and numpy will add the array to 0 to produce the first stored array.
    spectrum_hists = defaultdict(int)
    # Combine result histograms
    for (run, site, ad), result in zip(run_keys, results):
        spectrum_hists[site, ad] += np.array(result[2], dtype=int)
    # Normalize
    for key, spectrum in spectrum_hists.items():
        spectrum_hists[key] = spectrum / sum(spectrum)

    if update_db:
        spectrum_rows = []
        for (site, ad), spectrum in spectrum_hists.items():
            for i, bin_value in enumerate(spectrum):
                spectrum_rows.append((
                    spectrum_label,
                    site,
                    ad,
                    binning_id,
                    i,
                    bin_value,
                ))
        print(f'Loading in {len(spectrum_rows)} rows')
        print(f'First row: {spectrum_rows[0]}')
        with common.get_db(bgs_database) as conn:
            cursor = conn.cursor()
            cursor.executemany('''
                INSERT OR REPLACE INTO
                    accidentals_spectrum
                VALUES
                    (?, ?, ?, ?, ?, ?)
                ''',
                spectrum_rows
            )
    else:
        print(results[:10])
        print(spectrum_hists)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bgs_database', help='Database for storing output')
    parser.add_argument('main_database', help='Database for runs & for results')
    parser.add_argument('--acc-file-template',
        help='Template for run-by-run acc files, e.g.'
        '/home/skohn/EH{site}/acc_ad{ad}/acc_ad{ad}_{run}.root'
    )
    parser.add_argument('--binning-id', type=int)
    parser.add_argument('--input-label', required=True)
    parser.add_argument('--spec-label', required=True)
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(
        args.bgs_database,
        args.main_database,
        args.acc_file_template,
        args.binning_id,
        args.input_label,
        args.spec_label,
        args.update_db,
    )

