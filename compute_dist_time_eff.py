import argparse
import math
import os
import sqlite3

import common

def get_stat_error(n_passes_cut, n_total):
    """Binomial error = sqrt(p_passes * p_fails/n_total)."""
    p_passes = n_passes_cut/n_total
    p_fails = 1 - p_passes
    return math.sqrt(p_passes * p_fails / n_total)

def main(file_template, database, label, run_crosscheck, update_db):
    import ROOT
    with common.get_db(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                Hall,
                DetNo,
                Peak,
                Resolution
            FROM
                delayed_energy_fits
            WHERE
                Hall > 0
                AND DetNo > 0
                AND Source = ?
            ''',
            (label,)
        )
        results = cursor.fetchall()
    cut_lookup = {(row['Hall'], row['DetNo']): row for row in results}
    effs = {}
    stat_errors = {}
    binning_errors = {}
    DT_LOW_BIN = 1
    DT_CUT_UP_BIN = 16  # 800mm
    DT_ALL_UP_BIN = 60  # 3000mm
    for (site, det), cuts in cut_lookup.items():
        path = file_template.format(site=site, ad=det)
        print(path)
        infile = ROOT.TFile(path, 'READ')
        delayed_vs_DT = infile.Get('ed_DT_sub')
        low_limit = cuts['Peak'] - 3 * cuts['Resolution']
        up_limit = cuts['Peak'] + 3 * cuts['Resolution']
        print(low_limit, up_limit)
        energy_low_bin = delayed_vs_DT.GetYaxis().FindBin(low_limit)
        energy_up_bin = delayed_vs_DT.GetYaxis().FindBin(up_limit)
        cut_error = ROOT.Double()
        only_excluded_error = ROOT.Double()
        cut_integral = delayed_vs_DT.IntegralAndError(DT_LOW_BIN, DT_CUT_UP_BIN,
                energy_low_bin, energy_up_bin, cut_error)
        all_integral = delayed_vs_DT.Integral(DT_LOW_BIN, DT_ALL_UP_BIN,
                energy_low_bin, energy_up_bin)
        stat_error = get_stat_error(cut_integral, all_integral)
        bin_width_error_integral = delayed_vs_DT.Integral(DT_LOW_BIN,
                DT_CUT_UP_BIN, energy_low_bin+1, energy_up_bin-1)
        efficiency = cut_integral/all_integral
        bin_width_changes = bin_width_error_integral/all_integral
        bin_width_error = (efficiency - bin_width_changes)/2

        if run_crosscheck:
            DT_3m_bin = 61
            DT_5m_bin = 100
            crosscheck_error_with_Ed_cut = ROOT.Double()
            crosscheck_error_without_Ed_cut = ROOT.Double()
            crosscheck_with_Ed_cut = delayed_vs_DT.IntegralAndError(
                DT_3m_bin, DT_5m_bin, energy_low_bin, energy_up_bin,
                crosscheck_error_with_Ed_cut
            )
            crosscheck_without_Ed_cut = delayed_vs_DT.IntegralAndError(
                DT_3m_bin, DT_5m_bin, 1, delayed_vs_DT.GetNbinsY(),
                crosscheck_error_without_Ed_cut
            )
            print(f'EH{site}-AD{det}')
            print(f'Delayed energy range: {low_limit:.4f} -> {up_limit:.4f} MeV')
            print(f'Delayed range bins: {energy_low_bin} -> {energy_up_bin}')
            print(f'Integral: 3m->5m, with Ed cut: {crosscheck_with_Ed_cut:.1f} +/- '
                f'{crosscheck_error_with_Ed_cut:.1f}'
            )
            print(f'Integral: 3m->5m, without Ed cut: {crosscheck_without_Ed_cut:.1f} +/- '
                f'{crosscheck_error_without_Ed_cut:.1f}'
            )
            DT_sub = infile.Get('DT_sub')
            plain_DT_error = ROOT.Double()
            plain_DT_integral = DT_sub.IntegralAndError(
                DT_3m_bin, DT_5m_bin, plain_DT_error
            )
            print(f'DT distribution integral, 3m->5m: {plain_DT_integral:.1f} +/- '
                f'{plain_DT_error:.1f}'
            )


        if not update_db:
            print(f'EH{site}-AD{det}')
            print(f'Nominal: {100*efficiency:.02f}%')
            print(f'Deviation due to 0.005MeV bin width: {100*bin_width_error:.03f}%')
            print(f'Statistical error: {100*stat_error:.03f}%')
            print(f'N(passes): {cut_integral} +/- {stat_error * cut_integral}')
            print(f'N(total):  {all_integral}')
        infile.Close()

        effs[site, det] = efficiency
        stat_errors[site, det] = stat_error
        binning_errors[site, det] = bin_width_error
    if update_db:
        with common.get_db(database) as conn:
            cursor = conn.cursor()
            for (site, det), eff in effs.items():
                cursor.execute('''INSERT OR REPLACE INTO
                distance_time_cut_efficiency VALUES
                (?, ?, ?, ?, ?, ?)''', (site, det, label, eff, stat_errors[site, det],
                    binning_errors[site, det]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_template', help='template for .../EH{site}/sub_ad{ad}/etc.')
    parser.add_argument('database')
    parser.add_argument('--update-db', action='store_true')
    parser.add_argument('--label',
        help='database label for delayed fits and to store efficiency values'
    )
    parser.add_argument('--crosscheck', action='store_true',
        help='print cross-check value of integral from DT=3m to 5m'
    )
    args = parser.parse_args()
    main(args.file_template, args.database, args.label, args.crosscheck, args.update_db)
