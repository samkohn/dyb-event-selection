import argparse
import math
import os
import sqlite3

import common

def get_stat_error(n_passes_cut, n_fails_cut, error_passes, error_fails):
    """Compute the statistical error on the efficiency.

    This function performs error propagation on the expression:

        Efficiency = N(passes DT cut) / N(total)

    including the partial correlation between the numerator and the
    denominator. This is accomplished by assigning:

        N(total) = N(passes) + N(fails)

    Then rewriting the expression as:

        Efficiency = 1 / (1 + N(fails) / N(passes))

    This expression now has no hidden statistical correlations due to the same
    histogram bin being counted twice, which the original expression does have.
    The final mathematical expression is:

        Error = ...
            N(fails)/N(passes) * Error(N(fails)/N(passes)) / (1 + fails/passes)^2
    """
    value_ratio = n_fails_cut / n_passes_cut
    passes_relative_error_sq = pow(error_passes / n_passes_cut, 2)
    fails_relative_error_sq = pow(error_fails / n_fails_cut, 2)
    quotient_propagation = math.sqrt(passes_relative_error_sq +
            fails_relative_error_sq)
    power_term = pow(1 + value_ratio, 2)
    return value_ratio * quotient_propagation/power_term

def main(input_basepath, database, update_db):
    import ROOT
    with common.get_db(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('''SELECT Hall, DetNo, Peak, Resolution
        FROM delayed_energy_fits WHERE Hall > 0 AND DetNo > 0''')
        results = cursor.fetchall()
    cut_lookup = {(row['Hall'], row['DetNo']): row for row in results}
    effs = {}
    stat_errors = {}
    binning_errors = {}
    DT_LOW_BIN = 1
    DT_CUT_UP_BIN = 16  # 800mm
    DT_ALL_UP_BIN = 60  # 3000mm
    for (site, det), cuts in cut_lookup.items():
        site_dir = 'EH{}'.format(site)
        sub_dir = 'sub_ad{}'.format(det)
        name = 'sub_ad{}.root'.format(det)
        path = os.path.join(input_basepath, site_dir, sub_dir, name)
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
        only_excluded_integral = delayed_vs_DT.IntegralAndError(DT_CUT_UP_BIN+1,
                DT_ALL_UP_BIN, energy_low_bin, energy_up_bin,
                only_excluded_error)
        stat_error = get_stat_error(cut_integral, only_excluded_integral,
                cut_error.real, only_excluded_error.real)
        bin_width_error_integral = delayed_vs_DT.Integral(DT_LOW_BIN,
                DT_CUT_UP_BIN, energy_low_bin+1, energy_up_bin-1)
        efficiency = cut_integral/all_integral
        bin_width_changes = bin_width_error_integral/all_integral
        bin_width_error = (efficiency - bin_width_changes)/2

        print(f'EH{site}-AD{det}')
        print(f'Nominal: {100*efficiency:.02f}%')
        print(f'Deviation due to 0.005MeV bin width: {100*bin_width_error:.03f}%')
        print(f'Statistical error: {100*stat_error:.03f}%')
        print(f'N(passes): {cut_integral} +/- {cut_error.real}')
        print(f'N(fails):  {only_excluded_integral} +/- {only_excluded_error.real}')
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
                (?, ?, ?, ?, ?)''', (site, det, eff, stat_errors[site, det],
                    binning_errors[site, det]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='directory containing EH1, EH2, EH3')
    parser.add_argument('database')
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(args.input, args.database, args.update_db)
