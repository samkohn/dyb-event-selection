import argparse
import sqlite3
import os
from array import array
import math

import common

MAGIC_LOW_BOUND = 1.6
MAGIC_UP_BOUND = 2.8

def main(input_basepath, database, update_db):
    import ROOT
    with common.get_db(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('''SELECT Hall, DetNo, Peak, Resolution FROM
            delayed_energy_fits
            WHERE
                Hall > 0
                AND DetNo > 0
                AND Source = "nominal"''')  #TODO hardcoded
        results = cursor.fetchall()
    cut_lookup = {(row['Hall'], row['DetNo']): row for row in results}
    values_lookup = {}
    errors_lookup = {}
    for (site, det), cuts in cut_lookup.items():
        site_dir = 'EH{}'.format(site)
        sub_dir = 'sub_nominal_ad{}'.format(det)  #TODO hardcoded
        name = 'sub_ad{}.root'.format(det)
        path = os.path.join(input_basepath, site_dir, sub_dir, name)
        print(path)
        infile = ROOT.TFile(path, 'READ')
        spectrum_2d = infile.Get('final')
        delayed_spectrum = spectrum_2d.ProjectionY()
        low_limit = cuts['Peak'] - 3 * cuts['Resolution']
        up_limit = cuts['Peak'] + 3 * cuts['Resolution']
        low_bin = delayed_spectrum.FindBin(low_limit)
        high_bin = delayed_spectrum.FindBin(up_limit)
        nominal_value = delayed_spectrum.Integral(low_bin, high_bin)
        nominal_error = math.sqrt(sum(delayed_spectrum.GetBinError(index)**2 for index in
                range(low_bin, high_bin + 1)))
        print(nominal_value)
        low_bin = delayed_spectrum.FindBin(MAGIC_LOW_BOUND)
        high_bin = delayed_spectrum.FindBin(MAGIC_UP_BOUND)
        wide_value = delayed_spectrum.Integral(low_bin, high_bin)
        wide_error = math.sqrt(sum(delayed_spectrum.GetBinError(index)**2 for index in
                range(low_bin, high_bin + 1)))
        print(wide_value)
        values_lookup[site, det] = (nominal_value, wide_value)
        errors_lookup[site, det] = (nominal_error, wide_error)
        infile.Close()
    xvals = array('f', [x[1] for x in values_lookup.values()])
    yvals = array('f', [x[0] for x in values_lookup.values()])
    xerrs = array('f', [x[1] for x in errors_lookup.values()])
    yerrs = array('f', [x[0] for x in errors_lookup.values()])
    print(errors_lookup)
    graph_to_fit = ROOT.TGraphErrors(len(values_lookup), xvals, yvals, xerrs,
            yerrs)
    graph_to_fit.Draw('A*')
    fit_result = graph_to_fit.Fit('pol1', 'QS')
    y_intercept, slope = fit_result.Parameter(0), fit_result.Parameter(1)
    if update_db:
        with common.get_db(database) as conn:
            cursor = conn.cursor()
            for (site, det), (nominal, wide) in values_lookup.items():
                model_value = y_intercept + slope * wide
                relative_deviation = 1 - model_value / nominal
                print(site, det, model_value, relative_deviation)
                error = errors_lookup[site, det][0] / nominal
                cursor.execute('''INSERT OR REPLACE INTO delayed_energy_uncertainty_1
                    VALUES (?, ?, ?, ?)''', (site, det, relative_deviation, error))
        ROOT.gPad.Print('eff_uncertainty.pdf')
    else:
        intercept_error = fit_result.ParError(0)
        slope_error = fit_result.ParError(1)
        print(f'y intercept = {y_intercept} +/- {intercept_error}')
        print(f'slope = {slope} +/- {slope_error}')
        for (site, det), (nominal, wide) in values_lookup.items():
            model_value = y_intercept + slope * wide
            relative_deviation = 1 - model_value / nominal
            print(site, det, model_value, relative_deviation)


if __name__ == '__main__':
    raise RuntimeError("You better double check the Source and file path"
        " and comment out this line before running for real.")
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('database')
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(args.input, args.database, args.update_db)
