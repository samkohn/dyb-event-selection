"""Run the "normalized" version of accidentals subtraction."""

import argparse
import math
import sqlite3

import common

NORMALIZATION_BIN_LOW = 61  # 3m
NORMALIZATION_BIN_UP = 200  # 10m

def main(database, sub_file_template, label, new_label, update_db):
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
    rows = []
    for (site, ad), cuts in cut_lookup.items():
        infilename = sub_file_template.format(site=site, ad=ad)
        infile = ROOT.TFile(infilename, 'READ')
        ed_DT_data = infile.Get('ed_DT_data')
        ed_DT_bg = infile.Get('ed_DT_bg')
        low_limit = cuts['Peak'] - 3 * cuts['Resolution']
        up_limit = cuts['Peak'] + 3 * cuts['Resolution']
        energy_low_bin = ed_DT_data.GetYaxis().FindBin(low_limit)
        energy_up_bin = ed_DT_data.GetYaxis().FindBin(up_limit)
        data_integral = ed_DT_data.Integral(
            NORMALIZATION_BIN_LOW,
            NORMALIZATION_BIN_UP,
            energy_low_bin,
            energy_up_bin,
        )
        bg_integral = ed_DT_bg.Integral(
            NORMALIZATION_BIN_LOW,
            NORMALIZATION_BIN_UP,
            energy_low_bin,
            energy_up_bin,
        )
        correction_factor = data_integral / bg_integral
        print(f'EH{site}-AD{ad} correction factor: {correction_factor-1:.5f}')
        ed_DT_sub = ed_DT_bg.Clone("ed_DT_sub")
        ed_DT_bg.Scale(correction_factor)
        ed_DT_sub.Add(ed_DT_data, ed_DT_bg, 1, -1)
        new_integral = ed_DT_sub.Integral(
            NORMALIZATION_BIN_LOW,
            NORMALIZATION_BIN_UP,
            energy_low_bin,
            energy_up_bin,
        )
        print(f'EH{site}-AD{ad} new integral: {new_integral:.3e}')
        eff_cut_integral = ed_DT_sub.Integral(1, 16, energy_low_bin, energy_up_bin)
        all_integral = ed_DT_sub.Integral(1, 200, energy_low_bin, energy_up_bin)
        efficiency = eff_cut_integral/all_integral
        print(f'EH{site}-AD{ad} efficiency: {efficiency*100:.2f}%')
        p = eff_cut_integral/all_integral
        q = 1 - p
        n = all_integral
        rows.append((site, ad, new_label, efficiency, math.sqrt(p*q/n), None))
    if update_db:
        with common.get_db(database) as conn:
            cursor = conn.cursor()
            cursor.executemany('''
                INSERT INTO
                    distance_time_cut_efficiency
                VALUES
                    (?, ?, ?, ?, ?, ?)
                ''',
                rows,
            )



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('sub_template')
    parser.add_argument('--label')
    parser.add_argument('--update-db', action='store_true')
    parser.add_argument('--db-new-label')
    args = parser.parse_args()
    main(args.database, args.sub_template, args.label, args.db_new_label, args.update_db)
