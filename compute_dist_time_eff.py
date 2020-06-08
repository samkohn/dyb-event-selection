import argparse
import os
import sqlite3

def main(input_basepath, database, update_db):
    import ROOT
    with sqlite3.Connection(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('''SELECT Hall, DetNo, Peak, Resolution
        FROM delayed_energy_fits WHERE Hall > 0 AND DetNo > 0''')
        results = cursor.fetchall()
    cut_lookup = {(row['Hall'], row['DetNo']): row for row in results}
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
        cut_integral = delayed_vs_DT.Integral(DT_LOW_BIN, DT_CUT_UP_BIN,
                energy_low_bin, energy_up_bin)
        all_integral = delayed_vs_DT.Integral(DT_LOW_BIN, DT_ALL_UP_BIN,
                energy_low_bin, energy_up_bin)
        only_excluded_integral = delayed_vs_DT.Integral(DT_CUT_UP_BIN+1,
                DT_ALL_UP_BIN, energy_low_bin, energy_up_bin)
        bin_width_error_integral = delayed_vs_DT.Integral(DT_LOW_BIN,
                DT_CUT_UP_BIN, energy_low_bin+1, energy_up_bin-1)
        efficiency = cut_integral/all_integral
        bin_width_changes = bin_width_error_integral/all_integral
        deviation = bin_width_changes - efficiency
        print(f'EH{site}-AD{det}')
        print(f'Nominal: {100*efficiency:.02f}%')
        print(f'Deviation due to 0.005MeV bin width: {100*deviation:.03f}%')
        infile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='directory containing EH1, EH2, EH3')
    parser.add_argument('database')
    parser.add_argument('--update-db')
    args = parser.parse_args()
    main(args.input, args.database, args.update_db)
