import argparse
import math
import sqlite3

import delayeds

def main(infile, update_db, pairing, run, detector, selection_str):
    import ROOT
    accfile = ROOT.TFile(infile, 'READ')
    all_pairs = accfile.Get('all_pairs')
    num_pass_DT_cut = all_pairs.Draw('energy[0]',
            f'({delayeds._NH_THU_DIST_TIME_CUT_STR}) && ({selection_str})',
            'goff')
    num_pairs = all_pairs.Draw('energy[0]', selection_str, 'goff')
    efficiency = num_pass_DT_cut / num_pairs
    error = math.sqrt(num_pass_DT_cut * (1 - efficiency)) / num_pairs  # binomial
    try:
        percent_error = 100 * error / efficiency
    except ZeroDivisionError:
        percent_error = 0
    if update_db is None:
        print(f'Pairing type: {pairing}')
        print(f'Efficiency: {efficiency:.6f} +/- {error:.6f} ({percent_error:.1f}%)')
        print(f'Total pairs: {num_pairs}')
        print(f'Passed DT cut: {num_pass_DT_cut}')
    else:
        all_pairs.GetEntry(0)
        with sqlite3.Connection(update_db) as conn:
            cursor = conn.cursor()
            cursor.execute('''INSERT OR REPLACE INTO distance_time_eff_study
                VALUES (?, ?, ?, ?, ?, ?)''',
                (run, detector, pairing, efficiency, error, num_pairs))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--update-db')
    parser.add_argument('--pairing')
    parser.add_argument('-r', '--run', type=int)
    parser.add_argument('-d', '--detector', type=int, choices=(1, 2, 3, 4))
    parser.add_argument('--selection', default='1')
    args = parser.parse_args()
    main(args.infile, args.update_db, args.pairing, args.run, args.detector,
            args.selection)
