"""Compute and dump the bg normalizations from the MC spec file."""
import argparse

import common

def main(infilename, database, bg_index, label):

    livetime_index = '2'
    muon_eff_index = '3'
    mult_eff_index = '4'
    with open(infilename, 'r') as infile:
        for line in infile:
            if line[0] == '#' or len(line) < 10:
                continue
            splits = line.split()
            if splits[1] == livetime_index:
                livetimes = [float(val) for val in splits[2:]]
            if splits[1] == muon_eff_index:
                muon_effs = [float(val) for val in splits[2:]]
            if splits[1] == mult_eff_index:
                mult_effs = [float(val) for val in splits[2:]]
            if splits[1] == str(bg_index):
                bg_rates = [float(val) for val in splits[2:]]
                break
        else:  # nobreak
            raise ValueError(f"Couldn't find background index {bg_index}")
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        for (hall, det), livetime, mu_eff, mult_eff, bg_rate in zip(
            common.all_ads, livetimes, muon_effs, mult_effs, bg_rates
        ):
            cursor.execute('''
                INSERT OR REPLACE INTO
                    accidentals_counts
                VALUES
                    (?, ?, ?, ?)
                ''',
                (label, hall, det, livetime * bg_rate * mu_eff * mult_eff),
            )
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Theta13-etc. ToyMC file')
    parser.add_argument('database')
    parser.add_argument('--bg-index', type=int, help='"row" number in input file')
    parser.add_argument('--label')
    args = parser.parse_args()
    main(args.infile, args.database, args.bg_index, args.label)


