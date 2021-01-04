"""Compute and dump the bg normalizations from the MC spec file."""
import argparse

import common

BG_LOOKUP = {
    '11': 'accidental',
    '13': 'li9',
    '15': 'fast-neutron',
    '17': 'amc',
}

BG_ERROR_LOOKUP = {
    '12': 'accidental',
    '14': 'li9',
    '16': 'fast-neutron',
    '18': 'amc',
}

def main(infilename, database, label):

    livetime_index = '2'
    muon_eff_index = '3'
    mult_eff_index = '4'
    bg_values = {}
    bg_errors = {}
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
            if splits[1] in BG_LOOKUP:
                bg_rates = [float(val) for val in splits[2:]]
                bg_values[BG_LOOKUP[splits[1]]] = bg_rates
            if splits[1] in BG_ERROR_LOOKUP:
                bg_error = [float(val) for val in splits[2:]]
                bg_errors[BG_ERROR_LOOKUP[splits[1]]] = bg_error
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        for i, ((hall, det), livetime, mu_eff, mult_eff) in enumerate(zip(
            common.all_ads, livetimes, muon_effs, mult_effs
        )):
            for bg_type in bg_values:
                bg_rate = bg_values[bg_type][i]
                bg_error = bg_errors[bg_type][i]
                cursor.execute('''
                    INSERT OR REPLACE INTO
                        bg_counts
                    VALUES
                        (?, ?, ?, ?, ?, ?)
                    ''',
                    (
                        label,
                        hall,
                        det,
                        bg_type,
                        livetime * bg_rate * mu_eff * mult_eff,
                        livetime * bg_error * mu_eff * mult_eff
                    ),
                )
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Theta13-etc. ToyMC file')
    parser.add_argument('database')
    parser.add_argument('--label')
    args = parser.parse_args()
    main(args.infile, args.database, args.label)


