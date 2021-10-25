"""Compute and dump bg counts given a JSON file with the rates.

The JSON file format should be:

{
  "bg-name1": {
    "EH1-AD1": {
      "rate": rate-per-AD-per-day,
      "error": error-on-rate
    },
    "EH1-AD2": {
      ... etc
    }
  },
  "bg-name2": {
    ... etc
  }
}
"""
import argparse
import json

import numpy as np

import common
import compute_dataset_summary as summary

def parse_hall_det(eh_ad_str):
    """Convert "EHX-ADY" into (X, Y) (as ints)."""
    return (int(eh_ad_str[2]), int(eh_ad_str[6]))

def main(infilename, bg_database, main_database, label, rates_label):
    # Get muon and multiplicity efficiencies and livetimes
    daq_livetime_days = summary.daq_livetime_days(main_database, rates_label)
    muon_effs = summary.muon_efficiency(main_database, rates_label)
    mult_effs = summary.multiplicity_efficiency(main_database, rates_label)

    with open(infilename, 'r') as infile:
        bgs_rates = json.load(infile)
    # Convert verbose JSON format into a more-useful array/list format.
    bg_values = {}
    bg_errors = {}
    for bg_type, bg_dict in bgs_rates.items():
        bg_values[bg_type] = []
        bg_errors[bg_type] = []
        for (hall, det) in common.all_ads:
            bg = bg_dict[f'EH{hall}-AD{det}']
            bg_values[bg_type].append(bg['rate'])
            bg_errors[bg_type].append(bg['error'])

    with common.get_db(bg_database) as conn:
        cursor = conn.cursor()
        for i, ((hall, det), livetime, mu_eff, mult_eff) in enumerate(zip(
            common.all_ads, daq_livetime_days, muon_effs, mult_effs
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
    parser.add_argument('infile', help='JSON file with rates & errors')
    parser.add_argument('bg_database', help='Destination')
    parser.add_argument('main_database', help='Source for livetimes & effs')
    parser.add_argument('--label')
    parser.add_argument('--rates-label',
        help='Label for lookup of singles and muon rates'
    )
    args = parser.parse_args()
    main(args.infile, args.bg_database, args.main_database, args.label, args.rates_label)
