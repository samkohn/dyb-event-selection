"""Compute the bin-by-bin shape change due to relative energy scale uncertainty.

The MC prompt spectrum is distorted by +/- the relative uncertainty.
The fractional change in the number of events for each prompt energy bin
(specified by the binning_id)
is saved to the database.
Each bin's value can be scaled (correlated-ly) by the pull parameter
for a given AD to get the approximate impact of a change in energy scale
by the given fraction of the nominal uncertainty.
"""

import argparse
from array import array
import sqlite3

import numpy as np

def main(database, binning_id, mc_infile, source, rel_uncertainty, update_db):
    import ROOT
    infile = ROOT.TFile(mc_infile, 'READ')
    mc_data = infile.Get('toy')
    # Read in the prompt energy bins
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                BinEdgeEnergy_keV/1000.0
            FROM
                reco_binnings
            WHERE
                Id = ?
            ORDER BY
                BinEdgeIndex
            ''',
            (binning_id,)
        )
        bins = cursor.fetchall()
    bins = array('f', [bin_tuple[0] for bin_tuple in bins])
    nbins = len(bins) - 1  # bins array includes upper bin edge
    nominal_histogram = ROOT.TH1F('nominal', 'nominal', nbins, bins[0], bins[nbins])
    nominal_histogram.GetXaxis().Set(nbins, bins)
    escaled_plus_histogram = ROOT.TH1F(
        'escaled_plus', 'escaled_plus', nbins, bins[0], bins[nbins]
    )
    escaled_plus_histogram.GetXaxis().Set(nbins, bins)
    escaled_minus_histogram = ROOT.TH1F(
        'escaled_minus', 'escaled_minus', nbins, bins[0], bins[nbins]
    )
    escaled_minus_histogram.GetXaxis().Set(nbins, bins)
    mc_data.Draw("res_p >> nominal", "target == 1", "goff")
    mc_data.Draw(f"res_p * (1 + {rel_uncertainty}) >> escaled_plus", "target == 1", "goff")
    mc_data.Draw(f"res_p * (1 - {rel_uncertainty}) >> escaled_minus", "target == 1",
            "goff")
    # Load the bin contents into numpy arrays for easy manipulation
    nominal_values = np.zeros((nbins,))
    escaled_plus_values = np.zeros_like(nominal_values)
    escaled_minus_values = np.zeros_like(nominal_values)
    for i in range(nbins):
        nominal_values[i] = nominal_histogram.GetBinContent(i + 1)
        escaled_plus_values[i] = escaled_plus_histogram.GetBinContent(i + 1)
        escaled_minus_values[i] = escaled_minus_histogram.GetBinContent(i + 1)
    plus_coeffs = escaled_plus_values/nominal_values - 1
    minus_coeffs = escaled_minus_values/nominal_values - 1
    rows = []
    for i in range(nbins):
        rows.append(
            (source, binning_id, i, plus_coeffs[i], minus_coeffs[i], rel_uncertainty)
        )
    if update_db:
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.executemany('''
                INSERT OR REPLACE
                INTO
                    rel_energy_scale_shape
                VALUES
                    (?, ?, ?, ?, ?, ?)
                ''',
                rows
            )
    else:
        print(np.array(rows)[:, 2:])







if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('binning_id')
    parser.add_argument('mc_infile')
    parser.add_argument('source')
    parser.add_argument('--rel-uncertainty', type=float, default=0.005)
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(
        args.database,
        args.binning_id,
        args.mc_infile,
        args.source,
        args.rel_uncertainty,
        args.update_db
    )
