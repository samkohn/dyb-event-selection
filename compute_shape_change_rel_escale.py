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

import numpy as np

import common
from prediction import survival_probability, default_osc_params

# Approximate effective baselines for each hall
effective_baselines = {
    1: 365,
    2: 500,
    3: 1540,
}

def get_p_sur_helpers(sin2_2theta13, m2_ee):
    """Return a bunch of expressions to make the TTree::Draw expression sane."""
    theta13 = 0.5 * np.arcsin(np.sqrt(sin2_2theta13))
    theta12 = default_osc_params.theta12
    m2_21 = default_osc_params.m2_21
    # Hierarchy factor = 2 * int(hierarchy) - 1 (+1 if True, else -1)
    hierarchy = 2 * int(default_osc_params.hierarchy) - 1
    m2_32 = m2_ee - hierarchy * default_osc_params.m2_ee_conversion
    m2_31 = m2_32 + hierarchy * m2_21
    cos4theta13 = np.power(np.cos(theta13), 4)
    cos2theta12 = np.power(np.cos(theta12), 2)
    sin2theta12 = np.power(np.sin(theta12), 2)
    sin2_2theta12 = np.power(np.sin(2*theta12), 2)
    return (
        cos4theta13,
        cos2theta12,
        sin2theta12,
        sin2_2theta12,
        m2_21,
        m2_31,
        m2_32,
    )


def main(database, binning_id, mc_infile, source, rel_uncertainty, sin2_2theta13, m2_ee,
        site, update_db):
    import ROOT
    infile = ROOT.TFile(mc_infile, 'READ')
    mc_data = infile.Get('toy')
    # Read in the prompt energy bins
    with common.get_db(database) as conn:
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
    # Prepare for survival probability computation
    (
        cos4theta13,
        cos2theta12,
        sin2theta12,
        sin2_2theta12,
        m2_21,
        m2_31,
        m2_32,
    ) = get_p_sur_helpers(sin2_2theta13, m2_ee)
    baseline = effective_baselines[site]
    energy_var = 'res_p'
    p_sur_string = (
        '1'
        f'- {cos4theta13 * sin2_2theta12} * TMath::Power('
        f'    TMath::Sin({1.267 * m2_21 * baseline} / {energy_var}),'
        '     2'
        ')'
        f'- {cos2theta12 * sin2_2theta13} * TMath::Power('
        f'    TMath::Sin({1.267 * m2_31 * baseline} / {energy_var}),'
        '     2'
        ')'
        f'- {sin2theta12 * sin2_2theta13} * TMath::Power('
        f'    TMath::Sin({1.267 * m2_32 * baseline} / {energy_var}),'
        '     2'
        ')'
    )
    selection_string = f'{p_sur_string} * (target == 1)'
    mc_data.Draw("res_p >> nominal", selection_string, "goff")
    mc_data.Draw(f"res_p * (1 + {rel_uncertainty}) >> escaled_plus", selection_string,
        "goff")
    mc_data.Draw(f"res_p * (1 - {rel_uncertainty}) >> escaled_minus", selection_string,
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
            (source, site, sin2_2theta13, m2_ee, binning_id, i, plus_coeffs[i], minus_coeffs[i], rel_uncertainty)
        )
    if update_db:
        with common.get_db(database) as conn:
            cursor = conn.cursor()
            cursor.executemany('''
                INSERT OR REPLACE
                INTO
                    rel_energy_scale_shape
                VALUES
                    (?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''',
                rows
            )
    else:
        print(rows)







if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('binning_id')
    parser.add_argument('mc_infile')
    parser.add_argument('source')
    parser.add_argument('--rel-uncertainty', type=float, default=0.005)
    parser.add_argument('--sin2', type=float, required=True)
    parser.add_argument('--dm2', type=float, required=True)
    parser.add_argument('--site', type=int, required=True)
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(
        args.database,
        args.binning_id,
        args.mc_infile,
        args.source,
        args.rel_uncertainty,
        args.sin2,
        args.dm2,
        args.site,
        args.update_db
    )
