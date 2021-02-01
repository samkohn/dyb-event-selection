"""Compute the change in prompt energy cut efficiency due to oscillation effects.

The MC prompt spectrum is distorted by oscillation.
The fraction of MC events passing the prompt energy cut is computed
for no oscillation and for a variety of theta13 values
and for baselines for each reactor-AD pair.
An appropriate adjustment in detection efficiency should be made
to compensate for these changes:
    N_comparable = N_AD/(1 + correction)
"""

import argparse
from array import array
import itertools
import multiprocessing
import sqlite3

import numpy as np
import tenacity

import common
from prediction import survival_probability, default_osc_params
import prediction as pred

cores = range(1, 7)

NUM_EVENTS = 100_000_000

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


def main(
    database,
    binning_id,
    mc_infile,
    source,
    m2_ee_values,
    sin2_values,
    nominal_eff,
    toymc_extracted,
    update_db,
):
    import ROOT
    if toymc_extracted is None:
        infile = ROOT.TFile(mc_infile, 'READ')
        mc_data = infile.Get('toy')
        mc_data.SetBranchStatus('*', 0)
        mc_data.SetBranchStatus('Ev', 1)
        mc_data.SetBranchStatus('res_p', 1)
        mc_data.SetBranchStatus('target', 1)
        # Extract each event's true antineutrino and prompt reco energy
        events = []
        n_entries = mc_data.GetEntries()
        for i in range(min(n_entries, NUM_EVENTS)):
            mc_data.GetEntry(i)
            if mc_data.target == 1:
                events.append((
                    mc_data.Ev,
                    mc_data.res_p,
                ))
        events = np.array(events)
        np.save('prompt_eff_osc_toymc.npy', events)
        infile.Close()
    else:
        events = np.load(toymc_extracted)
    print("Loaded up event array")
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
        fine_bins = cursor.fetchall()
    fine_bins = array('f', [bin_tuple[0] for bin_tuple in fine_bins])
    n_fine_bins = len(fine_bins) - 1  # bins array includes upper bin edge
    # Coarse bins split into "too small," "accepted," and "too large"
    bins = array('f', [0, fine_bins[0], fine_bins[n_fine_bins], 20])
    nbins = len(bins) - 1
    if nominal_eff is None:
        nominal_values, _ = np.histogram(events[:, 1], bins=bins)
        nominal_eff = nominal_values[1]/sum(nominal_values)
        print("no osc", nominal_eff)
        #nominal_hist = ROOT.TH1F("nominal", "nominal", nbins, bins[0], bins[nbins])
        #nominal_hist.GetXaxis().Set(nbins, bins)
        # load the bin contents into numpy arrays for easy manipulation
        #for i in range(nbins):
            #nominal_hist.SetBinContent(i + 1, nominal_values[i])
    print("Filled nominal histogram")
    if update_db:
        with common.get_db(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                INSERT OR REPLACE
                INTO
                    prompt_eff_osc_corrections
                VALUES
                    ("no osc", 0, 0, 0, 0, 0, ?, ?, 0)
                ''',
                (binning_id, nominal_eff)
            )
    #osc_histograms = {}
    osc_effs = {}
    osc_corrections = {}
    rows = []
    for core, halldet, sin2_2theta13, m2_ee in itertools.product(
        cores, pred.all_ads, sin2_values, m2_ee_values
    ):
        print('starting', core, halldet, sin2_2theta13, m2_ee)
        hall, det = halldet
        baseline = pred.distances[core][f'EH{hall}'][det - 1]  # Lookup in prediction
        theta13 = 0.5 * np.arcsin(np.sqrt(sin2_2theta13))
        weights = survival_probability(
            baseline,
            events[:, 0],
            theta13,
            m2_ee,
            default_osc_params,
        )
        osc_values, _ = np.histogram(events[:, 1], bins=bins, weights=weights)
        name = f'eh{hall}_ad{det}_core{core}_{sin2_2theta13:.4f}'
        #osc_hist = ROOT.TH1F(name, name, nbins, bins[0], bins[nbins])
        #osc_hist.GetXaxis().Set(nbins, bins)
        #osc_histograms[hall, det, core, sin2_2theta13] = osc_hist
        #for i in range(nbins):
            #osc_hist.SetBinContent(i + 1, osc_values[i])
        efficiency = osc_values[1]/sum(osc_values)
        correction = efficiency/nominal_eff - 1
        osc_effs[hall, det, core, sin2_2theta13] = efficiency
        osc_corrections[hall, det, core, sin2_2theta13] = correction
        row = (
            source,
            hall,
            det,
            core,
            sin2_2theta13,
            m2_ee,
            binning_id,
            efficiency,
            correction,
        )
        rows.append(row)
        if not update_db:
            print(row)
    if update_db:
        print('updating db')
        fill_db(database, rows)
    print('finished')
    return


@tenacity.retry(
    reraise=True,
    wait=tenacity.wait_random_exponential(max=60),
    retry=tenacity.retry_if_exception_type(sqlite3.Error),
)
def fill_db(database, rows):
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.executemany('''
            INSERT OR REPLACE
            INTO
                prompt_eff_osc_corrections
            VALUES
                (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''',
            rows
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('binning_id', type=int, help='Used to determine energy bounds')
    parser.add_argument('mc_infile')
    parser.add_argument('source')
    parser.add_argument('--sin2-values', type=float, nargs='+')
    parser.add_argument('--dm2-values', type=float, nargs='+')
    parser.add_argument('--no-osc-eff', type=float,
        help='Use the provided efficiency rather than computing it again'
    )
    parser.add_argument('--toymc-npy')
    parser.add_argument('--update-db', action='store_true')
    parser.add_argument('--multiprocessing', action='store_true')
    args = parser.parse_args()
    if args.multiprocessing:
        with multiprocessing.Pool(10) as pool:
            pool.starmap(main, [(
                args.database,
                args.binning_id,
                args.mc_infile,
                args.source,
                [dm2],
                [sin2],
                args.no_osc_eff,
                args.toymc_npy,
                args.update_db,
            ) for sin2, dm2 in itertools.product(args.sin2_values, args.dm2_values)
            ])

    main(
        args.database,
        args.binning_id,
        args.mc_infile,
        args.source,
        args.dm2_values,
        args.sin2_values,
        args.no_osc_eff,
        args.toymc_npy,
        args.update_db
    )
