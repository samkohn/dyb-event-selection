"""Compute the li9 rate given the fit output and all the various efficiencies.

"""
import argparse
from math import exp
import sqlite3

from li9_common import TAU_LI9, TAU_HE8, R_HE8

def main(database, site, label, ntag):
    params = retrieve_db_params(database, site, label, ntag)

    extension_pre_20ms_factor = (
            R_HE8 / exp(-0.02/TAU_LI9) + (1 - R_HE8) / exp(-0.02/TAU_HE8)
    )
    extension_pre_20ms = params['N_Li9_He8'] * extension_pre_20ms_factor
    err_extension = params['N_Li9_He8_error'] * extension_pre_20ms_factor


    mult_cut_eff = 0.95  # TODO
    DT_cut_eff = 0.7  # TODO
    prompt_energy_eff = 0.7834  # Source: Chris's spectrum bins 15->48 / 7->48
    WSMuon_veto_eff = params['VetoEfficiency']
    if ntag:
        ntag_eff = 0.6  # TODO propagate uncertainty
    else:
        ntag_eff = 1

    efficiency_corrected = (
            extension_pre_20ms
            / mult_cut_eff
            / DT_cut_eff
            / prompt_energy_eff
            / WSMuon_veto_eff
            / ntag_eff
    )
    efficiency_corrected_err = (
            err_extension
            / mult_cut_eff
            / DT_cut_eff
            / prompt_energy_eff
            / WSMuon_veto_eff
            / ntag_eff
    )

    S_PER_DAY = 60*60*24
    rate_preveto = efficiency_corrected / (params['DAQLivetime_ns']/1e9) * S_PER_DAY
    rate_preveto_err = efficiency_corrected_err / (params['DAQLivetime_ns']/1e9) * S_PER_DAY

    veto_corrections = {
            'low': exp(-0.0008/TAU_LI9),
            'mid': exp(-0.0008/TAU_LI9),
            'high': exp(-1/TAU_LI9),
    }
    rate_postveto = rate_preveto * veto_corrections[label]
    rate_postveto_err = rate_preveto_err * veto_corrections[label]
    print(f'Rate pre-veto = {rate_preveto} +/- {rate_preveto_err} per AD per day')
    print(f'Final rate = {rate_postveto} +/- {rate_postveto_err} per AD per day')


def retrieve_db_params(database, eh, label, ntag):
    with sqlite3.Connection(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('SELECT N_Li9_He8, N_Li9_He8_error, DAQLivetime_ns, VetoEfficiency '
                'FROM li9_fits INNER JOIN li9_muon_rates '
                'USING (Hall, EnergyClass, NeutronTag) '
                'WHERE Hall = ? AND EnergyClass = ? '
                'AND NeutronTag = ?', (eh, label, ntag))
        return cursor.fetchone()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('--site', choices=[1, 2, 3], type=int)
    parser.add_argument('--energy', choices=['low', 'mid', 'high'])
    parser.add_argument('--ntag', action='store_true')
    args = parser.parse_args()
    main(args.database, args.site, args.energy, args.ntag)
