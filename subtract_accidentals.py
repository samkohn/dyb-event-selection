from __future__ import print_function

import math
import argparse
import json
import logging
import os.path
import sqlite3

import tenacity

import common
from adevent import _EMAX_THU
from delayeds import (_NH_THU_DIST_TIME_CUT_STR as DT_CUT_LITERAL, _NH_THU_DIST_TIME_STR
        as DT_VALUE_LITERAL, _NH_THU_MAX_TIME, _NH_THU_MIN_TIME)

DR_VALUE_LITERAL = '(dr_to_prompt[1])'

def coinc_rate(rs, rmu, tc):
    rsum = rs + rmu
    sum_exp = math.exp(-rsum*tc)
    mu_exp = math.exp(-rmu*tc)
    term1 = rmu/rsum * (1 - sum_exp) + sum_exp
    term2 = rs/rsum * mu_exp * (1 - sum_exp)
    term3 = -rs/(2*rs + rmu) * mu_exp * (1 - math.exp(-(2*rs+rmu)*tc))
    prefactor = rs * rs * tc * math.exp(-rs*tc)
    return prefactor * (term1 + term2 + term3)

@tenacity.retry(
    reraise=True,
    wait=tenacity.wait_random_exponential(max=60),
    retry=tenacity.retry_if_exception_type(sqlite3.Error),
)
def subtract(outfilename, datafilename, accfilename, ad, rs, rmu, livetime,
        acc_rate, run_number, database, label):
    import ROOT
    if acc_rate is None:
        window_size_s = (_NH_THU_MAX_TIME - _NH_THU_MIN_TIME)/1e9
        base_rate = coinc_rate(rs, rmu, window_size_s)
    else:
        base_rate = acc_rate
    if 'adtime' in label:
        DT_VALUE = '(dr_to_prompt_AdTime[1] + 1000/600e3 * dt_to_prompt[1])'
        DT_CUT = '(dr_to_prompt_AdTime[1] + 1000/600e3 * dt_to_prompt[1] < 800)'
        DR_VALUE = '(dr_to_prompt_AdTime[1])'
    else:
        DT_VALUE = DT_VALUE_LITERAL
        DT_CUT = DT_CUT_LITERAL
        DR_VALUE = DR_VALUE_LITERAL
    EMAX_CUT = f'(energy[0] < {_EMAX_THU} && energy[1] < {_EMAX_THU})'
    hist_parameters = (2100, 1.5, 12, 2100, 1.5, 12)
    datafile = ROOT.TFile(datafilename, 'READ')
    raw_spectrum = ROOT.TH2F('raw_spec', 'raw_spec', *hist_parameters)
    raw_spectrum.Sumw2()
    ad_events = datafile.Get('ad_events')
    ad_events.Draw('energy[1]:energy[0] >> raw_spec',
            f'multiplicity == 2 && {DT_CUT} && {DR_VALUE} >= 0 && {EMAX_CUT}',
            'goff')
    accfile = ROOT.TFile(accfilename, 'READ')
    acc_entries = accfile.Get('accidentals')
    acc_spectrum = ROOT.TH2F('acc_spectrum', 'acc_spectrum', *hist_parameters)
    acc_spectrum.Sumw2()
    acc_entries.Draw('energy[1]:energy[0] >> acc_spectrum', '', 'goff')
    acc_entries.Draw('energy[0]:energy[1] >>+acc_spectrum', '', 'goff')
    num_acc_events = acc_spectrum.GetEntries()
    DT_cut_fails = accfile.Get('DT_cut_fails')
    eps_distance = num_acc_events/(
            DT_cut_fails.GetEntries() + num_acc_events)
    outfile = ROOT.TFile(outfilename, 'RECREATE')
    if num_acc_events == 0:
        logging.info('Found run with 0 acc events passing DT cut: Run %d, %s', run_number, accfilename)
        if database is not None:
            with common.get_db(database, timeout=60) as conn:
                c = conn.cursor()
                # RunNo, DetNo, BaseRate, DistanceEff, AccScaleFactor,
                # DistanceCrossCheck, DistanceCrossCheck_error
                if label is None:
                    c.execute('''INSERT OR REPLACE INTO accidental_subtraction
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                    (run_number, ad, base_rate, 0, 0, 0, 0, 0, 0))
                else:
                    c.execute('''INSERT OR REPLACE INTO accidental_subtraction
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NULL, NULL)''',
                    (run_number, ad, label, base_rate, 0, 0, 0, 0, 0, 0))
                    conn.commit()
        outfile.Write()
        datafile.Close()
        return
    final_spectrum = ROOT.TH2F('final', 'final', *hist_parameters)
    final_spectrum.Sumw2()
    datafile.cd()
    final_spectrum.Add(raw_spectrum, acc_spectrum, 1,
        -base_rate*eps_distance*livetime/num_acc_events)
    final_spectrum.Rebin2D(2, 2)
    print(base_rate)
    print(eps_distance)
    print(livetime)
    print(num_acc_events)
    print(raw_spectrum.GetEntries())
    outfile.cd()
    distance_axis_parameters = (100, 0, 5000)
    energy_axis_parameters = (2100, 1.5, 12)
    dr_spectrum_actual = ROOT.TH1F('dr_data', 'dr_data', *distance_axis_parameters)
    dr_spectrum_actual.Sumw2()
    dr_spectrum_bg = ROOT.TH1F('dr_bg', 'dr_bg', *distance_axis_parameters)
    dr_spectrum_bg.Sumw2()
    dr_spectrum_sub = ROOT.TH1F('dr_sub', 'dr_sub', *distance_axis_parameters)
    ed_vs_dr_actual = ROOT.TH2F('ed_dr_data', 'ed_dr_data',
            *distance_axis_parameters, *energy_axis_parameters)
    ep_vs_dr_actual = ROOT.TH2F('ep_dr_data', 'ep_dr_data',
            *distance_axis_parameters, *energy_axis_parameters)
    ep_vs_dr_bg = ROOT.TH2F('ep_dr_bg', 'ep_dr_bg', *distance_axis_parameters,
            *energy_axis_parameters)
    ep_vs_dr_sub = ROOT.TH2F('ep_dr_sub', 'ep_dr_sub', *distance_axis_parameters,
            *energy_axis_parameters)
    ed_vs_dr_bg = ROOT.TH2F('ed_dr_bg', 'ed_dr_bg', *distance_axis_parameters,
            *energy_axis_parameters)
    ed_vs_dr_sub = ROOT.TH2F('ed_dr_sub', 'ed_dr_sub', *distance_axis_parameters,
            *energy_axis_parameters)
    DT_spectrum_actual = ROOT.TH1F('DT_data', 'DT_data', *distance_axis_parameters)
    DT_spectrum_actual.Sumw2()
    DT_spectrum_bg = ROOT.TH1F('DT_bg', 'DT_bg', *distance_axis_parameters)
    DT_spectrum_bg.Sumw2()
    DT_spectrum_sub = ROOT.TH1F('DT_sub', 'DT_sub', *distance_axis_parameters)
    ed_vs_DT_actual = ROOT.TH2F('ed_DT_data', 'ed_DT_data',
            *distance_axis_parameters, *energy_axis_parameters)
    ep_vs_DT_actual = ROOT.TH2F('ep_DT_data', 'ep_DT_data',
            *distance_axis_parameters, *energy_axis_parameters)
    ep_vs_DT_bg = ROOT.TH2F('ep_DT_bg', 'ep_DT_bg', *distance_axis_parameters,
            *energy_axis_parameters)
    ep_vs_DT_sub = ROOT.TH2F('ep_DT_sub', 'ep_DT_sub',
            *distance_axis_parameters, *energy_axis_parameters)
    ed_vs_DT_bg = ROOT.TH2F('ed_DT_bg', 'ed_DT_bg', *distance_axis_parameters,
            *energy_axis_parameters)
    ed_vs_DT_sub = ROOT.TH2F('ed_DT_sub', 'ed_DT_sub',
            *distance_axis_parameters, *energy_axis_parameters)
    ad_events.Draw(f'{DR_VALUE} >> dr_data', 'multiplicity == 2 && '
            f'{EMAX_CUT} && {DR_VALUE} < 5000 && {DR_VALUE} >= 0', 'goff')
    scale_factor = base_rate * eps_distance * livetime / num_acc_events
    bg_pairs = accfile.Get('all_pairs')
    bg_pairs.Draw(f'{DR_VALUE_LITERAL} >> dr_bg', (str(scale_factor) +
            f' * 2 * ({DR_VALUE_LITERAL} < 5000 && ' +
            f'{DR_VALUE_LITERAL} >= 0)'), 'same goff')
    dr_spectrum_sub.Add(dr_spectrum_actual, dr_spectrum_bg, 1, -1)
    ad_events.Draw(f'{DT_VALUE} >> DT_data',
        f'multiplicity == 2 && {EMAX_CUT} && {DT_VALUE} < 5000 && {DR_VALUE} >= 0',
        'goff')
    bg_pairs.Draw(f'{DT_VALUE_LITERAL} >> DT_bg',
        f'{scale_factor} * 2 * ({DT_VALUE_LITERAL} < 5000 && {DR_VALUE_LITERAL} >= 0)', 'goff')
    DT_spectrum_sub.Add(DT_spectrum_actual, DT_spectrum_bg, 1, -1)
    if database is not None:
        # Then do the fit to get the parameters for the database
        fit_result = dr_spectrum_sub.Fit('pol0', 'QN0S', '', 2000, 5000)
        distance_cross_check_value = fit_result.Parameter(0)
        distance_cross_check_error = fit_result.ParError(0)
        fit_result = DT_spectrum_sub.Fit('pol0', 'QN0S', '', 2000, 5000)
        DT_cross_check_value = fit_result.Parameter(0)
        DT_cross_check_error = fit_result.ParError(0)
    ad_events.Draw(f'energy[1]:{DR_VALUE} >> ed_dr_data',
            f'multiplicity == 2 && {EMAX_CUT} && {DR_VALUE} < 5000 && {DR_VALUE} >= 0',
            'goff')
    ad_events.Draw(f'energy[0]:{DR_VALUE} >> ep_dr_data',
            f'multiplicity == 2 && {EMAX_CUT} && {DR_VALUE} < 5000 && {DR_VALUE} >= 0',
            'goff')
    bg_pairs.Draw(f'energy[0]:{DR_VALUE_LITERAL} >> ed_dr_bg',
            f'2 * ({DR_VALUE_LITERAL} < 5000 && multiplicity == 2 && '
            f'{DR_VALUE_LITERAL} >= 0)', 'goff')
    ep_vs_dr_sub.Add(ep_vs_dr_actual, ep_vs_dr_bg, 1,
            -base_rate*eps_distance*livetime/num_acc_events)
    bg_pairs.Draw(f'energy[1]:{DR_VALUE_LITERAL} >> ed_dr_bg',
            f'2 * ({DR_VALUE_LITERAL} < 5000 && multiplicity == 2 && '
            f'{DR_VALUE_LITERAL} >= 0)', 'goff')
    ed_vs_dr_sub.Add(ed_vs_dr_actual, ed_vs_dr_bg, 1,
            -base_rate*eps_distance*livetime/num_acc_events)
    ad_events.Draw(f'energy[1]:{DT_VALUE} >> ed_DT_data',
            f'multiplicity == 2 && {EMAX_CUT} && {DT_VALUE} < 5000 && {DR_VALUE} >= 0',
            'goff')
    ad_events.Draw(f'energy[0]:{DT_VALUE} >> ep_DT_data',
            f'multiplicity == 2 && {EMAX_CUT} && {DT_VALUE} < 5000 && {DR_VALUE} >= 0',
            'goff')
    bg_pairs.Draw(f'energy[0]:{DT_VALUE_LITERAL} >> ed_DT_bg',
            f'2 * ({DT_VALUE_LITERAL} < 5000 && multiplicity == 2 && '
            f'{DR_VALUE_LITERAL} >= 0)', 'goff')
    ep_vs_DT_sub.Add(ep_vs_DT_actual, ep_vs_DT_bg, 1,
            -base_rate*eps_distance*livetime/num_acc_events)
    bg_pairs.Draw(f'energy[1]:{DT_VALUE_LITERAL} >> ed_DT_bg',
            f'2 * ({DT_VALUE_LITERAL} < 5000 && multiplicity == 2 && '
            f'{DR_VALUE_LITERAL} >= 0)', 'goff')
    ed_vs_DT_sub.Add(ed_vs_DT_actual, ed_vs_DT_bg, 1,
            -base_rate*eps_distance*livetime/num_acc_events)
    outfile.Write()
    datafile.Close()
    if database is not None:
        with common.get_db(database) as conn:
            c = conn.cursor()
            # RunNo, DetNo, BaseRate, DistanceEff, AccScaleFactor,
            # DistanceCrossCheck, DistanceCrossCheck_error
            if label is None:
                c.execute('''INSERT OR REPLACE INTO accidental_subtraction
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                (run_number, ad, base_rate, eps_distance, scale_factor,
                    DT_cross_check_value, DT_cross_check_error,
                    distance_cross_check_value, distance_cross_check_error))
            else:
                c.execute('''INSERT OR REPLACE INTO accidental_subtraction
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NULL, NULL)''',
                (run_number, ad, label, base_rate, eps_distance, scale_factor,
                    DT_cross_check_value, DT_cross_check_error,
                    distance_cross_check_value, distance_cross_check_error))
                conn.commit()


def main(output, datafile, accfile, database, ad, override_acc_rate, label, update_db):
    try:
        with open(os.path.splitext(datafile)[0] + '.json', 'r') as f:
            stats = json.load(f)
            # livetime = stats['usable_livetime']/1e9
            run_number = stats['run']
            site = stats['site']
    except FileNotFoundError:
        import ROOT
        infile = ROOT.TFile(datafile, 'READ')
        ad_events = infile.Get('ad_events')
        ad_events.GetEntry(0)
        run_number = ad_events.run
        site = ad_events.site

    with common.get_db(database) as conn:
        c = conn.cursor()
        if override_acc_rate:
            singles_rate = None
        else:
            c.execute('SELECT Rate_Hz FROM singles_rates WHERE '
                    'RunNo = ? AND DetNo = ?', (run_number, ad))
            singles_rate, = c.fetchone()
        c.execute('SELECT Rate_Hz, Livetime_ns/1e9 FROM muon_rates WHERE '
                'RunNo = ? AND DetNo = ?', (run_number, ad))
        muon_rate, livetime = c.fetchone()
    database = database if update_db else None
    subtract(output, datafile, accfile, ad, singles_rate, muon_rate, livetime,
            override_acc_rate, run_number, database, label)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('datafile')
    parser.add_argument('accfile')
    parser.add_argument('database')
    parser.add_argument('ad', type=int, choices=[1, 2, 3, 4])
    parser.add_argument('--override-acc-rate', type=float, default=None)
    parser.add_argument('-o', '--output')
    parser.add_argument('--update-db', action='store_true')
    parser.add_argument('--label',
        help='Label in DB - if absent, will assume no Label column in schema')
    args = parser.parse_args()
    main(
        args.output,
        args.datafile,
        args.accfile,
        args.database,
        args.ad,
        args.override_acc_rate,
        args.label,
        args.update_db,
    )
