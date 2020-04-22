from __future__ import print_function

import math
import argparse
import json
import os.path
import sqlite3

def coinc_rate(rs, rmu, tc):
    rsum = rs + rmu
    sum_exp = math.exp(-rsum*tc)
    mu_exp = math.exp(-rmu*tc)
    term1 = rmu/rsum * (1 - sum_exp) + sum_exp
    term2 = rs/rsum * mu_exp * (1 - sum_exp)
    term3 = -rs/(2*rs + rmu) * mu_exp * (1 - math.exp(-(2*rs+rmu)*tc))
    prefactor = rs * rs * tc * math.exp(-rs*tc)
    return prefactor * (term1 + term2 + term3)

def main(outfilename, datafilename, accfilename, ad, rs, rmu, livetime,
        acc_rate, run_number, database):
    import ROOT
    if acc_rate is None:
        base_rate = coinc_rate(rs, rmu, 0.000399)
    else:
        base_rate = acc_rate
    hist_parameters = (2100, 1.5, 12, 2100, 1.5, 12)
    datafile = ROOT.TFile(datafilename, 'READ')
    raw_spectrum = ROOT.TH2F('raw_spec', 'raw_spec', *hist_parameters)
    raw_spectrum.Sumw2()
    ad_events = datafile.Get('ad_events')
    ad_events.Draw('energy[1]:energy[0] >> raw_spec',
            'multiplicity == 2 && dr_to_prompt[1] < 500',
            'goff')
    accfile = ROOT.TFile(accfilename, 'READ')
    acc_entries = accfile.Get('accidentals')
    acc_spectrum = ROOT.TH2F('acc_spectrum', 'acc_spectrum', *hist_parameters)
    acc_spectrum.Sumw2()
    acc_entries.Draw('energy[1]:energy[0] >> acc_spectrum', '', 'goff')
    acc_entries.Draw('energy[0]:energy[1] >>+acc_spectrum', '', 'goff')
    num_acc_events = acc_spectrum.GetEntries()
    distance_fails = accfile.Get('distance_cut_fails')
    eps_distance = num_acc_events/(
            distance_fails.GetEntries() + num_acc_events)
    outfile = ROOT.TFile(outfilename, 'RECREATE')
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
    dr_spectrum_actual = ROOT.TH1F('dr_data', 'dr_data', 100, 0, 5000)
    dr_spectrum_actual.Sumw2()
    dr_spectrum_bg = ROOT.TH1F('dr_bg', 'dr_bg', 100, 0, 5000)
    dr_spectrum_bg.Sumw2()
    dr_spectrum_sub = ROOT.TH1F('dr_sub', 'dr_sub', 100, 0, 5000)
    ed_vs_dr_actual = ROOT.TH2F('ed_dr_data', 'ed_dr_data', 100, 0, 5000, 210,
            1.5, 12)
    ep_vs_dr_actual = ROOT.TH2F('ep_dr_data', 'ep_dr_data', 100, 0, 5000, 210,
            1.5, 12)
    ep_vs_dr_bg = ROOT.TH2F('ep_dr_bg', 'ep_dr_bg', 100, 0, 5000, 210,
            1.5, 12)
    ep_vs_dr_sub = ROOT.TH2F('ep_dr_sub', 'ep_dr_sub', 100, 0, 5000, 210, 1.5,
            12)
    ed_vs_dr_bg = ROOT.TH2F('ed_dr_bg', 'ed_dr_bg', 100, 0, 5000, 210,
            1.5, 12)
    ed_vs_dr_sub = ROOT.TH2F('ed_dr_sub', 'ed_dr_sub', 100, 0, 5000, 210, 1.5,
            12)
    ad_events.Draw('dr_to_prompt[1] >> dr_data', 'multiplicity == 2 && '
            'dr_to_prompt[1] < 5000', 'goff')
    scale_factor = base_rate * eps_distance * livetime / num_acc_events
    bg_pairs = accfile.Get('all_pairs')
    bg_pairs.Draw('dr_to_prompt[1] >> dr_bg', (str(scale_factor) +
            ' * 2 * (dr_to_prompt[1] < 5000 && ' +
            'dr_to_prompt[1] >= 0)'), 'same goff')
    dr_spectrum_sub.Add(dr_spectrum_actual, dr_spectrum_bg, 1, -1)
    if database is not None:
        # Then do the fit to get the parameters for the database
        fit_result = dr_spectrum_sub.Fit('pol0', 'QN0S', '', 2000, 5000)
        cross_check_value = fit_result.Parameter(0)
        cross_check_error = fit_result.ParError(0)
    ad_events.Draw('energy[1]:dr_to_prompt[1] >> ed_dr_data',
            'multiplicity == 2 && dr_to_prompt[1] < 5000 && dr_to_prompt[1] >= 0',
            'goff')
    ad_events.Draw('energy[0]:dr_to_prompt[1] >> ep_dr_data',
            'multiplicity == 2 && dr_to_prompt[1] < 5000 && dr_to_prompt[1] >= 0',
            'goff')
    bg_pairs.Draw('energy[0]:dr_to_prompt[1] >> ed_dr_bg',
            '2 * (dr_to_prompt[1] < 5000 && multiplicity == 2 && '
            'dr_to_prompt[1] >= 0)', 'goff')
    ep_vs_dr_sub.Add(ep_vs_dr_actual, ep_vs_dr_bg, 1,
            -base_rate*eps_distance*livetime/num_acc_events)
    bg_pairs.Draw('energy[1]:dr_to_prompt[1] >> ed_dr_bg',
            '2 * (dr_to_prompt[1] < 5000 && multiplicity == 2 && '
            'dr_to_prompt[1] >= 0)', 'goff')
    ed_vs_dr_sub.Add(ed_vs_dr_actual, ed_vs_dr_bg, 1,
            -base_rate*eps_distance*livetime/num_acc_events)
    outfile.Write()
    datafile.Close()
    if database is not None:
        with sqlite3.Connection(database) as conn:
            c = conn.cursor()
            # RunNo, DetNo, BaseRate, DistanceEff, AccScaleFactor,
            # DistanceCrossCheck, DistanceCrossCheck_error
            c.execute('''INSERT OR REPLACE INTO accidental_subtraction
            VALUES (?, ?, ?, ?, ?, ?, ?)''',
            (run_number, ad, base_rate, eps_distance, scale_factor,
                cross_check_value, cross_check_error))
            conn.commit()




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('datafile')
    parser.add_argument('accfile')
    parser.add_argument('database')
    parser.add_argument('ad', type=int, choices=[1, 2, 3, 4])
    parser.add_argument('--override-acc-rate', type=float, default=None)
    parser.add_argument('-o', '--output')
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()

    with open(os.path.splitext(args.datafile)[0] + '.json', 'r') as f:
        stats = json.load(f)
        livetime = stats['usable_livetime']/1e9
        run_number = stats['run']
        site = stats['site']

    with sqlite3.Connection(args.database) as conn:
        c = conn.cursor()
        if args.override_acc_rate:
            singles_rate = None
        else:
            c.execute('SELECT Rate_Hz FROM singles_rates WHERE '
                    'RunNo = ? AND DetNo = ?', (run_number, args.ad))
            singles_rate, = c.fetchone()
        c.execute('SELECT Rate_Hz FROM muon_rates WHERE '
                'RunNo = ? AND DetNo = ?', (run_number, args.ad))
        muon_rate, = c.fetchone()

    database = args.database if args.update_db else None
    main(args.output, args.datafile, args.accfile, args.ad, singles_rate, muon_rate, livetime,
            args.override_acc_rate, run_number, database)
