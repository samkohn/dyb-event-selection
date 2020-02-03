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
    term3 = rs/(2*rs + rmu) * mu_exp * (1 - math.exp(-(2*rs+rmu)*tc))
    prefactor = rs * rs * tc * math.exp(-rs*tc)
    return prefactor * (term1 + term2 + term3)

def main(datafilename, accfilename, ad, rs, rmu, livetime, acc_rate):
    import ROOT
    if acc_rate is None:
        base_rate = coinc_rate(rs, rmu, 0.0004)
    else:
        base_rate = acc_rate
    datafile = ROOT.TFile(datafilename, 'READ')
    raw_spectrum = ROOT.TH2F('raw_spec', 'raw_spec', 210, 1.5, 12, 210, 1.5, 12)
    ad_events = datafile.Get('ad_events')
    ad_events.Draw('energy[1]:energy[0] >> raw_spec',
            'dr_to_prompt[1] < 500',
            'goff')
    accfile = ROOT.TFile(accfilename, 'READ')
    acc_spectrum = accfile.Get('acc_spectrum')
    num_acc_events = acc_spectrum.GetEntries()
    distance_fails = accfile.Get('distance_cut_fails')
    eps_distance = num_acc_events/(
            distance_fails.GetEntries() + num_acc_events)
    outfile = ROOT.TFile('acc_spectrum_test22011.root', 'RECREATE')
    final_spectrum = ROOT.TH2F('final', 'final', 210, 1.5, 12, 210, 1.5, 12)
    datafile.cd()
    final_spectrum.Add(raw_spectrum, acc_spectrum, 1,
        -base_rate*eps_distance*livetime/num_acc_events)
    final_spectrum.Rebin2D(2, 2)
    print(base_rate)
    print(eps_distance)
    print(livetime)
    print(num_acc_events)
    print(raw_spectrum.GetEntries())
    final_spectrum.Draw('colz')
    outfile.cd()
    dr_spectrum_actual = ROOT.TH1F('dr_data', 'dr_data', 100, 0, 5000)
    dr_spectrum_bg = ROOT.TH1F('dr_bg', 'dr_bg', 100, 0, 5000)
    ed_vs_dr_actual = ROOT.TH2F('ed_dr_data', 'ed_dr_data', 100, 0, 5000, 210,
            1.5, 12)
    ep_vs_dr_actual = ROOT.TH2F('ep_dr_data', 'ep_dr_data', 100, 0, 5000, 210,
            1.5, 12)
    ed_vs_dr_bg = ROOT.TH2F('ed_dr_bg', 'ed_dr_bg', 100, 0, 5000, 210,
            1.5, 12)
    ad_events.Draw('dr_to_prompt[1] >> dr_data', 'dr_to_prompt[1] < 5000')
    scale_factor = base_rate * eps_distance * livetime / num_acc_events
    bg_pairs = accfile.Get('all_pairs')
    bg_pairs.Draw('dr_to_prompt[1] >> dr_bg', (str(scale_factor) +
            ' * 2 * (dr_to_prompt[1] < 5000 && ' +
            'dr_to_prompt[1] >= 0)'), 'same')
    ROOT.gPad.Print('test.pdf')
    ad_events.Draw('energy[1]:dr_to_prompt[1] >> ed_dr_data',
            'dr_to_prompt[1] < 5000 && dr_to_prompt[1] >= 0')
    ad_events.Draw('energy[0]:dr_to_prompt[1] >> ep_dr_data',
            'dr_to_prompt[1] < 5000 && dr_to_prompt[1] >= 0')
    bg_pairs.Draw('energy[1]:dr_to_prompt[1] >> ed_dr_bg', (str(scale_factor) +
            ' * 2 * (dr_to_prompt[1] < 5000)'))
    outfile.Write()
    datafile.Close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('datafile')
    parser.add_argument('accfile')
    parser.add_argument('database')
    parser.add_argument('ad', type=int, choices=[1, 2, 3, 4])
    parser.add_argument('--override-acc-rate', type=float, default=None)
    args = parser.parse_args()

    with open(os.path.splitext(args.datafile)[0] + '.json', 'r') as f:
        stats = json.load(f)
        livetime = stats['usable_livetime']/1e9
        run_number = stats['run']

    with sqlite3.Connection(args.database) as conn:
        c = conn.cursor()
        c.execute('SELECT Rate_Hz FROM singles_rates WHERE '
                'RunNo = ? AND DetNo = ?', (run_number, args.ad))
        singles_rate, = c.fetchone()
        #c.execute('SELECT Rate_Hz FROM muon_rates WHERE '
                #'RunNo = ? AND DetNo = ?', (run_number, args.ad))
        #muon_rate, = c.fetchone()
        muon_rate = 214

    main(args.datafile, args.accfile, args.ad, singles_rate, muon_rate, livetime,
            args.override_acc_rate)
