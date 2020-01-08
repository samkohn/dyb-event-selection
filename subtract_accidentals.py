from __future__ import print_function

import math
import argparse
import json

def coinc_rate(rs, rmu, tc):
    rsum = rs + rmu
    sum_exp = math.exp(-rsum*tc)
    mu_exp = math.exp(-rmu*tc)
    term1 = rmu/rsum * (1 - sum_exp) + sum_exp
    term2 = rs/rsum * mu_exp * (1 - sum_exp)
    term3 = rs/(2*rs + rmu) * mu_exp * (1 - math.exp(-(2*rs+rmu)*tc))
    prefactor = rs * rs * tc * math.exp(-rs*tc)
    return prefactor * (term1 + term2 + term3)

def main(datafilename, accfilename, rs, rmu, livetime):
    import ROOT
    base_rate = coinc_rate(rs, rmu, 0.0004)
    datafile = ROOT.TFile(datafilename, 'READ')
    raw_spectrum = ROOT.TH2F('raw', 'raw', 210, 1.5, 12, 210, 1.5, 12)
    ad_events = datafile.Get('ad_events')
    ad_events.Draw('energy:energy_previous_PromptLike >> raw',
            'coincidence_number == 2 && dr_previous_PromptLike < 500 &&'
            '!tag_AnyMuonVeto && (tag_flasher == 0'
            '|| tag_flasher == 2)', 'goff')
    accfile = ROOT.TFile(accfilename, 'READ')
    acc_spectrum = accfile.Get('acc_spectrum')
    num_acc_events = acc_spectrum.GetEntries()
    distance_fails = accfile.Get('distance_cut_fails')
    eps_distance = num_acc_events/(
            distance_fails.GetEntries() + num_acc_events)
    outfile = ROOT.TFile('acc_spectrum_test.root', 'RECREATE')
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
    outfile.Write()
    ROOT.gPad.Print('test.pdf')
    datafile.Close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('datafile')
    parser.add_argument('livetimefile')
    parser.add_argument('accfile')
    parser.add_argument('ad', type=int, choices=[1, 2, 3, 4])
    parser.add_argument('rs', type=int)
    parser.add_argument('rmu', type=int)
    args = parser.parse_args()

    with open(args.livetimefile, 'r') as f:
        stats = json.load(f)
        livetime = stats['usable_livetime'][str(args.ad)]/1e9

    main(args.datafile, args.accfile, args.rs, args.rmu, livetime)
