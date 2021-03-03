"""Print out a variety of accidentals cross checks."""

import argparse
import os

import common

def main(sub_dir_template, label, rates_label):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    bin_2m = 41
    bin_3m = 61
    bin_up_to_5m = 100
    bin_5m_and_up = 101
    bin_6m_and_up = 121
    bin_10m = 200

    # print integrals of dr from 2m to 10m and DT from 3m to 10m
    for site, ad in common.all_ads:
        infilename = os.path.join(
            sub_dir_template.format(site=site, ad=ad),
            f'sub_ad{ad}.root'
        )
        infile = ROOT.TFile(infilename, 'READ')
        dr_sub = infile.Get('dr_sub')
        dr_data = infile.Get('dr_data')
        dr_bg = infile.Get('dr_bg')
        DT_sub = infile.Get('DT_sub')
        DT_data = infile.Get('DT_data')
        DT_bg = infile.Get('DT_bg')
        dr_error = ROOT.Double()
        DT_error = ROOT.Double()
        dr_integral = dr_sub.IntegralAndError(bin_2m, bin_10m, dr_error)
        dr_data_integral = dr_data.Integral(bin_2m, bin_10m)
        dr_bg_integral = dr_bg.Integral(bin_2m, bin_10m)
        dr_rel_diff = dr_bg_integral/dr_data_integral - 1
        dr_fit_result = dr_sub.Fit("pol0", "QN0S", "", 2000, 5000)
        DT_integral = DT_sub.IntegralAndError(bin_3m, bin_10m, DT_error)
        DT_data_integral = DT_data.Integral(bin_3m, bin_10m)
        DT_bg_integral = DT_bg.Integral(bin_3m, bin_10m)
        DT_rel_diff = DT_bg_integral/DT_data_integral - 1
        DT_fit_result = DT_sub.Fit("pol0", "QN0S", "", 3000, 6000)
        print('--------')
        print(f'EH{site}-AD{ad}')
        print(
            f'dr 2m->10m: {dr_integral:.1f} +/- {dr_error:.1f} '
            f'({dr_integral/dr_error:.1f} sigma) '
            f'({100*abs(dr_rel_diff):.3f}% '
            f'{"OVER" if dr_rel_diff > 0 else "UNDER"}-subtraction)'
        )
        print(
            f'dr 2m->5m fit: {dr_fit_result.Parameter(0):.1f} +/- '
            f'{dr_fit_result.ParError(0):.1f} '
            f'({dr_fit_result.Parameter(0)/dr_fit_result.ParError(0):.1f} sigma)'
        )
        print(
            f'DT 3m->10m: {DT_integral:.1f} +/- {DT_error:.1f} '
            f'({DT_integral/DT_error:.1f} sigma) '
            f'({100*abs(DT_rel_diff):.3f}% '
            f'{"OVER" if DT_rel_diff > 0 else "UNDER"}-subtraction)'
        )
        print(
            f'DT 3m->6m fit: {DT_fit_result.Parameter(0):.1f} +/- '
            f'{DT_fit_result.ParError(0):.1f} '
            f'({DT_fit_result.Parameter(0)/DT_fit_result.ParError(0):.1f} sigma)'
        )
        dr_integral = dr_sub.IntegralAndError(bin_5m_and_up, bin_10m, dr_error)
        DT_integral = DT_sub.IntegralAndError(bin_6m_and_up, bin_10m, DT_error)
        print(
            f'dr 5m->10m: {dr_integral:.1f} +/- {dr_error:.1f} '
            f'({dr_integral/dr_error:.1f} sigma)'
        )
        print(
            f'DT 6m->10m: {DT_integral:.1f} +/- {DT_error:.1f} '
            f'({DT_integral/DT_error:.1f} sigma)'
        )
        infile.Close()
    return





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('sub_dir_template')
    parser.add_argument('--label')
    parser.add_argument('--rates-label')
    args = parser.parse_args()
    main(args.sub_dir_template, args.label, args.rates_label)
