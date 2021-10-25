"""Scan acc-subtracted dr efficiency."""

import argparse

import numpy as np

import common

def main(sub_file_template, outfilename):
    import ROOT
    efficiencies = {}
    for hall, det in common.all_ads:
        infilename = sub_file_template.format(site=hall, ad=det)
        sub_file = ROOT.TFile(infilename, 'READ')
        dr_sub = sub_file.Get('dr_sub')
        full_integral = dr_sub.Integral()
        cumulative_sums = [0]
        for i in range(1, dr_sub.GetNbinsX() + 1):
            bin_value = dr_sub.GetBinContent(i)
            cumulative_sums.append(bin_value + cumulative_sums[-1])
        efficiencies[hall, det] = np.array(cumulative_sums) / full_integral
    return efficiencies





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('sub_file_template')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    main(args.sub_file_template, args.output)
