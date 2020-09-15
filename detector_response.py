"""Build and save the detector response matrix.

Converts between prompt reconstructed energy
and true antineutrino energy.
"""

import argparse
from array import array
import json
import sqlite3
import time


def main(toymc_infile, update_db, image_outfile):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    infile = ROOT.TFile(toymc_infile, 'READ')
    toy_data = infile.Get('toy')
    nbins = (164, 1.8, 12, 27, 1, 12)
    # Bins: 0.25MeV from 1.5 to 8MeV, then 8-12MeV as 1 bin
    hardcoded_reco_bins = array('f', [1.5 + 0.25 * i for i in range(27)] + [12])
    hardcoded_true_bins = array('f', [1.8 + 0.05 * i for i in range(164)] + [12])
    detector_response_hist = ROOT.TH2F('det_response', 'det_response', *nbins)
    detector_response_hist.GetXaxis().Set(nbins[0], hardcoded_true_bins)
    detector_response_hist.GetYaxis().Set(nbins[3], hardcoded_reco_bins)
    if image_outfile:
        canvas = ROOT.TCanvas('c1', 'c1', 800, 800)
        margin = 0.14
        canvas.SetRightMargin(margin)
        canvas.SetLeftMargin(margin)
        canvas.SetTopMargin(margin)
        canvas.SetBottomMargin(margin)
    toy_data.Draw("res_p:Ev >> det_response", "", "colz")
    if image_outfile:
        detector_response_hist.GetXaxis().SetTitle("E_{#nu,true} [MeV]")
        detector_response_hist.GetYaxis().SetTitle("E_{p,reco} [MeV]")
        detector_response_hist.SetLabelSize(0.05, 'xy')
        detector_response_hist.SetLabelSize(0.05, 'z')
        detector_response_hist.SetTitleSize(0.05, 'z')
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gPad.SetLogz()
        ROOT.gPad.Print(image_outfile)

    # Turn the histogram bins into a JSONify-able list of lists.
    response_matrix = []
    for y_index in range(1, nbins[3]+1):
        row = []
        for x_index in range(1, nbins[0] + 1):
            content = detector_response_hist.GetBinContent(x_index, y_index)
            row.append(content)
        response_matrix.append(row)

    if update_db:
        with sqlite3.Connection(update_db) as conn:
            cursor = conn.cursor()
            cursor.execute('''INSERT OR REPLACE INTO detector_response
                VALUES ("THU ToyMC res_p:Ev No Cuts Better binning", ?, ?, ?)''',
                (
                    json.dumps(response_matrix),
                    json.dumps(hardcoded_reco_bins.tolist()),
                    json.dumps(hardcoded_true_bins.tolist()),
                )
            )




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('--update-db')
    parser.add_argument('--pdf')
    args = parser.parse_args()
    main(args.input, args.update_db, args.pdf)
