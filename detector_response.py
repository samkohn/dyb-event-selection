"""Build and save the detector response matrix.

Converts between prompt reconstructed energy
and true antineutrino energy.
"""

import argparse
import json
import time

import numpy as np

import common

def main(toymc_infile, update_db, image_outfile, binning):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    infile = ROOT.TFile(toymc_infile, 'READ')
    toy_data = infile.Get('toy')
    keV = 1000
    hardcoded_true_bins = np.concatenate((np.linspace(1.8*keV, 9.95*keV, 164), [12*keV]))
    if binning == 'nominal' or binning is None:
        # Bins: 1.5, 1.6, 0.2MeV from 1.6 to 8MeV, then 8-12MeV as 1 bin
        hardcoded_reco_bins = np.concatenate(([1.5*keV], np.linspace(1.6*keV, 8*keV, 33), [12*keV]))
        nbins = (164, 1.8*keV, 12*keV, 34, 1.5*keV, 12*keV)
    elif binning == 'rate-only':
        hardcoded_reco_bins = np.array([1.5*keV, 12*keV])
        nbins = (164, 1.8*keV, 12*keV, 1, 1.5*keV, 12*keV)
    elif binning == 'nH modified 1':
        # Bins: 1.6MeV, 0.2MeV till 8MeV, then 8-12MeV as 1 bin
        hardcoded_reco_bins = np.concatenate((np.linspace(1.6*keV, 8*keV, 33), [12*keV]))
        nbins = (164, 1.8*keV, 12*keV, 33, 1.6*keV, 12*keV)
    #hardcoded_reco_bins = array('f', [1.5 + 0.25 * i for i in range(27)] + [12])
    #hardcoded_true_bins = array('f', [1.8 + 0.05 * i for i in range(164)] + [12])
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
    toy_data.Draw("res_p*1000:Ev*1000 >> det_response", "", "colz")
    if image_outfile:
        detector_response_hist.GetXaxis().SetTitle("E_{#nu,true} [keV]")
        detector_response_hist.GetYaxis().SetTitle("E_{p,reco} [keV]")
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
        with common.get_db(update_db) as conn:
            cursor = conn.cursor()
            label = f"THU ToyMC res_p:Ev No Cuts {binning} binning"
            cursor.execute('''INSERT OR REPLACE INTO detector_response
                VALUES (?, ?, ?, ?)''',
                (
                    label,
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
    parser.add_argument('--binning')
    args = parser.parse_args()
    main(args.input, args.update_db, args.pdf, args.binning)
