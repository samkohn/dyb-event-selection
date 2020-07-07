import argparse
from glob import iglob

import ROOT

if __name__ == '__main__':
    ROOT.gROOT.SetBatch(True)
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('--pdf', action='store_true')
    parser.add_argument('--root', action='store_true')
    args = parser.parse_args()

    events = ROOT.TChain('ad_events')

    infiles = iglob(args.infile)
    for infile in infiles:
        events.Add(infile)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    bins = [1050, 1.5, 12]
    hist = ROOT.TH1F('singles', 'singles',
            *(bins))

    canvas = ROOT.TCanvas('c1', 'c1', 800, 800)
    margin = 0.12
    canvas.SetRightMargin(margin)
    canvas.SetLeftMargin(margin)
    canvas.SetTopMargin(margin)
    canvas.SetBottomMargin(margin)

    events.Draw('energy[0] >> singles',
            'multiplicity == 1 && '
            'dt_cluster_to_prev_ADevent > 1500e3')

    hist.GetXaxis().SetTitle('Energy [MeV]')
    hist.GetYaxis().SetTitle('Number of events per 0.01 MeV')

    if args.pdf:
        canvas.SetLogy()
        canvas.Print(args.outfile + '.pdf')
    canvas.Close()
    if args.root:
        out = ROOT.TFile(args.outfile + '.root', 'RECREATE')
        hist.SetDirectory(out)
        hist.Write()
        out.Close()
