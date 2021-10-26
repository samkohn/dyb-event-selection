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
    dr_bins = [500, 0, 5]
    dt_bins = [750, 0, 1500]
    hist = ROOT.TH2F('dr_vs_dt', 'dr_vs_dt',
            *(dt_bins + dr_bins))

    canvas = ROOT.TCanvas('c1', 'c1', 800, 800)
    margin = 0.12
    canvas.SetRightMargin(margin)
    canvas.SetLeftMargin(margin)
    canvas.SetTopMargin(margin)
    canvas.SetBottomMargin(margin)

    events.Draw('(dr_to_prompt[1]/1e3):(dt_to_prompt[1]/1e3) >> dr_vs_dt',
            'multiplicity == 2', 'colz')

    hist.GetXaxis().SetTitle('Coincidence time [#mus]')
    hist.GetYaxis().SetTitle('Coincidence distance [m]')

    if args.pdf:
        canvas.SetLogz()
        canvas.Print(args.outfile + '.pdf')
    canvas.Close()
    if args.root:
        out = ROOT.TFile(args.outfile + '.root', 'RECREATE')
        hist.SetDirectory(out)
        hist.Write()
        out.Close()
