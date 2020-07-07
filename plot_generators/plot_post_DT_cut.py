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
    hist = ROOT.TH2F('double_coincidences', 'double_coincidences',
            *(bins*2))

    canvas = ROOT.TCanvas('c1', 'c1', 800, 800)
    margin = 0.12
    canvas.SetRightMargin(margin)
    canvas.SetLeftMargin(margin)
    canvas.SetTopMargin(margin)
    canvas.SetBottomMargin(margin)

    events.Draw('energy[1]:energy[0] >> double_coincidences',
            'multiplicity == 2 && '
            'dr_to_prompt[1] + 1000/600e3 * dt_to_prompt[1] < 800'
            , 'colz')

    hist.GetXaxis().SetTitle('Prompt energy [MeV]')
    hist.GetYaxis().SetTitle('Delayed energy [MeV]')

    if args.pdf:
        canvas.SetLogz()
        canvas.Print(args.outfile + '.pdf')
    canvas.Close()
    if args.root:
        out = ROOT.TFile(args.outfile + '.root', 'RECREATE')
        hist.SetDirectory(out)
        hist.Write()
        out.Close()
