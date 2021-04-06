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
    bins = [500, -2.5, 2.5]
    hist = ROOT.TH2F('single_xy', 'single_xy',
            *(bins + bins))
    hist2 = ROOT.TH2F('double_xy', 'double_xy',
            *(bins + bins))

    canvas = ROOT.TCanvas('c1', 'c1', 800, 800)
    margin = 0.12
    canvas.SetRightMargin(margin)
    canvas.SetLeftMargin(margin)
    canvas.SetTopMargin(margin)
    canvas.SetBottomMargin(margin)

    events.Draw('y/1e3:x/1e3 >> single_xy',
            'multiplicity == 1', 'colz')
    hist.GetXaxis().SetTitle('Reconstructed x [m]')
    hist.GetYaxis().SetTitle('Reconstructed y [m]')

    canvas2 = ROOT.TCanvas('c2', 'c2', 800, 800)
    canvas2.SetRightMargin(margin)
    canvas2.SetLeftMargin(margin)
    canvas2.SetTopMargin(margin)
    canvas2.SetBottomMargin(margin)

    events.Draw('y/1e3:x/1e3 >> double_xy',
            'multiplicity == 2', 'colz')
    hist2.GetXaxis().SetTitle('Reconstructed x [m]')
    hist2.GetYaxis().SetTitle('Reconstructed y [m]')


    if args.pdf:
        #canvas.SetLogz()
        canvas.Print(args.outfile + '_singles.pdf')
        canvas2.Print(args.outfile + '_doubles.pdf')
    canvas.Close()
    if args.root:
        out = ROOT.TFile(args.outfile + '.root', 'RECREATE')
        hist.SetDirectory(out)
        hist.Write()
        hist2.SetDirectory(out)
        hist2.Write()
        out.Close()
