import argparse
import os
import re

import ROOT

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infilename')
    parser.add_argument('histname')
    parser.add_argument('--ehad', type=int, nargs=2)
    parser.add_argument('-o', '--output')
    parser.add_argument('--logy', action='store_true')
    parser.add_argument('--logz', action='store_true')
    parser.add_argument('--colz', action='store_true')
    parser.add_argument('--rebin', type=int)
    parser.add_argument('--rebin2d', type=int, nargs = 2)
    parser.add_argument('--xrange', type=float, nargs=2)
    parser.add_argument('--yrange', type=float, nargs=2)
    parser.add_argument('--xlabel')
    parser.add_argument('--ylabel')
    parser.add_argument('--smallxnumbers', action='store_true')
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)

    if args.ehad is None:
        eh_ad_name = re.search('EH._AD.', args.infilename).group().replace('_', '-')
    else:
        eh_ad_name = f'EH{args.ehad[0]}-AD{args.ehad[1]}'

    infile = ROOT.TFile(args.infilename, 'READ')
    hist = infile.Get(args.histname)

    if args.rebin is not None:
        hist.Rebin(args.rebin)
    if args.rebin2d is not None:
        hist.Rebin2D(*args.rebin2d)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    hist.SetTitle(args.infilename)
    hist.SetLabelSize(0.05, 'xy')
    hist.SetLabelSize(0.05, 'z')
    hist.SetTitleSize(0.05, 'xyz')
    if args.smallxnumbers:
        hist.SetLabelSize(0.035, 'x')
    #hist.SetTitleOffset(0.95, 'y')
    canvas = ROOT.TCanvas('c1', 'c1', 800, 800)
    margin = 0.14
    canvas.SetRightMargin(margin)
    canvas.SetLeftMargin(margin)
    canvas.SetTopMargin(margin)
    canvas.SetBottomMargin(margin)


    hist.Draw('colz' if args.colz else '')
    text = ROOT.TPaveText(0.58, 0.73, 0.85, 0.82, 'NB NDC')
    text.AddText(eh_ad_name)
    text.SetFillColorAlpha(0, 0.5)
    text.Draw()

    if args.logy:
        canvas.SetLogy()
    if args.logz:
        canvas.SetLogz()

    if args.xlabel is not None:
        hist.GetXaxis().SetTitle(args.xlabel)
    if args.ylabel is not None:
        hist.GetYaxis().SetTitle(args.ylabel)
    if args.xrange is not None:
        hist.GetXaxis().SetRangeUser(*args.xrange)
    if args.yrange is not None:
        hist.GetYaxis().SetRangeUser(*args.yrange)
    if args.output is None:
        output = os.path.splitext(args.infilename)[0] + '.pdf'
    else:
        output = args.output
    canvas.Print(output)
    canvas.Close()
    infile.Close()
