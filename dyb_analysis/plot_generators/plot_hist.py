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
    parser.add_argument('--invert', action='store_true', help='Reciprocal of all bins')
    parser.add_argument('--mm-to-m-hack', action='store_true',
        help='Xbins in mm, convert to m')
    parser.add_argument('--errors', action='store_true',
        help='If --mm-to-m-hack, also transfer bin errors')
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

    if args.invert:
        for binx in range(1, hist.GetNbinsX() + 1):
            for biny in range(1, hist.GetNbinsY() + 1):
                for binz in range(1, hist.GetNbinsZ() + 1):
                    global_bin = hist.GetBin(binx, biny, binz)
                    value = hist.GetBinContent(global_bin)
                    hist.SetBinContent(global_bin, 1/value)

    if args.mm_to_m_hack:
        xaxis = hist.GetXaxis()
        nbins = xaxis.GetNbins()
        xmin = xaxis.GetBinLowEdge(1)
        xmax = xaxis.GetBinLowEdge(nbins + 1)
        new_hist = ROOT.TH1F("new_hist", "new_hist", nbins, xmin / 1000, xmax / 1000)
        for bin_index in range(1, nbins+1):
            bin_content = hist.GetBinContent(bin_index)
            new_hist.SetBinContent(bin_index, bin_content)
        if args.errors:
            for bin_index in range(1, nbins+1):
                bin_error = hist.GetBinError(bin_index)
                new_hist.SetBinError(bin_index, bin_error)
        hist = new_hist

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPalette(57)
    hist.SetTitle(args.infilename)
    hist.SetLabelSize(0.05, 'xy')
    hist.SetLabelSize(0.03, 'z')
    hist.SetTitleSize(0.05, 'xyz')
    if args.smallxnumbers:
        hist.SetLabelSize(0.035, 'x')
    hist.SetTitleOffset(1.9, 'y')
    hist.SetLabelOffset(0.010, 'y')
    canvas = ROOT.TCanvas('c1', 'c1', 800, 800)
    margin = 0.18
    canvas.SetRightMargin(margin/2)
    canvas.SetLeftMargin(margin)
    canvas.SetTopMargin(margin)
    canvas.SetBottomMargin(margin)


    hist.Draw('colz' if args.colz else '')
    text = ROOT.TPaveText(0.58, 0.73, 0.85, 0.82, 'NB NDC')
    #text = ROOT.TPaveText(0.63, 0.77, 0.87, 0.89, 'NB NDC')
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
