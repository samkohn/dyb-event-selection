import argparse
from array import array
import math

def main(infiles, outfile):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    data = ROOT.TChain('ad_events')
    print(f'Processing {len(infiles)} files')
    for filename in sorted(infiles):
        data.Add(filename)

    # hardcoded bin boundaries
    bins = array('f', [2.00e-03, 2.00e-02, 4.00e-02, 6.00e-02, 8.00e-02, 1.00e-01,
    2.00e-01, 3.00e-01, 4.00e-01, 5.00e-01, 7.50e-01, 1.00e+00,
    1.25e+00, 1.50e+00, 1.75e+00, 2.00e+00, 2.50e+00, 3.00e+00,
    4.00e+00, 5.00e+00, 6.00e+00, 7.00e+00, 8.00e+00, 9.00e+00, 1.00e+01])


    label = 'high'
    t_last_muon_hist = ROOT.TH1F(f'{label}energy', f'{label}energy', 24, 0.002, 10)
    t_last_muon_hist.GetXaxis().Set(24, bins)

    draw_expr = f'dt_{label}energy_muon/1e9 >> {t_last_muon_hist.GetName()}'
    selection_expr = (f'multiplicity == 2 && dt_{label}energy_muon < 1e10'
            f' && dt_{label}energy_muon > 2e6 && energy[0] > 3.5')

    print('creating TCanvas')
    canvas = ROOT.TCanvas('c1', 'c1', 800, 800)
    margin = 0.18
    canvas.SetRightMargin(margin)
    canvas.SetLeftMargin(margin)
    canvas.SetTopMargin(margin)
    canvas.SetBottomMargin(margin)

    print('drawing histogram')
    data.Draw(draw_expr, selection_expr)
    density(t_last_muon_hist)

    print('performing fit')
    fitter = ROOT.TF1("li9fitter", li9fit, 0.002, 10, 4)
    fitter.SetParameters(50, 100, 50, 0.01)
    results = t_last_muon_hist.Fit(fitter, "QS", "", 0.002, 10)
    for name, value, error in zip(['N_li9', 'N_IBD', 'N_BB', 'R_mu'],
            results.GetParams(), results.GetErrors()):
        print(f'{name}: {value} +/- {error}')
    print('Chi2:', results.Chi2())
    print('NDF:', results.Ndf())

    only_ibds = ROOT.TF1("only_ibds", only_ibd_term, 0.002, 10, 2)
    only_ibds.SetParameters(results.Parameter(1), results.Parameter(3))
    only_ibds.SetLineColor(3)
    only_ibds.Draw("same")
    ibds_li9 = ROOT.TF1("ibds_li9", ibd_and_li9, 0.002, 10, 3)
    ibds_li9.SetParameters(results.Parameter(0), results.Parameter(1),
            results.Parameter(3))
    ibds_li9.SetLineColor(4)
    ibds_li9.Draw("same")
    legend = ROOT.TLegend(0.58, 0.63, 0.85, 0.82, '', 'NB NDC')
    legend.AddEntry(t_last_muon_hist, "Data")
    legend.AddEntry(only_ibds, "IBD")
    legend.AddEntry(ibds_li9, "IBD + {}^{9}Li")
    legend.AddEntry(fitter, "IBD + {}^{9}Li + {}^{8}He + {}^{12}B")
    legend.Draw("same")

    ROOT.gStyle.SetOptStat(0)
    t_last_muon_hist.SetLabelSize(0.03, 'xy')
    t_last_muon_hist.SetTitleSize(0.04, 'xy')
    t_last_muon_hist.SetLineWidth(3)
    t_last_muon_hist.GetXaxis().SetTitle('Time since last muon [s]')
    t_last_muon_hist.GetYaxis().SetTitle('IBD candidates')
    ROOT.gPad.SetLogx()



    print('updating TCanvas')
    ROOT.gPad.Update()
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    ROOT.gPad.Modified()

    ROOT.gPad.Print(outfile)
    canvas.Close()

    writeout = ROOT.TFile('test_li9_fit.root', 'RECREATE')
    t_last_muon_hist.SetDirectory(writeout)
    t_last_muon_hist.Write()
    writeout.Write()
    writeout.Close()

def li9fit(x, par):
    """x[0] = time since last muon; par = [nli9he8, nibd, nbb, rmu]"""
    n_li9he8, n_ibd, n_bb, r_mu = par
    t_mu = x[0]
    lam_li9 = r_mu + 1/0.25723
    lam_he8 = r_mu + 1/0.17160
    lam_bb = r_mu + 2/0.0202
    li9_ratio = 0.85
    bkg_term = n_li9he8 * (li9_ratio * lam_li9 * math.exp(-lam_li9 * t_mu)
            + (1 - li9_ratio) * lam_he8 * math.exp(-lam_he8 * t_mu))
    ibd_term = n_ibd * r_mu * math.exp(-r_mu * t_mu)
    bb_term = n_bb * lam_bb * math.exp(-lam_bb * t_mu)
    return bkg_term + ibd_term + bb_term


def only_ibd_term(x, par):
    """x[0] = time since last muon; par = [nibd, rmu]"""
    n_ibd, r_mu = par
    t_mu = x[0]
    return n_ibd * r_mu * math.exp(-r_mu * t_mu)

def ibd_and_li9(x, par):
    """x[0] = time since last muon; par = [nli9he8, nibd, rmu]"""
    n_li9he8, n_ibd, r_mu = par
    t_mu = x[0]
    lam_li9 = r_mu + 1/0.25723
    li9_ratio = 0.85
    li9_term = n_li9he8 * li9_ratio * lam_li9 * math.exp(-lam_li9 * t_mu)
    ibd_term = only_ibd_term(x, [n_ibd, r_mu])
    return li9_term + ibd_term

def ibd_li9_he8(x, par):
    """x[0] = time since last muon; par = [nli9he8, nibd, rmu]"""
    n_li9he8, n_ibd, r_mu = par
    t_mu = x[0]
    lam_li9 = r_mu + 1/0.25723
    lam_he8 = r_mu + 1/0.17160
    li9_ratio = 0.85
    bkg_term = n_li9he8 * (li9_ratio * lam_li9 * math.exp(-lam_li9 * t_mu)
            + (1 - li9_ratio) * lam_he8 * math.exp(-lam_he8 * t_mu))
    ibd_term = n_ibd * r_mu * math.exp(-r_mu * t_mu)
    return bkg_term + ibd_term



def density(h):
    """Divide each bin content and error by its width bin width."""
    h.Sumw2()
    for bin_index in range(1, h.GetNbinsX() + 1):
        val = h.GetBinContent(bin_index)
        width = h.GetBinWidth(bin_index)
        error = h.GetBinError(bin_index)
        new_val = val / width
        new_error = error / width
        h.SetBinContent(bin_index, new_val)
        h.SetBinError(bin_index, new_error)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile')
    parser.add_argument('infiles', nargs='+')
    args = parser.parse_args()
    infiles = args.infiles
    outfile = args.outfile
    main(infiles, outfile)
