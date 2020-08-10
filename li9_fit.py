import argparse
from array import array
import math
import os
import sqlite3

def main(infiles, outfile, database, label, ntag, site, update_db):
    import ROOT
    ROOT.gROOT.SetBatch(True)

    ntag_str = '_ntag' if ntag else ''
    name = f'{label}energy{ntag_str}'

    with sqlite3.Connection(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('SELECT Rate_Hz FROM li9_muon_rates '
                'WHERE Hall = ? AND EnergyClass = ? AND NeutronTag = ?',
                (site, label, int(ntag)))
        row = cursor.fetchone()
        rmu = row['Rate_Hz']
        cursor.execute('SELECT Peak, Resolution, DetNo FROM delayed_energy_fits '
                'WHERE Hall = ?', (site,))
        bounds = {}
        rows = cursor.fetchall()
        for row in rows:
            peak = row['Peak']
            width = row['Resolution']
            upper = peak + 3 * width
            lower = peak - 3 * width
            bounds[row['DetNo']] = (lower, upper)
    energy_cut_strs = [f'(detector[0] == {det} && energy[1] > {lower} '
            f' && energy[1] < {upper})'
            for det, (lower, upper) in bounds.items()]
    energy_cut_str = '(' + '||'.join(energy_cut_strs) + ')'

    print('creating TCanvas')
    canvas = ROOT.TCanvas('c1', 'c1', 800, 800)
    margin = 0.18
    canvas.SetRightMargin(margin)
    canvas.SetLeftMargin(margin)
    canvas.SetTopMargin(margin)
    canvas.SetBottomMargin(margin)


    # hardcoded bin boundaries
    bins = array('f', [2.00e-03, 2.00e-02, 4.00e-02, 6.00e-02, 8.00e-02, 1.00e-01,
    2.00e-01, 3.00e-01, 4.00e-01, 5.00e-01, 7.50e-01, 1.00e+00,
    1.25e+00, 1.50e+00, 1.75e+00, 2.00e+00, 2.50e+00, 3.00e+00,
    4.00e+00, 5.00e+00, 6.00e+00, 7.00e+00, 8.00e+00, 9.00e+00, 1.00e+01])

    test_file = ROOT.TFile(infiles[0], 'READ')
    test_histogram = test_file.Get(name)
    if not test_histogram:
        test_file.Close()
        data = ROOT.TChain('ad_events')
        print(f'Processing {len(infiles)} files')
        for filename in sorted(infiles):
            data.Add(filename)
        data.SetBranchStatus('*', 0)

        t_last_muon_hist = ROOT.TH1F(f'{label}energy{ntag_str}',
                f'{label}energy{ntag_str}', 24, 0.002, 10)
        t_last_muon_hist.GetXaxis().Set(24, bins)

        data.SetBranchStatus(f'dt_{label}energy_muon{ntag_str}', 1)
        data.SetBranchStatus('energy', 1)
        data.SetBranchStatus('dr_to_prompt', 1)
        data.SetBranchStatus('dt_to_prompt', 1)
        data.SetBranchStatus('detector', 1)
        draw_expr = f'dt_{label}energy_muon{ntag_str}/1e9 >> {t_last_muon_hist.GetName()}'
        selection_expr = (f'dt_{label}energy_muon{ntag_str} < 1e10'
                f' && dt_{label}energy_muon{ntag_str} > 2e6 && energy[0] > 3.5'
                f' && ({energy_cut_str})'
                ' && dr_to_prompt[1] + 1000/600e3 * dt_to_prompt[1] < 800')

        print('drawing histogram')
        data.Draw(draw_expr, selection_expr)
        density(t_last_muon_hist)
    else:
        t_last_muon_hist = test_histogram

    print('performing fit')
    fitter = ROOT.TF1("li9fitter", li9fit, 0.002, 10, 4)
    n_IBD_guess = t_last_muon_hist.GetBinContent(12)/rmu  # 1s bin
    n_bkg_guess = abs(
            t_last_muon_hist.GetBinContent(6)  # 0.1s bin
            - t_last_muon_hist.GetBinContent(12)
    ) / (rmu + 4)
    n_bb_guess = abs(
            t_last_muon_hist.GetBinContent(1)
            - t_last_muon_hist.GetBinContent(2)
    ) / (rmu + 50)
    fitter.SetParameters(n_bkg_guess, n_IBD_guess, n_bb_guess)
    fitter.FixParameter(3, rmu)
    fitter.SetParLimits(0, 0, 1e9)
    fitter.SetParLimits(1, 0, 1e9)
    fitter.SetParLimits(2, 0, 1e9)
    results = t_last_muon_hist.Fit(fitter, "QS", "", 0.002, 10)
    for name, value, error in zip(['N_li9', 'N_IBD', 'N_BB', 'R_mu'],
            results.GetParams(), results.GetErrors()):
        print(f'{name}: {value} +/- {error}')
    print('Chi2:', results.Chi2())
    print('NDF:', results.Ndf())

    if update_db:
        print('Updating database')
        params = [x[0] for x in zip(results.GetParams(), range(4))]
        errors = [x[0] for x in zip(results.GetErrors(), range(4))]
        N_Li9_He8, N_IBD, N_BB, _ = params
        N_Li9_He8_error, N_IBD_error, N_BB_error, _ = errors
        chi2 = results.Chi2()
        num_bins = fitter.GetNumberFitPoints()
        num_params = fitter.GetNpar()
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''
            INSERT OR REPLACE INTO li9_fits
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
            (site, label, int(ntag), N_Li9_He8, N_Li9_He8_error, N_IBD, N_IBD_error,
                N_BB, N_BB_error, chi2, num_bins, num_params))
            conn.commit()

    only_ibds = ROOT.TF1("only_ibds", only_ibd_term, 0.002, 10, 2)
    only_ibds.SetParameters(results.Parameter(1), results.Parameter(3))
    only_ibds.SetLineColor(3)
    only_ibds.Draw("same")
    ibds_li9 = ROOT.TF1("ibds_li9", ibd_and_li9, 0.002, 10, 3)
    ibds_li9.SetParameters(results.Parameter(0), results.Parameter(1),
            results.Parameter(3))
    ibds_li9.SetLineColor(4)
    ibds_li9.Draw("same")
    ibds_li9_he8 = ROOT.TF1("ibds_li9_he8", ibd_li9_he8, 0.002, 10, 3)
    ibds_li9_he8.SetParameters(results.Parameter(0), results.Parameter(1),
            results.Parameter(3))
    ibds_li9_he8.SetLineColor(6)
    ibds_li9_he8.Draw("same")
    legend = ROOT.TLegend(0.58, 0.63, 0.85, 0.82, '', 'NB NDC')
    legend.AddEntry(t_last_muon_hist, "Data")
    legend.AddEntry(only_ibds, "IBD")
    legend.AddEntry(ibds_li9, "IBD + {}^{9}Li")
    legend.AddEntry(ibds_li9_he8, "IBD + {}^{9}Li + {}^{8}He")
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

    outroot = os.path.splitext(outfile)[0] + '.root'
    writeout = ROOT.TFile(outroot, 'RECREATE')
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
    parser.add_argument('database')
    parser.add_argument('infiles', nargs='+')
    parser.add_argument('--site', type=int)
    parser.add_argument('--energy', choices=['low', 'mid', 'high'], required=True)
    parser.add_argument('--ntag', action='store_true')
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    infiles = args.infiles
    outfile = args.outfile
    main(infiles, outfile, args.database, args.energy, args.ntag, args.site,
            args.update_db)
