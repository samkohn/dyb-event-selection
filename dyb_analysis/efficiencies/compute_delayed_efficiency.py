import argparse

import common
import fit_delayed

def main(infilename, outfilename, database):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    data = infile.Get('toy')
    delayed_energy = 'res_d'
    selection_all = 'target >= 0 && Volume_nCap >= 0 && res_p > 1.5'
    delayed_spectrum = ROOT.TH1F('spectrum', 'spectrum', 525, 1.5, 6.75)
    delayed_spectrum.Sumw2()
    data.Draw(delayed_energy + '>>spectrum', selection_all, 'goff')
    mu0 = 2.3
    sigma0 = 0.14
    alpha0 = 0.8
    scale0 = 2
    approx_integral = delayed_spectrum.Integral(1,
            delayed_spectrum.GetNbinsX())
    bin_width = delayed_spectrum.GetBinWidth(1)
    norm0 = approx_integral * bin_width
    fitter = ROOT.TF1("calo_fitter", fit_delayed.calorimeter_fn, 1.8, 3.1, 5)
    fitter.SetParameters(mu0, sigma0, scale0, alpha0, norm0)
    if outfilename is None:
        options = 'QN0S'
    else:
        options = 'QS'
    fit_result = delayed_spectrum.Fit(fitter, options)
    mu, sigma, scale, alpha, norm = [fit_result.Parameter(i) for i in range(5)]
    mu_err, sigma_err, scale_err, alpha_err, norm_err = [fit_result.ParError(i)
            for i in range(5)]
    delayed_spectrum.GetXaxis().SetRangeUser(1.5, 3.4)
    if outfilename is not None:
        ROOT.gPad.Print(outfilename)
    if database is None:
        print('mu, sigma, scale, alpha, norm')
        print(mu, sigma, scale, alpha, norm)
    else:
        # Special values for site and ad for MC results
        site = 0
        ad = 0
        with common.get_db(database) as conn:
            c = conn.cursor()
            c.execute('''INSERT OR REPLACE INTO delayed_energy_fits
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
            (site, ad, mu, mu_err, sigma, sigma_err, scale, scale_err, alpha,
                alpha_err, norm, norm_err))
            conn.commit()
    lower_energy = mu - 3*sigma
    upper_energy = mu + 3*sigma
    # Now compute the efficiencies
    selection_cut = '{} && res_d >= {} && res_d <= {}'.format(selection_all,
            lower_energy, upper_energy)
    selection_all_in_volume = 'target >= 0 && Volume_nCap == {} && res_p > 1.5'
    selection_cut_in_volume = '{} && res_d >= {} && res_d <= {}'.format(
            selection_all_in_volume, lower_energy, upper_energy)
    def efficiency(selection, cut_selection):
        total = data.Draw(delayed_energy, selection, 'goff')
        reduced = data.Draw(delayed_energy, cut_selection, 'goff')
        return reduced/total
    total_efficiency = efficiency(selection_all, selection_cut)
    efficiencies = {'total': total_efficiency}
    volumes = {
            'GdLS': 1,
            'IAV': 2,
            'LS': 3,
            'OAV': 4,
            }
    for name, volume in volumes.items():
        selection = selection_all_in_volume.format(volume)
        cut_selection = selection_cut_in_volume.format(volume)
        efficiencies[name] = efficiency(selection, cut_selection)
    if database is None:
        print(efficiencies)
    else:
        with common.get_db(database) as conn:
            c = conn.cursor()
            for name, eff_value in efficiencies.items():
                c.execute('''INSERT OR REPLACE INTO delayed_energy_efficiency
                VALUES (?, ?)''', (name, eff_value))

    infile.Close()
    return




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('--update-db', default=None)
    args = parser.parse_args()
    main(args.input, args.output, args.update_db)
