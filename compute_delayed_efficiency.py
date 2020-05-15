import argparse
import sqlite3

import fit_delayed

def main(infilename, outfilename, database):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    data = infile.Get('toy')
    delayed_energy = 'res_d'
    selection_all = 'target >= 0 && Volume_nCap >= 0'
    delayed_spectrum = ROOT.TH1F('spectrum', 'spectrum', 400, 1, 5)
    data.Draw(delayed_energy + '>>spectrum', selection_all, 'goff')
    mu0 = 2.3
    sigma0 = 0.14
    alpha0 = 1.5
    n0 = 2
    norm0 = delayed_spectrum.GetMaximum()
    fitter = ROOT.TF1("cb_fitter", fit_delayed.crystal_ball, 1.8, 3.1, 5)
    fitter.SetParameters(mu0, sigma0, alpha0, n0, norm0)
    if outfilename is None:
        options = 'QN0S'
    else:
        options = 'QS'
    fit_result = delayed_spectrum.Fit(fitter, options)
    mu, sigma, alpha, n, norm = [fit_result.Parameter(i) for i in range(5)]
    delayed_spectrum.GetXaxis().SetRangeUser(1.5, 3.4)
    if outfilename is not None:
        ROOT.gPad.Print(outfilename)
    if database is None:
        print('mu, sigma, alpha, n, norm')
        print(mu, sigma, alpha, n, norm)
    else:
        # Special values for site and ad for MC results
        site = 0
        ad = 0
        with sqlite3.Connection(database) as conn:
            c = conn.cursor()
            c.execute('''INSERT OR REPLACE INTO delayed_energy_fits
            VALUES (?, ?, ?, ?, ?, ?, ?)''',
            (site, ad, mu, sigma, alpha, n, norm))
            conn.commit()
    lower_energy = mu - 3*sigma
    upper_energy = mu + 3*sigma
    # Now compute the efficiencies
    selection_cut = '{} && res_d >= {} && res_d <= {}'.format(selection_all,
            lower_energy, upper_energy)
    selection_all_in_volume = 'target >= 0 && Volume_nCap == {}'
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
        with sqlite3.Connection(database) as conn:
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
