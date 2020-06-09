import argparse
import math
import sqlite3
import time

def crystal_ball(x, p):
    """Crystal Ball function.

    x (input data) is 1-dimensional, i.e. len(x) == 1

    p (parameters) has 5 entries:
      p[0] = mean of gaussian
      p[1] = variance of gaussian
      p[2] = cutoff between peak and tail
      p[3] = power of the power-law tail
      p[4] = normalization
    """
    inval = x[0]
    mu, sigma, alpha, n, norm = p
    if alpha < 0:
        alpha = -alpha
    scaled = (inval - mu)/sigma
    if scaled > -alpha:
        return norm * math.exp(-scaled*scaled/2)
    else:
        return norm * (n/alpha)**n * math.exp(-alpha*alpha/2) * (
                n/alpha - alpha - scaled)**(-n)

def calorimeter_fn(x, p):
    """Calorimeter function.

    x (input data) is 1-dimensional, i.e. len(x) == 1, x == [value].

    p (parameters) has 5 entries:
      p[0] = true peak energy (~mu)
      p[1] = energy resolution (~sigma)
      p[2] = scale of exponential (~lambda, i.e. lam*exp(lam*E))
      p[3] = fraction of events in true peak (~alpha)
      p[4] = normalization
    """
    try:
        inval = x[0]
        peak, resolution, exp_scale, peak_fraction, norm = p
        erf = math.erf
        exp = math.exp
        sqrt2 = math.sqrt(2)
        sqrt2pi = math.sqrt(2*math.pi)
        peak_norm = peak_fraction / (resolution * sqrt2pi)
        peak_exp = exp(-(inval - peak)**2/(2 * resolution**2))
        peak_term = peak_norm * peak_exp

        tail_coeff = (1 - peak_fraction) * exp_scale / (exp(exp_scale * peak) - 1)
        tail_exp = exp(((resolution * exp_scale)**2 + 2 * exp_scale * inval)/2)
        tail_erf1 = erf((peak - inval - resolution**2 * exp_scale) / (sqrt2 *
            resolution))
        tail_erf2 = erf((-inval - resolution**2 * exp_scale) / (sqrt2 *
            resolution))
        tail_term = tail_coeff * tail_exp * (tail_erf1 - tail_erf2)
        return norm * (peak_term + tail_term)
    except:
        print([y for y in x])
        print([y for y in p])
        raise


def main(infilename, outfilename, site, ad, database):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    spectrum_2d = infile.Get('final')
    delayed_spectrum = spectrum_2d.ProjectionY()
    mu0 = 2.3
    sigma0 = .14
    alpha0 = 0.8
    scale0 = 2
    # Approximate the normalization as the area under the lower half of the
    # x-range of the histogram (1.5-6.75 MeV)
    approx_integral = delayed_spectrum.Integral(1,
            delayed_spectrum.GetNbinsX()//2)
    bin_width = delayed_spectrum.GetBinWidth(1)
    norm0 = approx_integral * bin_width

    fitter = ROOT.TF1("calo_fitter", calorimeter_fn, 1.5, 12, 5)
    fitter.SetParameters(mu0, sigma0, scale0, alpha0, norm0)

    if outfilename is None:
        options = 'QN0S'
    else:
        options = 'QS'
    try:
        fit_result = delayed_spectrum.Fit(fitter, options, '', 1.6, 2.8)
    except:
        delayed_spectrum.Draw()
        fitter.Draw()
        time.sleep(10)
        ROOT.gPad.Print('error.pdf')
        raise
    mu, sigma, scale, alpha, norm = [fit_result.Parameter(i) for i in range(5)]
    mu_err, sigma_err, scale_err, alpha_err, norm_err = [fit_result.ParError(i)
            for i in range(5)]
    chi_square = fit_result.Chi2()
    num_bins = fitter.GetNumberFitPoints()
    num_params = fitter.GetNpar()
    delayed_spectrum.GetXaxis().SetRangeUser(1.5, 3.4)
    if outfilename is not None:
        ROOT.gPad.Print(outfilename)
    infile.Close()
    if database is not None:
        with sqlite3.Connection(database) as conn:
            c = conn.cursor()
            c.execute('''INSERT OR REPLACE INTO delayed_energy_fits
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
            (site, ad, mu, mu_err, sigma, sigma_err, scale, scale_err, alpha,
                alpha_err, norm, norm_err, chi_square, num_bins, num_params))
            conn.commit()
    else:
        print(mu, sigma, scale, alpha, norm)
        print(mu_err, sigma_err, scale_err, alpha_err, norm_err)
        print(f'Chi2/NDF: {chi_square:.03f} / ({num_bins} - {num_params})')
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('-o', '--output', default=None)
    parser.add_argument('--site', type=int)
    parser.add_argument('--ad', type=int)
    parser.add_argument('--update-db', default=None)
    args = parser.parse_args()
    main(args.input, args.output, args.site, args.ad, args.update_db)
