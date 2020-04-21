import argparse
import math

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


def main(infilename, outfilename):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    spectrum_2d = infile.Get('final')
    delayed_spectrum = spectrum_2d.ProjectionY()
    delayed_spectrum.Sumw2()
    mu0 = 2.3
    sigma0 = .14
    alpha0 = 1.5
    n0 = 2
    norm0 = delayed_spectrum.GetMaximum()

    fitter = ROOT.TF1("cb_fitter", crystal_ball, 1.8, 3.1, 5)
    fitter.SetParameters(mu0, sigma0, alpha0, n0, norm0)

    delayed_spectrum.Fit(fitter)
    delayed_spectrum.GetXaxis().SetRangeUser(1.5, 3.4)
    ROOT.gPad.Print(outfilename)
    infile.Close()
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    args = parser.parse_args()
    main(args.input, args.output)
