import argparse

import numpy as np

import prediction


def chi_square(constants, fit_params):
    """Compute the chi-square value for a single set of parameters.

    Starting with just comparing EH1 to EH3, 1 set of reco bins.
    """
    observed = np.zeros_like(constants.observed_spectra[(3, 1)])
    for halldet in prediction.far_ads:
        observed += constants.observed_spectra[halldet]

    predicted = prediction.predict_halls(constants, fit_params)[1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    args = parser.parse_args()
    constants = prediction.default_constants(args.database)
    starting_params = prediction.FitParams(0.15, 2.5e-3)
    print(chi_square(constants, starting_params))
