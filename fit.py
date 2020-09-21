import argparse
from pprint import pprint

import numpy as np

import prediction as pred


def chi_square(constants, fit_params):
    """Compute the chi-square value for a single set of parameters.
    """
    chi_square = 0
    observed = {ad: constants.observed_candidates[ad] for ad in pred.far_ads}
    predicted, reco_bins = pred.predict_ad_to_ad_obs(constants, fit_params)
    pprint(observed)
    pprint(predicted)
    #Main part
    for (far_halldet, near_halldet), n_predicted in predicted.items():
        n_observed = observed[far_halldet]
        sigma_observed = np.sqrt(n_observed)  # TODO placeholder error
        numerator = np.power(n_observed - n_predicted, 2)
        denominator = np.power(sigma_observed, 2)
        chi_square += np.sum(numerator/denominator)  # sum over energy bins

    # Pull terms
    for halldet, pull in fit_params.pull_bg.items():
        numerator = pull*pull
        denominator = 245 * 245  # TODO placeholder error
        chi_square += numerator/denominator
    return chi_square




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    args = parser.parse_args()
    constants = pred.default_constants(args.database)
    starting_params = pred.FitParams(0.15, 2.5e-3, pred.ad_dict(0))
    print(chi_square(constants, starting_params))
