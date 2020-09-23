import argparse
from pprint import pprint

import numpy as np
from scipy.optimize import minimize

import prediction as pred


def chi_square(constants, fit_params):
    """Compute the chi-square value for a single set of parameters.
    """
    chi_square = 0
    observed = {ad: constants.observed_candidates[ad] for ad in pred.far_ads}
    predicted, reco_bins = pred.predict_ad_to_ad_obs(constants, fit_params)
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
        denominator = 0.010*0.010  # TODO relative error on accidentals
        chi_square += numerator/denominator
    for halldet, pull in fit_params.pull_near_stat.items():
        numerator = pull * pull
        denominator = np.power(1/np.sqrt(constants.observed_candidates[halldet]), 2)
        chi_square += numerator/denominator
    return chi_square

def objective_fn(x, constants):
    """Convert arguments from scipy to desired format and return chi_square.
    """
    fit_params = pred.FitParams.from_list(x)
    return chi_square(constants, fit_params)

def fit(starting_params, constants):
    """Perform the fit with the given starting parameters.
    """
    positive = (0, None)
    nobound = (None, None)
    bounds = ([positive, (2.5e-3, 2.5e-3)]
        + [nobound] * 8
        + [nobound] * 4
    )

    result = minimize(objective_fn, np.array(starting_params.to_list()),
            args=(constants,), bounds=bounds)
    return result



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    args = parser.parse_args()
    constants = pred.default_constants(args.database)
    starting_params = pred.FitParams(
            0.15,
            2.5e-3,
            pred.ad_dict(0),
            pred.ad_dict(0, halls='near')
    )
    print(chi_square(constants, starting_params))
    result = fit(starting_params, constants)
    print(result)
    cov = np.array(result.hess_inv.todense())
    print(np.sqrt(cov[0,0]))
    print(np.sqrt(cov[1, 1]))
