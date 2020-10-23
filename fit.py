import argparse
from pprint import pprint
import sqlite3
import sys

import numpy as np
from scipy.optimize import least_squares

import prediction as pred

def chi_square(constants, fit_params, return_array=False, debug=False, near_ads=None,
        rate_only=False):
    """Compute the chi-square value for a single set of parameters.

    Set return_array=True to return an array of terms rather than the sum.
    """
    chi_square = 0
    if near_ads is None:
        near_ads = pred.near_ads
    far_ads = pred.far_ads
    if rate_only:
        observed = {ad: np.sum(constants.observed_candidates[ad], keepdims=True) for ad in far_ads}
    else:
        observed = {ad: constants.observed_candidates[ad] for ad in far_ads}
    predicted, reco_bins = pred.predict_ad_to_ad_obs(constants, fit_params)
    if rate_only:
        for key, val in predicted.items():
            predicted[key] = np.sum(val, keepdims=True)
        num_bins = 1
        num_pulls = 26
    else:
        num_bins = 37
        num_pulls = 26
    if debug:
        pprint(observed)
        pprint(predicted)
    return_array_values = np.zeros((len(far_ads) * len(near_ads) * num_bins + num_pulls,))
    term_index = 0

    #Main part
    for (far_halldet, near_halldet), n_predicted in predicted.items():
        if far_halldet not in far_ads or near_halldet not in near_ads:
            continue
        n_observed = observed[far_halldet]
        sigma_observed = np.sqrt(n_observed)  # TODO placeholder error
        numerator = np.power(n_observed - n_predicted, 2)
        denominator = np.power(sigma_observed, 2)
        ratio = numerator/denominator
        chi_square += np.sum(numerator/denominator)  # sum over energy bins
        return_array_values[term_index : term_index + len(ratio)] = ratio
        term_index += len(ratio)

    # Pull terms
    for halldet, pull in fit_params.pull_bg.items():
        numerator = pull*pull
        denominator = 0.010*0.010  # TODO relative error on accidentals
        chi_square += numerator/denominator
        return_array_values[term_index] = numerator/denominator
        term_index += 1
    for halldet, pull in fit_params.pull_near_stat.items():
        numerator = pull * pull
        denominator = np.power(np.sum(constants.observed_candidates[halldet],
            keepdims=True), -1)
        chi_square += numerator/denominator
        return_array_values[term_index] = numerator/denominator
        term_index += 1
    for core in range(1, 7):
        pull = fit_params.pull_reactor[core]
        numerator = pull * pull
        denominator = 0.008 * 0.008
        chi_square += numerator/denominator
        return_array_values[term_index] = numerator/denominator
        term_index += 1
    for halldet, pull in fit_params.pull_efficiency.items():
        numerator = pull * pull
        denominator = 0.003 * 0.003
        chi_square += numerator/denominator
        return_array_values[term_index] = numerator/denominator
        term_index += 1
    if return_array:
        return return_array_values
    else:
        return sum(chi_square)

def residual_fn(x, constants, near_ads=None, rate_only=False):
    """Convert arguments from scipy to desired format and return the residuals.
    """
    fit_params = pred.FitParams.from_list(x)
    # Take the square root because the fitter wants the linear residuals
    # and squares them on its own...
    residuals = np.sqrt(chi_square(constants, fit_params, return_array=True,
        near_ads=near_ads, rate_only=rate_only))
    return residuals

def residual_frozen_param(frozen_dict, near_ads, rate_only):
    """Return a residuals function but with certain parameters frozen.

    Frozen parameters should be expressed as ``index: value`` pairs,
    where the index is relative to FitParams.to_list() ordering.
    """
    num_frozen = len(frozen_dict)
    def residual(x, constants):
        f"""Convert arguments from scipy to desired format and return the residuals.

        The first parameter should be an array/list of the non-frozen params.
        The second parameter should be a FitConstants object.
        The frozen params will be inserted automatically.

        Frozen params:

        {frozen_dict}
        """
        param_list = []
        x_index = 0
        for i in range(len(x) + num_frozen):
            if i in frozen_dict:
                param_list.append(frozen_dict[i])
            else:
                param_list.append(x[x_index])
                x_index += 1
        return residual_fn(param_list, constants, near_ads=near_ads, rate_only=rate_only)
    return residual

def fit_lsq_frozen(starting_params, constants, frozen_params, near_ads, rate_only):
    """Perform the fit with certain parameters frozen."""
    frozen_params_dict = {}
    all_params = starting_params.to_list()
    x0 = []
    for i, param in enumerate(all_params):
        if i in frozen_params:
            frozen_params_dict[i] = param
        else:
            x0.append(param)
    residual = residual_frozen_param(frozen_params_dict, near_ads, rate_only)
    result = least_squares(residual, x0, args=(constants,), method='trf')
    return result

def sigma_searcher(fit_params, constants, side='both'):
    """Find the +/- 1-sigma values given the best fit parameters.

    ``side`` could be 'upper', 'lower', or 'both'

    Returns the value(s) of theta13 at the error boundary/ies in a list
    of length 1 or 2.

    If side is 'both', the 2 entries are always returned [upper, lower]
    """
    best_theta13 = fit_params.theta13
    best_sin2 = np.power(np.sin(2 * best_theta13), 2)
    guess_1sigma_sin2 = 0.0045
    guess_upper_sin2 = best_sin2 + guess_1sigma_sin2
    guess_lower_sin2 = best_sin2 - guess_1sigma_sin2
    guess_upper = 0.5 * np.arcsin(np.sqrt(guess_upper_sin2))
    guess_lower = 0.5 * np.arcsin(np.sqrt(guess_lower_sin2))
    near_ads = None
    min_chi2 = chi_square(constants, fit_params, near_ads=near_ads)
    to_return = []
    if side in ('upper', 'both'):
        result = least_squares(sigma_search_resid, guess_upper,
                args=(constants, min_chi2), method='trf', xtol=1e-3)
        upper_limit = result.x[0]
        to_return.append(upper_limit)
    if side in ('lower', 'both'):
        result = least_squares(sigma_search_resid, guess_lower,
                args=(constants, min_chi2), method='trf', xtol=1e-3)
        lower_limit = result.x[0]
        to_return.append(lower_limit)
    return to_return

def sigma_search_resid(x, constants, min_chi2):
    """Return the difference between the delta chi-square and 1.

    x is just [theta13_search]
    """
    print('Running fit')
    starting_params = pred.FitParams(
            x[0],
            pred.ad_dict(0),
            pred.ad_dict(0, halls='near'),
            pred.core_dict(0),
            pred.ad_dict(0),
    )
    fit_result = fit_lsq_frozen(starting_params, constants, range(9), None)
    chi2 = np.power(fit_result.fun, 2).sum()
    print(f'Trial chi2 = {chi2:.5f}')
    return chi2 - min_chi2 - 1


def chi_square_grid(starting_params, constants, theta13_values):
    """Return a grid of optimal chi-square values for the given parameter grid.
    """
    best_theta13 = starting_params.theta13
    result = np.zeros((len(theta13_values)))
    for i, theta13 in enumerate(theta13_values):
        starting_params.theta13 = theta13
        fit_result = fit_lsq_frozen(starting_params, constants, range(9), None)
        result[i] = np.power(fit_result.fun, 2).sum()
    starting_params.theta13 = best_theta13  # restore original value
    return result

def save_result(database, source, index, theta13_best, theta13_low_err, theta13_up_err, chi_square):
    """Save the specified results to the database.
    """
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        sin2_best, sin2_low, sin2_up = np.power(np.sin(2*np.array(
            [theta13_best, theta13_low_err, theta13_up_err]
        )), 2)
        cursor.execute('''INSERT OR REPLACE INTO toymc_tests
        VALUES (?, ?, ?, ?, ?, ?, ?)''',
        (source, index, sin2_best, chi_square, 3, sin2_up, sin2_low))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("--no-fit", action='store_true')
    parser.add_argument("--scan", action='store_true')
    parser.add_argument("--source")
    parser.add_argument("--source-index", type=int)
    parser.add_argument("--update-db")
    parser.add_argument("--shape", action='store_true')
    parser.add_argument("--dm2ee", type=float)
    parser.add_argument("--debug", action='store_true')
    args = parser.parse_args()
    rate_only = not args.shape
    constants = pred.load_constants(args.config)
    starting_params = pred.FitParams(
            0.15,
            2.48e-3,
            pred.ad_dict(0),
            pred.ad_dict(0, halls='near'),
            pred.core_dict(0),
            pred.ad_dict(0),
    )
    frozen_params = list(range(2, 10))
    if args.dm2ee is not None:
        starting_params.m2_ee = args.dm2ee
        frozen_params = [1] + frozen_params  # Don't modify m2_ee in fit
    if args.dm2ee is None and rate_only:
        raise ValueError("Can't fit dm2_ee in rate-only")
    print(starting_params)
    print(chi_square(constants, starting_params, rate_only=rate_only, debug=args.debug))
    if args.debug:
        print(chi_square(constants, starting_params, return_array=True,
            rate_only=rate_only))
    if not args.no_fit:
        near_ads = None
        result = fit_lsq_frozen(starting_params, constants, frozen_params,
                near_ads=near_ads, rate_only=rate_only)
        print(repr(result.x))
        print('sin22theta13 =', np.power(np.sin(2*result.x[0]), 2))
        print(result.success)
        print(result.message)
        if not result.success:
            sys.exit(0)
        if 1 in frozen_params:
            fit_params = pred.FitParams.from_list(
                [result.x[0], starting_params.m2_ee] + [0] * 8 + result.x[1:].tolist()
            )
        else:
            fit_params = pred.FitParams.from_list(
                [result.x[0], result.x[1]] + [0] * 8 + result.x[2:].tolist()
            )
        print(chi_square(constants, fit_params, return_array=False, near_ads=near_ads,
            rate_only=rate_only))
        if args.debug:
            print(chi_square(constants, fit_params, return_array=True, near_ads=near_ads,
                rate_only=rate_only))
        else:
            print(chi_square(constants, fit_params, return_array=True, near_ads=near_ads,
                rate_only=rate_only)[-26:])
        print(fit_params)
        if args.debug:
            print("Observed & Predicted & Avg Predicted")
        min_chi2 = chi_square(constants, fit_params, debug=args.debug, near_ads=near_ads,
                rate_only=rate_only)
        if args.scan:
            print("Chi-square scan")
            print(sigma_searcher(fit_params, constants))
