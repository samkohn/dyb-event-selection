import argparse
from pprint import pprint
import sys

import numpy as np
from scipy.optimize import least_squares

import prediction as pred

def chi_square(constants, fit_params, return_array=False, debug=False):
    """Compute the chi-square value for a single set of parameters.

    Set return_array=True to return an array of terms rather than the sum.
    """
    chi_square = 0
    observed = {ad: constants.observed_candidates[ad] for ad in pred.far_ads}
    predicted, reco_bins = pred.predict_ad_to_ad_obs(constants, fit_params)
    if debug:
        pprint(observed)
        pprint(predicted)
    return_array_values = np.zeros((30,))
    term_index = 0

    # Average the near-hall predictions
    predicted_avg = pred.ad_dict(0, halls='far')
    for (far_halldet, near_halldet), n_predicted in predicted.items():
        predicted_avg[far_halldet] += n_predicted/4
    if debug:
        pprint(predicted_avg)

    #Main part
    for far_halldet, n_predicted in predicted_avg.items():
        n_observed = observed[far_halldet]
        sigma_observed = np.sqrt(n_observed)  # TODO placeholder error
        numerator = np.power(n_observed - n_predicted, 2)
        denominator = np.power(sigma_observed, 2)
        chi_square += np.sum(numerator/denominator)  # sum over energy bins
        return_array_values[term_index] = np.sum(numerator/denominator)
        term_index += 1

    # Pull terms
    for halldet, pull in fit_params.pull_bg.items():
        numerator = pull*pull
        denominator = 0.010*0.010  # TODO relative error on accidentals
        chi_square += numerator/denominator
        return_array_values[term_index] = numerator/denominator
        term_index += 1
    for halldet, pull in fit_params.pull_near_stat.items():
        numerator = pull * pull
        denominator = np.power(constants.observed_candidates[halldet], -1)
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

def residual_fn(x, constants):
    """Convert arguments from scipy to desired format and return the residuals.
    """
    fit_params = pred.FitParams.from_list(x)
    # Take the square root because the fitter wants the linear residuals
    # and squares them on its own...
    residuals = np.sqrt(chi_square(constants, fit_params, return_array=True))
    return residuals

def residual_frozen_param(frozen_dict):
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
        return residual_fn(param_list, constants)
    return residual

def fit_lsq_frozen(starting_params, constants, frozen_params):
    """Perform the fit with certain parameters frozen."""
    frozen_params_dict = {}
    all_params = starting_params.to_list()
    x0 = []
    for i, param in enumerate(all_params):
        if i in frozen_params:
            frozen_params_dict[i] = param
        else:
            x0.append(param)
    residual = residual_frozen_param(frozen_params_dict)
    result = least_squares(residual, x0, args=(constants,), method='trf')
    return result

def chi_square_grid(starting_params, constants, theta13_values):
    """Return a grid of optimal chi-square values for the given parameter grid.
    """
    best_theta13 = starting_params.theta13
    result = np.zeros((len(theta13_values)))
    for i, theta13 in enumerate(theta13_values):
        starting_params.theta13 = theta13
        fit_result = fit_lsq_frozen(starting_params, constants, [0])
        result[i] = np.power(fit_result.fun, 2).sum()
    starting_params.theta13 = best_theta13  # restore original value
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("--no-fit", action='store_true')
    parser.add_argument("--scan", action='store_true')
    args = parser.parse_args()
    constants = pred.load_constants(args.config)
    starting_params = pred.FitParams(
            0.15,
            pred.ad_dict(0),
            pred.ad_dict(0, halls='near'),
            pred.core_dict(0),
            pred.ad_dict(0),
    )
    print(starting_params)
    print(chi_square(constants, starting_params))
    print(chi_square(constants, starting_params, return_array=True))
    if not args.no_fit:
        result = fit_lsq_frozen(starting_params, constants, {})
        print(repr(result.x))
        print('sin22theta13 =', np.power(np.sin(2*result.x[0]), 2))
        print(result.success)
        print(result.message)
        if not result.success:
            sys.exit(0)
        fit_params = pred.FitParams.from_list(result.x)
        print(chi_square(constants, fit_params, return_array=False))
        print(chi_square(constants, fit_params, return_array=True))
        print(fit_params)
        print("Observed & Predicted & Avg Predicted")
        chi_square(constants, fit_params, debug=True)
        if args.scan:
            print("Chi-square scan")
            theta13 = fit_params.theta13
            low_value = theta13 - 0.005
            up_value = theta13 + 0.005
            values = np.linspace(low_value, up_value, 15)
            grid = (chi_square_grid(fit_params, constants, values))
            pprint(dict(zip(values, grid)))
            pprint(dict(zip(np.power(np.sin(2*values), 2), np.sqrt(grid - np.min(grid)))))
