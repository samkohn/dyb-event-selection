"""Compute the error budget using the subtraction method."""

import argparse
import multiprocessing
from pprint import pprint

import numpy as np

import fit
import prediction as pred

def sigma_searcher(
    best_theta13,
    min_chi2,
    constants,
    starting_params,
    frozen_params,
    near_ads,
    rate_only,
    avg_near,
    side='both',
):
    """Find the +/- 1-sigma values given the best fit parameters.

    ``side`` could be 'upper', 'lower', or 'both'

    Returns the value(s) of theta13 at the error boundary/ies in a list
    of length 1 or 2.

    If side is 'both', the 2 entries are always returned [lower, upper].
    """
    best_sin2 = np.sin(2 * best_theta13)**2
    guess_1sigma_sin2 = 0.009
    guess_upper_sin2 = best_sin2 + guess_1sigma_sin2
    guess_lower_sin2 = best_sin2 - guess_1sigma_sin2
    guess_upper = 0.5 * np.arcsin(np.sqrt(guess_upper_sin2))
    guess_lower = 0.5 * np.arcsin(np.sqrt(guess_lower_sin2))
    min_chi2
    fit_args = (
        constants,
        starting_params.clone(),
        frozen_params,
        min_chi2,
        near_ads,
        rate_only,
        avg_near,
    )
    to_return = []
    if side in ('lower', 'both'):
        result = fit.least_squares(sigma_search_resid, guess_lower,
                args=fit_args, method='trf', xtol=1e-3, diff_step=0.00001)
        lower_limit = result.x[0]
        to_return.append(lower_limit)
    if side in ('upper', 'both'):
        result = fit.least_squares(sigma_search_resid, guess_upper,
                args=fit_args, method='trf', xtol=1e-3, diff_step=0.00001)
        upper_limit = result.x[0]
        to_return.append(upper_limit)
    return to_return

def sigma_search_resid(
    x,
    constants,
    starting_params,
    frozen_params,
    min_chi2,
    near_ads,
    rate_only,
    avg_near,
):
    """Return the difference between the delta chi-square and 1.

    x is just [theta13_search].

    The return value is *not* squared. The standard least_squares
    minimizer squares the residuals automatically.
    """
    #print('------')
    #print(f'Running fit with sin2={np.sin(2*x[0])**2}')
    fit_params = starting_params.clone()
    fit_params.theta13 = x[0]
    # Ensure that theta13 is frozen
    if 0 not in frozen_params:
        frozen_params = frozen_params.copy()
        frozen_params.append(0)
    if len(frozen_params) == len(starting_params.to_list()):
        # Fully-constrained chi-square, no parameters to fit
        chi2 = fit.chi_square(
            constants,
            fit_params,
            near_ads=near_ads,
            rate_only=rate_only,
            avg_near=avg_near,
        )
    else:
        _, fit_result = fit.fit_lsq_frozen(
            fit_params.clone(),
            constants,
            frozen_params,
            near_ads,
            rate_only,
            avg_near,
            raw_result=True,
        )
        chi2 = np.power(fit_result.fun, 2).sum()
    #print(f'Trial chi2 = {chi2:.5f}, difference={chi2 - min_chi2 - 1:.5f}')
    return chi2 - min_chi2 - 1

def run_search(
    label,
    m2_ee,
    fit_config_filename,
    frozen_params,
    near_ads,
    rate_only,
    avg_near,
    side,
):
    print(f'----Beginning computation for {label} ({side})----')
    constants = pred.load_constants(fit_config_filename)
    starting_params = pred.FitParams(0.1, m2_ee)
    # Initial fit to obtain min chi2 and best fit result
    frozen_params_free_theta13 = frozen_params.copy()
    del frozen_params_free_theta13[frozen_params.index(0)]
    best_fit_params, best_fit_results = fit.fit_lsq_frozen(
        starting_params,
        constants,
        frozen_params_free_theta13,
        near_ads=near_ads,
        rate_only=rate_only,
        avg_near=avg_near,
        raw_result=True,
    )
    best_theta13 = best_fit_params.theta13
    best_sin2 = best_fit_params.sin2_2theta13
    min_chi2 = sum((best_fit_results.fun)**2)
    bounds = sigma_searcher(
        best_theta13,
        min_chi2,
        constants,
        starting_params,
        frozen_params,
        near_ads,
        rate_only,
        avg_near,
        side=side,
    )
    print(f'----Finished computation for {label} ({side})----')
    return bounds

def effective_error(lower, upper):
    return (upper - lower) / 2

def main(fit_config_filename, m2_ee, rate_only, avg_near):
    if not rate_only:
        raise NotImplementedError("Haven't implemented rate+shape yet")
    constants = pred.load_constants(fit_config_filename)
    starting_params = pred.FitParams(0.1, m2_ee)
    near_ads = None
    frozen_params = fit.get_frozen_params(['all'], False, 'pulled')
    print('----Performing initial fit----')
    labels = [
        'All',
        'Stat. (EH3)',
        'Stat (EH1 & EH2)',
        'Reactor',
        'Efficiency',
        'Relative energy scale',
        'Accidentals',
        'Li9/He8',
        'Fast neutrons',
        'AmC',
        'Radiogenic neutrons',
        'Input mixing parameters',
    ]
    pulls = set(fit.pull_choices)
    frozen_params_list = [
        fit.get_frozen_params(['all'], True, 'pulled'),
        fit.get_frozen_params([], True, 'frozen'),
        fit.get_frozen_params(pulls - {'near-stat'}, True, 'frozen'),
        fit.get_frozen_params(pulls - {'reactor'}, True, 'frozen'),
        fit.get_frozen_params(pulls - {'efficiency'}, True, 'frozen'),
        fit.get_frozen_params(pulls - {'rel-escale'}, True, 'frozen'),
        fit.get_frozen_params(pulls - {'accidental'}, True, 'frozen'),
        fit.get_frozen_params(pulls - {'li9'}, True, 'frozen'),
        fit.get_frozen_params(pulls - {'fast-neutron'}, True, 'frozen'),
        fit.get_frozen_params(pulls - {'amc'}, True, 'frozen'),
        fit.get_frozen_params(pulls - {'rad-n'}, True, 'frozen'),
        fit.get_frozen_params(pulls - {'theta12', 'm2_21'}, True, 'pulled'),
    ]
    args_list = []
    for label, frozen_params in zip(labels, frozen_params_list):
        args_list.append([
            label,
            m2_ee,
            fit_config_filename,
            frozen_params,
            near_ads,
            rate_only,
            avg_near,
            'lower',
        ])
        args_list.append([
            label,
            m2_ee,
            fit_config_filename,
            frozen_params,
            near_ads,
            rate_only,
            avg_near,
            'upper',
        ])
    with multiprocessing.Pool() as pool:
        found_bounds_list = np.array(pool.starmap(run_search, args_list))
        print(found_bounds_list)
    # 1 row per scan: [[lower, upper], [lower, upper], ...]
    shaped_bounds = found_bounds_list.flatten().reshape((-1, 2))
    np.save('chi2_scan_error_budget_bounds_w_radn.npy', shaped_bounds)
    all_error_sq = effective_error(*shaped_bounds[0])**2
    stat_error_sq = effective_error(*shaped_bounds[1])**2
    sq_errors = {
        label: all_error_sq - effective_error(*bounds)**2
        for label, bounds in zip(labels, shaped_bounds)
    }
    sq_errors['All'] = all_error_sq
    sq_errors['Stat. (EH3)'] = stat_error_sq
    pprint(sq_errors)
    frac_errors = {
        label: err/all_error_sq for label, err in sq_errors.items()
    }
    pprint(frac_errors)

def get_simple_plus_minus(theta13, lower_bound, upper_bound, verbose=False):
    """Convert theta13 and the 68% CL bounds into xyz +abc -def for theta13 and sin2."""
    sin2 = np.sin(2 * theta13)**2
    sin2_lower = np.sin(2 * lower_bound)**2
    sin2_upper = np.sin(2 * upper_bound)**2
    theta13_pm = np.diff([lower_bound, theta13, upper_bound])
    theta13_pm[0] *= -1
    sin2_pm = np.diff([sin2_lower, sin2, sin2_upper])
    sin2_pm[0] *= -1
    if verbose:
        print(f'theta13 = {theta13}^{{+{theta13_pm[1]}}}_{{{theta13_pm[0]}}}')
        print(f'sin2(2theta13) = {sin2}^{{+{sin2_pm[1]}}}_{{{sin2_pm[0]}}}')
    return {
        'theta13': theta13,
        'theta13_pm': theta13_pm.tolist(),
        'theta13_CL': [lower_bound, upper_bound],
        'sin2': sin2,
        'sin2_pm': sin2_pm.tolist(),
        'sin2_CL': [sin2_lower, sin2_upper],
    }

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fit_config')
    parser.add_argument('--m2ee', type=float)
    parser.add_argument('--rate-only', action='store_true')
    parser.add_argument('--avg-near', action='store_true')
    args = parser.parse_args()
    main(args.fit_config, args.m2ee, args.rate_only, args.avg_near)


