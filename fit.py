import argparse
import multiprocessing
from pprint import pprint
import sqlite3
import sys

import numpy as np
from scipy.optimize import least_squares

import common
import prediction as pred

def chi_square(constants, fit_params, return_array=False, debug=False, near_ads=None,
        rate_only=False, avg_near=False, variant='poisson'):
    """Compute the chi-square value for a single set of parameters.

    Set return_array=True to return an array of terms rather than the sum.

    ``variant`` should be either ``'poisson'`` (for least-biased minimization), ``'pearson'``
    (for goodness-of-fit test), or ``'neyman'`` (included for completeness).
    """
    chi_square = 0
    if near_ads is None:
        near_ads = pred.near_ads
    far_ads = pred.far_ads
    if rate_only:
        observed = {ad: np.sum(constants.observed_candidates[ad], keepdims=True) for ad in far_ads}
    else:
        observed = {ad: constants.observed_candidates[ad] for ad in far_ads}
    predicted = pred.predict_ad_to_ad_obs(constants, fit_params)
    if rate_only:
        for key, val in predicted.items():
            predicted[key] = np.sum(val, keepdims=True)
        num_bins = 1
        num_shape_bins = 34  # TODO
    else:
        num_bins = 34
    num_pulls = fit_params.num_pulls
    if debug:
        pprint(observed)
        pprint(predicted)
    if avg_near:
        return_array_values = np.zeros((len(far_ads) * num_bins + num_pulls,))
    else:
        return_array_values = np.zeros((len(far_ads) * len(near_ads) * num_bins + num_pulls,))
    term_index = 0

    if avg_near:
        # Average the near-hall predictions
        denominator = len(near_ads)
        #  Near AD in the pairing doesn't matter but should be a member of near_ads
        predicted_avg = {(halldet, near_ads[0]): 0 for halldet in far_ads}
        for (far_halldet, near_halldet), n_predicted in predicted.items():
            if near_halldet in near_ads and far_halldet in far_ads:
                predicted_avg[far_halldet, near_ads[0]] += n_predicted/denominator
        if debug:
            pprint(predicted_avg)
        predicted = predicted_avg

    #Main part
    stat_lookup = {
        'poisson': poisson_stat_term,
        'pearson': pearson_stat_term,
        'neyman': neyman_stat_term,
    }
    stat_term_function = stat_lookup[variant]
    for (far_halldet, near_halldet), n_predicted in predicted.items():
        if far_halldet not in far_ads or near_halldet not in near_ads:
            continue
        n_observed = observed[far_halldet]
        stat_term = stat_term_function(n_observed, n_predicted)
        chi_square += np.sum(stat_term)  # sum over energy bins
        return_array_values[term_index : term_index + len(stat_term)] = stat_term
        term_index += len(stat_term)

    # Pull terms
    if 'accidental' in constants.nominal_bgs:
        acc_errors = constants.bg_errors['accidental']
        for halldet, pull in fit_params.pull_accidental.items():
            numerator = pull*pull
            denominator = acc_errors[halldet]**2
            chi_square += numerator/denominator
            return_array_values[term_index] = numerator/denominator
            term_index += 1
    if 'li9' in constants.nominal_bgs:
        li9_errors = constants.bg_errors['li9']
        for site, pull in fit_params.pull_li9.items():
            numerator = pull*pull
            # relative errors within site are the same so just pick EHX-AD1
            denominator = li9_errors[(site, 1)]**2
            chi_square += numerator/denominator
            return_array_values[term_index] = numerator/denominator
            term_index += 1
    if 'fast-neutron' in constants.nominal_bgs:
        fast_neutron_errors = constants.bg_errors['fast-neutron']
        for site, pull in fit_params.pull_fast_neutron.items():
            numerator = pull*pull
            # relative errors within site are the same so just pick EHX-AD1
            denominator = fast_neutron_errors[(site, 1)]**2
            chi_square += numerator/denominator
            return_array_values[term_index] = numerator/denominator
            term_index += 1
    if 'amc' in constants.nominal_bgs:
        # Only 1 AmC pull parameter
        amc_errors = constants.bg_errors['amc']
        numerator = fit_params.pull_amc**2
        denominator = amc_errors[1, 1]**2
        chi_square += numerator/denominator
        return_array_values[term_index] = numerator/denominator
        term_index += 1
    if 'alpha-n' in constants.nominal_bgs:
        alpha_n_errors = constants.bg_errors['alpha-n']
        for halldet, pull in fit_params.pull_alpha_n.items():
            numerator = pull*pull
            denominator = alpha_n_errors[halldet]**2
            chi_square += numerator/denominator
            return_array_values[term_index] = numerator/denominator
            term_index += 1
    for halldet, pull in fit_params.pull_near_stat.items():
        numerator = pull * pull
        denominator = 1 / constants.observed_candidates[halldet]
        pull_array = numerator/denominator
        chi_square += np.sum(pull_array)
        for element in pull_array:
            return_array_values[term_index] = element
            term_index += 1
    for core in range(1, 7):
        pull = fit_params.pull_reactor[core]
        numerator = pull * pull
        denominator = constants.reactor_err**2
        chi_square += numerator/denominator
        return_array_values[term_index] = numerator/denominator
        term_index += 1
    for halldet, pull in fit_params.pull_efficiency.items():
        numerator = pull * pull
        denominator = constants.efficiency_err**2
        chi_square += numerator/denominator
        return_array_values[term_index] = numerator/denominator
        term_index += 1
    for halldet, pull in fit_params.pull_rel_escale.items():
        numerator = pull * pull
        denominator = constants.rel_escale_err**2
        chi_square += numerator/denominator
        return_array_values[term_index] = numerator/denominator
        term_index += 1
    theta12_error = constants.theta12_err
    numerator = fit_params.pull_theta12**2
    denominator = theta12_error**2
    chi_square += numerator/denominator
    return_array_values[term_index] = numerator/denominator
    term_index += 1
    m2_21_error = constants.m2_21_err
    numerator = fit_params.pull_m2_21**2
    denominator = m2_21_error**2
    chi_square += numerator/denominator
    return_array_values[term_index] = numerator/denominator
    term_index += 1
    if return_array:
        return return_array_values
    else:
        return chi_square

def poisson_stat_term(n_observed, n_predicted):
    """Compute the Poisson maximum-likelihood chi-squared term, elementwise.

    Return a scalar if scalars are given, else return an array of the same
    shape as the inputs.

    The Poisson maximum-likelihood term is Pred - Obs + Obs * log(Obs / Pred).
    """
    diff_term = n_predicted - n_observed
    log_term = np.nan_to_num(n_observed * np.log(n_observed / n_predicted), nan=0,
            posinf=0, neginf=0)
    full_term = 2 * (diff_term + log_term)
    # Hedge against occasional floating-point errors where an element is
    # negative at the level of ~1e-13. Just make it 0 in that case.
    if any(full_term < 0):
        indices = np.nonzero(full_term < 0)
        # Print some debug info if desired
        #with np.printoptions(precision=16):
            #print(n_predicted[indices])
            #print(n_observed[indices])
            #print(diff_term[indices])
            #print(log_term[indices])
            #print(full_term[indices])
            #print(np.transpose(indices))
        full_term[indices] = 0
    return full_term

def pearson_stat_term(n_observed, n_predicted):
    """Compute the Pearson chi-squared term, elementwise.

    Return a scalar if scalars are given, else return an array of the same
    shape as the inputs.

    The Pearson chi-squared term is (Pred - Obs)**2 / Pred.
    """
    numerator = np.power(n_observed - n_predicted, 2)
    denominator = n_predicted
    ratio = numerator/denominator
    return ratio

def neyman_stat_term(n_observed, n_predicted):
    """Compute the Neyman chi-squared term, elementwise.

    Return a scalar if scalars are given, else return an array of the same
    shape as the inputs.

    The Neyman chi-squared term is (Pred - Obs)**2 / Obs.
    """
    numerator = np.power(n_observed - n_predicted, 2)
    denominator = n_observed
    ratio = numerator/denominator
    return ratio

def residual_fn(x, constants, near_ads=None, rate_only=False, avg_near=False):
    """Convert arguments from scipy to desired format and return the residuals.
    """
    fit_params = pred.FitParams.from_list(x)
    # Take the square root because the fitter wants the linear residuals
    # and squares them on its own...
    residuals = np.sqrt(chi_square(constants, fit_params, return_array=True,
        near_ads=near_ads, rate_only=rate_only, avg_near=avg_near,
        variant='poisson'))
    return residuals

def residual_frozen_param(frozen_dict, near_ads, rate_only, avg_near):
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
        return residual_fn(
            param_list, constants, near_ads=near_ads, rate_only=rate_only,
            avg_near=avg_near
        )
    return residual

def fit_lsq_frozen(starting_params, constants, frozen_params, near_ads, rate_only,
        avg_near, raw_result=False):
    """Perform the fit with certain parameters frozen."""
    frozen_params_dict = {}
    all_params = starting_params.to_list()
    x0 = []
    reverse_map = {}
    num_frozen_so_far = 0
    for i, param in enumerate(all_params):
        if i in frozen_params:
            frozen_params_dict[i] = param
        else:
            x0.append(param)
            reverse_map[i] = num_frozen_so_far
            num_frozen_so_far += 1
    residual = residual_frozen_param(frozen_params_dict, near_ads, rate_only, avg_near)
    result = least_squares(residual, x0, args=(constants,), method='trf')
    # Assemble best-fit FitParams object from fitter fitter output
    starting_param_list = starting_params.to_list()
    result_param_list = []
    for i, starting_param in enumerate(starting_param_list):
        if i in reverse_map:
            result_param_list.append(result.x[reverse_map[i]])
        else:
            result_param_list.append(starting_param)
    fit_params = pred.FitParams.from_list(result_param_list)

    if raw_result:
        return (fit_params, result)
    else:
        return fit_params


def sigma_searcher(fit_params, constants, side='both'):
    """Find the +/- 1-sigma values given the best fit parameters.

    ``side`` could be 'upper', 'lower', or 'both'

    Returns the value(s) of theta13 at the error boundary/ies in a list
    of length 1 or 2.

    If side is 'both', the 2 entries are always returned [upper, lower]
    """
    best_theta13 = fit_params.theta13
    best_sin2 = fit_params.sin2_2theta13
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

def save_result(
    database,
    description,
    constants,
    fit_params,
    rate_only,
    avg_near,
    sin2_error_plus=None,
    sin2_error_minus=None,
    m2_ee_error_plus=None,
    m2_ee_error_minus=None,
):
    """Save the specified results to the database.

    The FitConstants object is not saved to the db.
    It is only used to compute the chi-squareds.
    """
    chi2_poisson = chi_square(
        constants, fit_params, rate_only=rate_only, avg_near=avg_near
    )
    chi2_pearson = chi_square(
        constants, fit_params, rate_only=rate_only, avg_near=avg_near, variant='pearson'
    )
    fits_db_row = (
        description,
        int(rate_only),
        int(avg_near),
        fit_params.sin2_2theta13,
        fit_params.m2_ee,
        chi2_poisson,
        chi2_pearson,
        sin2_error_plus,
        sin2_error_minus,
        m2_ee_error_plus,
        m2_ee_error_minus,
        constants.input_osc_params.theta12,
        constants.input_osc_params.m2_21,
        int(constants.input_osc_params.hierarchy),
        constants.input_osc_params.m2_ee_conversion,
    )
    params_list = fit_params.to_list()
    params_rows = []
    for name, index in fit_params.index_map().items():
        if name in ('theta13', 'm2_ee'):
            # The best-fit values get saved separately
            continue
        if isinstance(index, int):
            # Only 1 value, not an array
            param_value = params_list[index]
            params_rows.append((name, None, param_value))
        else:
            # Multiple params in a list, index is a slice
            param_values = params_list[index]
            for param_index, param_value in enumerate(param_values):
                params_rows.append((name, param_index, param_value))
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute(f'''
            INSERT INTO
                fits (
                    Description,
                    IsRateOnly,
                    IsAvgNear,
                    FitSinSqT13,
                    FitDM2ee,
                    ChiSquareFit,
                    ChiSquareGof,
                    SinSqT13_ErrorPlus,
                    SinSqT13_ErrorMinus,
                    DM2ee_ErrorPlus,
                    DM2ee_ErrorMinus,
                    Theta12,
                    DM2_21,
                    Hierarchy,
                    DM2ee_conversion
                )
            VALUES
                ({", ".join("?"*15)}); -- 15 parameters
            ''',
            fits_db_row,
        )
        fit_id = cursor.lastrowid
        # Update params_rows to include the fit_id
        real_params_rows = []
        for row in params_rows:
            real_params_rows.append((fit_id, *row))
        cursor.executemany('''
            INSERT INTO
                pulls
            VALUES
                (?, ?, ?, ?)
            ''',
            real_params_rows,
        )
    return

def load_result(database, id_or_description):
    """Load saved results into (FitParams, dict of other attributes)."""
    if isinstance(id_or_description, int):
        use_id = True
    else:
        use_id = False
    with common.get_db(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute(f'''
            SELECT
                *
            FROM
                fits
            WHERE
                {"Id" if use_id else "Description"} = ?
            LIMIT
                1
            ''',
            (id_or_description,),
        )
        fit_info = dict(cursor.fetchone())
        fit_info['IsRateOnly'] = bool(fit_info['IsRateOnly'])
        fit_info['IsAvgNear'] = bool(fit_info['IsAvgNear'])
        fit_info['Hierarchy'] = bool(fit_info['Hierarchy'])
        fit_id = fit_info['Id']
        fit_theta13 = 0.5 * np.arcsin(np.sqrt(fit_info['FitSinSqT13']))
        cursor.execute('''
            SELECT
                ParamName,
                ParamIndex,
                ParamValue
            FROM
                pulls
            WHERE
                FitId = ?
            ORDER BY
                ParamName,
                ParamIndex
            ''',
            (fit_id,),
        )
        rows = cursor.fetchall()
        param_list = [None] * (2 + pred.FitParams.num_pulls)  # 2 for theta13, dm2
        index_map = pred.FitParams.index_map()
        param_list[index_map['theta13']] = fit_theta13
        param_list[index_map['m2_ee']] = fit_info['FitDM2ee']
        for row in rows:
            index = index_map[row['ParamName']]
            if isinstance(index, int):
                param_list[index] = row['ParamValue']
            else:
                offset = row['ParamIndex']
                real_index = index.start + offset
                param_list[real_index] = row['ParamValue']
        fit_params = pred.FitParams.from_list(param_list)
        return fit_params, fit_info


def grid(
    theta13_values,
    constants,
    starting_params,
    frozen_params,
    near_ads,
    rate_only,
    avg_near,
):
    """Run fits with a list of different theta13 starting values.

    The starting_params object is modified in-place but restored to its
    original state at the end of the function.

    fit_lsq_frozen_kwargs are passed directly to fit_lsq_frozen.
    """
    def fit_args_generator():
        for theta13_value in theta13_values:
            new_params = starting_params.clone()
            new_params.theta13 = theta13_value
            yield (
                new_params,
                constants,
                frozen_params,
                near_ads,
                rate_only,
                avg_near,
            )
        return
    def chi_square_args_generator(fit_params_list):
        for fit_params_value in fit_params_list:
            yield (
                constants,
                fit_params_value,
                False,
                False,
                near_ads,
                rate_only,
                avg_near,
                'poisson',
            )
        return
    with multiprocessing.Pool() as pool:
        fit_param_results = pool.starmap(fit_lsq_frozen, fit_args_generator())
        min_chisquares = pool.starmap(
            chi_square,
            chi_square_args_generator(fit_param_results),
        )
    return min_chisquares

def get_frozen_params(pulls, freeze_theta13, freeze_dm2ee, near_ads=None):
    """Get the list of frozen params given the fitter configuration.

    A list of near ADs can be provided, which will be used to freeze the
    near-stat pull parameters of any near AD that is not being used in
    the fit.
    """
    dummy_params = pred.FitParams(0, 0)
    index_map = dummy_params.index_map()
    # decide whether to freeze any of the pull parameters
    # (and, if rate-only, also freeze dm2_ee)
    frozen_params = []
    if freeze_theta13:
        frozen_params.append(index_map['theta13'])
    if freeze_dm2ee:
        frozen_params.append(index_map['m2_ee'])
    if 'all' in pulls:
        pass
    elif len(pulls) == 0:
        n_total_params = len(dummy_params.to_list())
        frozen_params.extend(list(range(2, n_total_params)))
    else:
        def freeze(name):
            if name in ('theta12', 'm2_21'):
                name = 'pull_' + name
            pull_slice = index_map[name]
            if isinstance(pull_slice, slice):  # range of values
                indices = list(range(pull_slice.start, pull_slice.stop))
                frozen_params.extend(indices)
            elif isinstance(pull_slice, int):  # just a single value
                frozen_params.append(pull_slice)
            else:
                raise ValueError(f"Can't parse starting_params.index_map()[{name!r}]")
        for pull_name in pull_choices:
            if pull_name not in pulls:
                freeze(pull_name)
        if 'near-stat' in pulls and near_ads is not None and len(near_ads) != 4:
            # Ensure that ADs not being counted are frozen out.
            # Normally I don't care but there are so many pulls that it
            # slows down the fitter.
            pulls_slice = index_map['near_stat']
            first = pulls_slice.start
            last = pulls_slice.stop
            n_near_pulls = last - first
            n_bins = n_near_pulls // 4
            for i, halldet in enumerate(pred.near_ads):
                if halldet not in near_ads:
                    first_for_ad = first + i * n_bins
                    last_for_ad = first + (i + 1) * n_bins
                    frozen_params.extend(list(range(first_for_ad, last_for_ad)))
    return frozen_params





pull_choices = (
    'reactor',
    'efficiency',
    'rel-escale',
    'accidental',
    'li9',
    'fast-neutron',
    'amc',
    'alpha-n',
    'near-stat',
    'theta12',
    'm2_21',
)
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
    parser.add_argument("--avg-near", action='store_true')
    parser.add_argument("--freeze-sin2-2theta13", type=float, default=None)
    parser.add_argument(
        "--pulls", nargs='*', choices=pull_choices+('all',), default=[]
    )
    parser.add_argument("--save-db", help="db file to save the results")
    parser.add_argument("--save-descr", help="verbose description to be saved")
    args = parser.parse_args()
    if (args.save_db is None) != (args.save_descr is None):
        raise ValueError("Must specify both or neither of --save-db and --save-descr")
    rate_only = not args.shape
    pulls = args.pulls
    near_ads = None
    constants = pred.load_constants(args.config)
    starting_theta13 = 0.15  # Near-ish to expected fit value
    if args.dm2ee is None and rate_only:
        raise ValueError("Can't fit dm2_ee in rate-only")
    elif rate_only or args.dm2ee is not None:
        starting_dm2 = args.dm2ee
        freeze_dm2ee = True
    else:
        starting_dm2 = 2.48e-3
        freeze_dm2ee = False
    if args.freeze_sin2_2theta13 is not None:
        starting_theta13 = 0.5 * np.arcsin(np.sqrt(args.freeze_sin2_2theta13))
        freeze_theta13 = True
    else:
        freeze_theta13 = False
    starting_params = pred.FitParams(
        starting_theta13,
        starting_dm2,
    )
    frozen_params = get_frozen_params(
        pulls,
        freeze_theta13,
        freeze_dm2ee,
        near_ads,
    )
    print(frozen_params)
    if args.debug:
        print(starting_params)
    print(chi_square(constants, starting_params, rate_only=rate_only, debug=args.debug,
        avg_near=args.avg_near, near_ads=near_ads, variant='poisson'))
    if args.debug:
        print(chi_square(constants, starting_params, return_array=True,
            rate_only=rate_only, avg_near=args.avg_near, near_ads=near_ads,
            variant='poisson'))
    if not args.no_fit:
        fit_params, result = fit_lsq_frozen(starting_params, constants, frozen_params,
                near_ads=near_ads, rate_only=rate_only, avg_near=args.avg_near,
                raw_result=True)
        #print(repr(result.x))
        print('sin22theta13 =', fit_params.sin2_2theta13)
        print('dm2_ee = ', fit_params.m2_ee)
        print(f'Num function evals: {result.nfev}')
        print(f'Num jacobian evals: {result.njev}')
        print(result.message)
        if not result.success:
            sys.exit(0)

        print('Min chi-square:', chi_square(constants, fit_params, return_array=False, near_ads=near_ads,
            rate_only=rate_only, avg_near=args.avg_near, variant='poisson'))
        if args.debug:
            print(chi_square(constants, fit_params, return_array=True, near_ads=near_ads,
                rate_only=rate_only, avg_near=args.avg_near, variant='poisson'))
        #print(fit_params)
        if args.debug:
            print("Observed & Predicted & Avg Predicted")
        min_chi2 = chi_square(constants, fit_params, debug=args.debug, near_ads=near_ads,
                rate_only=rate_only, avg_near=args.avg_near, variant='poisson')
        if args.scan:
            print("Chi-square scan")
            print(sigma_searcher(fit_params, constants))
        if args.save_descr is not None:
            save_result(
                args.save_db,
                args.save_descr,
                constants,
                fit_params,
                rate_only,
                args.avg_near,
            )
