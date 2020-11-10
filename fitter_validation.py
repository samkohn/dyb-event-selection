"""Script for running a bunch of fit configurations."""
import argparse
import itertools
import json
import multiprocessing
import os
import random
import sqlite3
import subprocess
import time

import numpy as np

import dump_LBNL_toyMC
import fit
import prediction as pred


def generate_toy(outfile_full, toy_code_dir, toy_config):
    current_dir = os.getcwd()
    os.chdir(toy_code_dir)
    with open('toymc_config_tmp.txt', 'w') as f:
        f.write(toy_config)
    if os.path.isfile(outfile_full):
        return outfile_full
    print(f'Generating toy: {outfile_full}')
    command = [
        'root', '-b', '-q', 'LoadClasses.C',
        f'genToySpectraTree.C+("toymc_config_tmp.txt", "{outfile_full}")'
    ]
    subprocess.run(command, check=True, capture_output=True)
    os.chdir(current_dir)
    return outfile_full


def main(database, label, source_category, toy_out_location, toy_code_dir,
        fit_config_template, config_template, find_errors, test_mode):
    rate_only = True
    toy_out_location = os.path.abspath(toy_out_location)
    if test_mode:
        sin2_values = np.array([0.065])
        dm2ee_values = np.array([0.0027])
    else:
        sin2_values = np.linspace(0.065, 0.09, 6)
        dm2ee_values = np.linspace(2.3e-3, 2.7e-3, 6)
    source_templates = {
        'no fluctuations default nGd binning': 'LBNL ToyMC 02 no fluctuations '
            'default nGd binning sin2={sin2} dm2ee={dm2ee}',
        'full fluctuations default nGd binning': 'LBNL ToyMC 05 w/ full fluctuations'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'far only fluctuations default nGd binning':
            'LBNL ToyMC 07 w/ stat fluctuations far, no fluctuations near (nominal),'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}'
    }
    toymc_out_numbers = {
        'no fluctuations default nGd binning': '02',
        'full fluctuations default nGd binning': '05',
        'far only fluctuations default nGd binning': '07',
    }

    source_template = source_templates[source_category]
    with open(fit_config_template, 'r') as fit_template_file:
        fit_config = json.load(fit_template_file)
    fit_file_name = 'fit_config_validation_tmp.json'
    # Read in the config file and modify it to include placeholders for
    # theta13 and dm2ee
    with open(config_template, 'r') as f:
        toy_config_template = f.read()
    toy_config_template = toy_config_template.replace(
        'sinSq2Theta13    0.084', 'sinSq2Theta13    {sin2}'
    ).replace(
        'deltaMSqee       2.48e-3', 'deltaMSqee       {dm2ee}'
    )

    NO_FLUCTUATIONS = ('02',)
    WITH_FLUCTUATIONS = ('05', '07')
    if toymc_out_numbers[source_category] in NO_FLUCTUATIONS:
        with multiprocessing.Pool(3) as pool:
            results = pool.starmap(run_validation_on_experiment, [(
                label,
                generate_toy(
                    os.path.join(
                        toy_out_location,
                        f'toy_{toymc_out_numbers[source_category]}_'
                        f'{sin2}_{dm2ee}.root'
                    ),
                    toy_code_dir,
                    toy_config_template.format(sin2=sin2, dm2ee=dm2ee),
                ),
                0,
                i,
                database,
                source_template,
                fit_config,
                fit_file_name,
                sin2,
                dm2ee,
                find_errors,
                rate_only
                ) for i, (sin2, dm2ee) in enumerate(
                    itertools.product(sin2_values, dm2ee_values)
                )
            ])
        load_to_database(database, results)
    elif toymc_out_numbers[source_category] in WITH_FLUCTUATIONS:
        for i, (sin2, dm2ee) in enumerate(
                itertools.product(sin2_values, dm2ee_values)
        ):
            print(sin2, dm2ee)
            toy_config = toy_config_template.format(sin2=sin2, dm2ee=dm2ee)
            toyfilename = generate_toy(
                os.path.join(
                    toy_out_location,
                    f'toy_{toymc_out_numbers[source_category]}_'
                    f'{sin2}_{dm2ee}.root'
                ),
                toy_code_dir,
                toy_config,
            )
            entries = range(0, 1000, 4)
            with multiprocessing.Pool(3) as pool:
                results = pool.starmap(run_validation_on_experiment, [(
                    label, toyfilename, entry, i*len(entries) + j, database,
                    source_template, fit_config, fit_file_name, sin2,
                    dm2ee, find_errors, rate_only) for j, entry in enumerate(entries)])
            load_to_database(database, results)


def load_to_database(database, results):
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        cursor.executemany('''INSERT OR REPLACE INTO fitter_validation_results
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', results)

def run_validation_on_experiment(label, toyfilename, entry, index, database,
        source_template, fit_config, fit_file_name, sin2, dm2ee, find_errors, rate_only):
    print(index)
    avg_near = True
    # Configure fit config file
    fit_config['num_coincs_source'] = source_template.format(sin2=sin2,
            dm2ee=dm2ee, entry=entry)
    fit_file_name = f'{fit_file_name}_{random.randint(0, 1000000000)}'
    with open(fit_file_name, 'w') as f:
        json.dump(fit_config, f)
    # Prepare fitter
    constants = pred.load_constants(fit_file_name)
    if rate_only:
        starting_dm2 = dm2ee
    else:
        starting_dm2 = 2.48e-3
    starting_params = pred.FitParams(
            0.15,
            starting_dm2,
            pred.ad_dict(0),
            pred.ad_dict(0, halls='near'),
            pred.core_dict(0),
            pred.ad_dict(0),
    )
    near_ads = None
    # decide whether to freeze any of the pull parameters
    # (and, if rate-only, also freeze dm2_ee)
    if rate_only:
        # frozen_params = range(1, 10)
        frozen_params = list(range(1, 29))  # 28 pull params
    else:
        # frozen_params = range(2, 10)
        frozen_params = range(2, 29)
    # Compute fit
    result = fit.fit_lsq_frozen(starting_params, constants, frozen_params,
            near_ads=near_ads, rate_only=rate_only, avg_near=avg_near)
    # Assemble best-fit FitParams object from fitter fitter output
    starting_param_list = starting_params.to_list()
    param_list = [result.x[0]]  # We know we can start with theta13
    if rate_only:
        param_list.append(dm2ee)
        first_pull_index = 1
    else:
        param_list.append(result.x[1])
        first_pull_index = 2
    for i, starting_param in enumerate(starting_param_list):
        if i < 2:
            continue
        if i in frozen_params:
            param_list.append(starting_param)
        else:
            param_list.append(result.x[i+first_pull_index])
    fit_params = pred.FitParams.from_list(param_list)
    # Compute min chi2 and sin2_2theta13 for posterity
    chi2_min = fit.chi_square(constants, fit_params, rate_only=rate_only,
            avg_near=avg_near)
    best_sin2 = np.power(np.sin(2*fit_params.theta13), 2)
    dm2ee_best = fit_params.m2_ee
    # Obtain the error on theta13 if desired
    if find_errors:
        upper, lower = fit.sigma_searcher(fit_params, constants)
        upper_sin2 = np.power(np.sin(2*upper), 2)
        lower_sin2 = np.power(np.sin(2*lower), 2)
        upper_error = upper_sin2 - best_sin2
        lower_error = best_sin2 - lower_sin2
        sin2_error = max(upper_error, lower_error)
        dm2ee_error = None
    else:
        sin2_error = None
        dm2ee_error = None
    # Delete the temporary fit config file
    os.remove(fit_file_name)
    return (label, index, sin2, best_sin2, sin2_error, dm2ee, dm2ee_best, dm2ee_error,
            chi2_min, rate_only)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('--label')
    parser.add_argument('--toy-output')
    parser.add_argument('--toy-scripts')
    parser.add_argument('--toy-config')
    parser.add_argument('--find-errors', action='store_true')
    parser.add_argument('--source-category', type=int,
        help='1: no fluctuations, 2: full fluctuations, 3: far stat only')
    parser.add_argument('--fit-config')
    parser.add_argument('--test-mode', action='store_true')
    args = parser.parse_args()
    source_categories = {
        1: 'no fluctuations default nGd binning',
        2: 'full fluctuations default nGd binning',
        3: 'far only fluctuations default nGd binning',
    }
    source_category = source_categories[args.source_category]
    main(
        args.database,
        args.label,
        source_category,
        args.toy_output,
        args.toy_scripts,
        args.fit_config,
        args.toy_config,
        args.find_errors,
        args.test_mode,
    )
