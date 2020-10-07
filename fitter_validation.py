"""Script for running a bunch of fit configurations."""
import argparse
import itertools
import json
import os
import sqlite3
import subprocess

import numpy as np

import dump_LBNL_toyMC
import fit
import prediction as pred


def generate_toy(outfilename, toy_code_dir, toy_config, sin2, dm2ee):
    current_dir = os.getcwd()
    os.chdir(toy_code_dir)
    with open('toymc_config_tmp.txt', 'w') as f:
        f.write(toy_config)
    outfile_full = os.path.join(outfilename, f'toy_01_{sin2}_{dm2ee}.root')
    if os.path.isfile(outfile_full):
        return outfile_full
    command = [
        'root', '-b', '-q', 'LoadClasses.C',
        f'genToySpectraTree.C+("toymc_config_tmp.txt", "{outfile_full}")'
    ]
    subprocess.run(command, check=True)
    os.chdir(current_dir)
    return outfile_full


def main(database, label, toy_out_location, toy_code_dir, config_template):
    sin2_values = np.linspace(0.065, 0.09, 11)
    dm2ee_values = np.linspace(2.3e-3, 2.7e-3, 11)
    source_template = "LBNL ToyMC no fluctuations sin2={sin2} dm2ee={dm2ee}"
    fit_config = {
        "database": "/home/skohn/parameters_new.db",
        "period": "8ad",
        "backgrounds": False,
        "backgrounds_source": "Nominal rate-only 9/17/2020",
        "mult_eff": [0.9770, 0.9774, 0.9782, 0.9781, 0.9785, 0.9787, 0.9785, 0.9784],
        "muon_eff": [0.8233, 0.8196, 0.8516, 0.8501, 0.9832, 0.9830, 0.9828, 0.9830],
        "masses": [19941, 19966, 19891, 19945, 19913, 19991, 19892, 19931],
        "num_coincs": True,
        "num_coincs_source": None, # This is what we'll replace in the loop
        "reco_bins": None,
        "det_response_source": "THU ToyMC res_p:Ev No Cuts nominal binning",
    }
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

    for i, (sin2, dm2ee) in enumerate(itertools.product(sin2_values, dm2ee_values)):
        toy_config = toy_config_template.format(sin2=sin2, dm2ee=dm2ee)
        toyfilename = generate_toy(toy_out_location, toy_code_dir, toy_config,
                sin2, dm2ee)
        dump_LBNL_toyMC.main(
            toyfilename,
            0,
            database,
            source_template.format(sin2=sin2, dm2ee=dm2ee),
            "nH nominal",
        )
        fit_config['num_coincs_source'] = source_template.format(sin2=sin2, dm2ee=dm2ee)
        with open(fit_file_name, 'w') as f:
            json.dump(fit_config, f)
        constants = pred.load_constants(fit_file_name)
        constants.input_osc_params.m2_ee = dm2ee
        starting_params = pred.FitParams(
                0.15,
                pred.ad_dict(0),
                pred.ad_dict(0, halls='near'),
                pred.core_dict(0),
                pred.ad_dict(0),
        )
        near_ads = None
        result = fit.fit_lsq_frozen(starting_params, constants, range(1, 9),
                near_ads=near_ads)
        fit_params = pred.FitParams.from_list(
                [result.x[0]] + [0] * 8 + result.x[1:].tolist()
        )
        chi2_min = fit.chi_square(constants, fit_params)
        upper, lower = fit.sigma_searcher(fit_params, constants)
        best_sin2 = np.power(np.sin(2*fit_params.theta13), 2)
        upper_sin2 = np.power(np.sin(2*upper), 2)
        lower_sin2 = np.power(np.sin(2*lower), 2)
        upper_error = upper_sin2 - best_sin2
        lower_error = best_sin2 - lower_sin2
        sin2_error = max(upper_error, lower_error)
        dm2ee_best = None
        dm2ee_error = None
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''INSERT OR REPLACE INTO fitter_validation_results
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
            (label, i, sin2, best_sin2, sin2_error, dm2ee, dm2ee_best, dm2ee_error,
                chi2_min, True))






if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('--label')
    parser.add_argument('--toy-output')
    parser.add_argument('--toy-scripts')
    parser.add_argument('--toy-config')
    args = parser.parse_args()
    main(args.database, args.label, args.toy_output, args.toy_scripts, args.toy_config)
