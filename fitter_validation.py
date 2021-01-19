"""Script for running a bunch of fit configurations."""
import argparse
import itertools
import json
import multiprocessing
import os
import random
import subprocess
import time

import numpy as np

import common
import dump_LBNL_toyMC
import fit
import prediction as pred


def generate_toy(outfile_full, toy_code_dir, toy_config, sin2, dm2ee):
    if os.path.isfile(outfile_full):
        return outfile_full
    current_dir = os.getcwd()
    os.chdir(toy_code_dir)
    with open('toymc_config_tmp.txt', 'w') as f:
        f.write(toy_config)
    print(f'Generating toy: {outfile_full}')
    command = [
        'root', '-b', '-q', 'LoadClasses.C',
        f'genToySpectraTree.C+("toymc_config_tmp.txt", "{outfile_full}", {sin2}, {dm2ee})'
    ]
    try:
        subprocess.run(command, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.cmd)
        print('-------STDOUT-------')
        print(e.stdout[-2000:])
        print('-------STDERR-------')
        print(e.stderr[-2000:])
        raise
    os.chdir(current_dir)
    return outfile_full


def main(database, label, source_category, toy_out_location, toy_code_dir,
        fit_config_template, config_template, find_errors, test_mode,
        gen_mc_only, dump_mc, rate_only, avg_near, pulls, stage
):
    num_multiprocessing = 64
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
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'near+far stat fluctuations default nGd binning':
            'LBNL ToyMC 08 w/ near+far stat fluctuations, default nGd binning '
            'sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'near only fluctuations default nGd binning':
            'LBNL ToyMC 07 w/ stat fluctuations near, no fluctuations far (nominal),'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'relative energy scale default nGd binning':
            'LBNL ToyMC 10 no stat fluctuations, fluctuate relative energy scale,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'core spectra default nGd binning':
            'LBNL ToyMC 11 no stat fluctuations, fluctuate core spectra,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'reactor power default nGd binning':
            'LBNL ToyMC 12 no stat fluctuations, fluctuate reactor power,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'no fluctuations w acc bg default nGd binning':
            'LBNL ToyMC 13 no stat fluctuations, w/ accidental background,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee}',
        'acc bg fluctuations default nGd binning':
            'LBNL ToyMC 14 accidental bg fluctuations, no IBD stat fluctuations,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'no fluctuations w li9 bg default nGd binning':
            'LBNL ToyMC 15 no stat fluctuations, w/ li9 background,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee}',
        'li9 bg fluctuations default nGd binning':
            'LBNL ToyMC 16 li9 bg fluctuations, no IBD stat fluctuations,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'no fluctuations w all all bg default nGd binning':
            'LBNL ToyMC 17 no stat fluctuations, w/ all backgrounds,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee}',
        'all bg fluctuations default nGd binning':
            'LBNL ToyMC 18 all bg fluctuations, no IBD stat fluctuations,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'efficiency default nGd binning':
            'LBNL ToyMC 19 no stat fluctuations, fluctuate detection efficiency,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'no fluctuations w fastn bg default nGd binning':
            'LBNL ToyMC 20 no stat fluctuations, w/ fastn background,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee}',
        'fastn bg fluctuations default nGd binning':
            'LBNL ToyMC 21 fastn bg fluctuations, no IBD stat fluctuations,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'no fluctuations w amc bg default nGd binning':
            'LBNL ToyMC 22 no stat fluctuations, w/ amc background,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee}',
        'amc bg fluctuations default nGd binning':
            'LBNL ToyMC 23 amc bg fluctuations, no IBD stat fluctuations,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'resolution default nGd binning':
            'LBNL ToyMC 24 no stat fluctuations, fluctuate energy resolution,'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'scint nonlinearity default nGd binning':
            'LBNL ToyMC 25 no stat fluctuations, fluctuate scint nonlinearity'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'all sys and stat default nGd binning':
            'LBNL ToyMC 26 all systematics and stat fluctuations'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'all sys and stat except solar osc default nGd binning':
            'LBNL ToyMC 27 all systematics and stat fluctuations except solar osc'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
        'all syst no stat default nGd binning':
            'LBNL ToyMC 28 all systematic fluctuations, no stat flucts'
            ' default nGd binning sin2={sin2} dm2ee={dm2ee} experiment #{entry}',
    }
    toymc_out_numbers = {
        'no fluctuations default nGd binning': '02',
        'full fluctuations default nGd binning': '05',
        'far only fluctuations default nGd binning': '07',
        'near+far stat fluctuations default nGd binning': '08',
        'near only fluctuations default nGd binning': '07',
        'relative energy scale default nGd binning': '10',
        'core spectra default nGd binning': '11',
        'reactor power default nGd binning': '12',
        'no fluctuations w acc bg default nGd binning': '13',
        'acc bg fluctuations default nGd binning': '14',
        'no fluctuations w li9 bg default nGd binning': '15',
        'li9 bg fluctuations default nGd binning': '16',
        'no fluctuations w all all bg default nGd binning': '17',
        'all bg fluctuations default nGd binning': '18',
        'efficiency default nGd binning': '19',
        'no fluctuations w fastn bg default nGd binning': '20',
        'fastn bg fluctuations default nGd binning': '21',
        'no fluctuations w amc bg default nGd binning': '22',
        'amc bg fluctuations default nGd binning': '23',
        'resolution default nGd binning': '24',
        'scint nonlinearity default nGd binning': '25',
        'all sys and stat default nGd binning': '26',
        'all sys and stat except solar osc default nGd binning': '27',
        'all syst no stat default nGd binning': '28',
    }
    mc_configurations = {
        # (Has far stat fluctuations, has near stats, has any systematic fluctuations)
        'no fluctuations default nGd binning': (False, False, False),
        'full fluctuations default nGd binning': (True, True, True),
        'far only fluctuations default nGd binning': (True, False, False),
        'near+far stat fluctuations default nGd binning': (True, True, False),
        'near only fluctuations default nGd binning': (False, True, False),
        'relative energy scale default nGd binning': (False, False, True),
        'core spectra default nGd binning': (False, False, True),
        'reactor power default nGd binning': (False, False, True),
        'no fluctuations w acc bg default nGd binning': (False, False, False),
        'acc bg fluctuations default nGd binning': (False, False, True),
        'no fluctuations w li9 bg default nGd binning': (False, False, False),
        'li9 bg fluctuations default nGd binning': (False, False, True),
        'no fluctuations w all all bg default nGd binning': (False, False, True),
        'all bg fluctuations default nGd binning': (False, False, True),
        'efficiency default nGd binning': (False, False, True),
        'no fluctuations w fastn bg default nGd binning': (False, False, False),
        'fastn bg fluctuations default nGd binning': (False, False, True),
        'no fluctuations w amc bg default nGd binning': (False, False, False),
        'amc bg fluctuations default nGd binning': (False, False, True),
        'resolution default nGd binning': (False, False, True),
        'scint nonlinearity default nGd binning': (False, False, True),
        'all sys and stat default nGd binning': (True, True, True),
        'all sys and stat except solar osc default nGd binning': (True, True, True),
        'all syst no stat default nGd binning': (False, False, True),
    }
    mc_configuration = mc_configurations[source_category]
    if mc_configuration == (False, False, True):
        nominal_near = False
        nominal_far = False
    else:
        nominal_near = not mc_configuration[1]
        nominal_far = not mc_configuration[0]
        if nominal_near and mc_configuration[2]:
            raise ValueError("Can't figure out how to suppress near fluctuations but add "
                    "systematics")


    source_template = source_templates[source_category]
    if config_template is not None:
        # Read in the config file and modify it to include placeholders for
        # theta13 and dm2ee
        with open(config_template, 'r') as f:
            toy_config_template = f.read()
        toy_config_template = toy_config_template.replace(
            'sinSq2Theta13    0.084', 'sinSq2Theta13    {sin2}'
        ).replace(
            'deltaMSqee       2.48e-3', 'deltaMSqee       {dm2ee}'
        )
    else:
        toy_config_template = ''

    NO_FLUCTUATIONS = ('02', '13', '15', '17', '20', '22')
    WITH_FLUCTUATIONS = ('05', '07', '08', '10', '11', '12', '14', '16', '18', '19', '21',
            '23', '24', '25', '26', '27', '28')
    # Get a fresh generator when needed
    grid_values = lambda: itertools.product(sin2_values, dm2ee_values)
    if toymc_out_numbers[source_category] in NO_FLUCTUATIONS:
        with multiprocessing.Pool(num_multiprocessing) as pool:
            toyfilenames = pool.starmap(generate_toymc_files, [(
                    toy_config_template,
                    toy_out_location,
                    sin2,
                    dm2ee,
                    toymc_out_numbers[source_category],
                    toy_code_dir,
                ) for sin2, dm2ee in grid_values()]
            )
        if dump_mc:
            # Not using multiprocessing so as not to overload db
            for toyfilename, (sin2, dm2ee) in zip(toyfilenames, grid_values()):
                dump_LBNL_toyMC.main(
                    toyfilename,
                    0,
                    database,
                    source_template.format(sin2=sin2, dm2ee=dm2ee),
                    "default",
                    nominal_near,
                    nominal_far,
                    stage,
                )
        if not gen_mc_only:
            with open(fit_config_template, 'r') as fit_template_file:
                fit_config = json.load(fit_template_file)
            fit_file_name = 'fit_config_validation_tmp.json'
            with multiprocessing.Pool(num_multiprocessing) as pool:
                results = pool.starmap(run_validation_on_experiment, [(
                    label,
                    toyfilename,
                    0,
                    i,
                    database,
                    source_template,
                    fit_config,
                    fit_file_name,
                    sin2,
                    dm2ee,
                    find_errors,
                    rate_only,
                    avg_near,
                    pulls,
                    ) for i, (toyfilename, (sin2, dm2ee)) in enumerate(
                        zip(toyfilenames, grid_values())
                    )
                ])
            load_to_database(database, results, mc_configuration)
    elif toymc_out_numbers[source_category] in WITH_FLUCTUATIONS:
        with multiprocessing.Pool(num_multiprocessing) as pool:
            toyfilenames = pool.starmap(generate_toymc_files, [(
                    toy_config_template,
                    toy_out_location,
                    sin2,
                    dm2ee,
                    toymc_out_numbers[source_category],
                    toy_code_dir,
                ) for sin2, dm2ee in grid_values()]
            )
        for i, (toyfilename, (sin2, dm2ee)) in enumerate(
            zip(toyfilenames, grid_values())
        ):
            entries = range(0, 1000)
            if dump_mc:
                sources = []
                for entry in entries:
                    sources.append(
                        source_template.format(sin2=sin2, dm2ee=dm2ee, entry=entry)
                    )
                dump_LBNL_toyMC.multiple(
                    toyfilename,
                    entries,
                    database,
                    sources,
                    "default",
                    nominal_near,
                    nominal_far,
                    stage,
                )
            print(sin2, dm2ee)
            if gen_mc_only:
                continue
            with open(fit_config_template, 'r') as fit_template_file:
                fit_config = json.load(fit_template_file)
            fit_file_name = 'fit_config_validation_tmp.json'
            with multiprocessing.Pool(num_multiprocessing) as pool:
                results = pool.starmap(run_validation_on_experiment, [(
                    label, toyfilename, entry, i*len(entries) + j, database,
                    source_template, fit_config, fit_file_name, sin2,
                    dm2ee, find_errors, rate_only, avg_near,
                    pulls) for j, entry in enumerate(entries)])
            load_to_database(database, results, mc_configuration)

def load_to_database(database, results, mc_configuration):
    extended_results = []
    for row in results:
        new_row = tuple(row) + tuple(mc_configuration)
        extended_results.append(new_row)
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.executemany('''INSERT OR REPLACE INTO fitter_validation_results
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', extended_results)
    return

def generate_toymc_files(toy_config_template, toy_out_location, sin2, dm2ee, toy_out_number, toy_code_dir):
    toy_config = toy_config_template.format(sin2=sin2, dm2ee=dm2ee)
    toy_outfile = os.path.join(
        toy_out_location,
        f'toy_{toy_out_number}_{sin2}_{dm2ee}.root'
    )
    toyfilename = generate_toy(
        toy_outfile,
        toy_code_dir,
        toy_config,
        sin2,
        dm2ee,
    )
    return toyfilename


def run_validation_on_experiment(label, toyfilename, entry, index, database,
        source_template, fit_config, fit_file_name, sin2, dm2ee, find_errors,
        rate_only, avg_near, pulls):
    # Configure fit config file
    fit_config['num_coincs_source'] = source_template.format(sin2=sin2,
            dm2ee=dm2ee, entry=entry)
    fit_file_name = f'{fit_file_name}_{random.randint(0, 1000000000)}'
    with open(fit_file_name, 'w') as f:
        json.dump(fit_config, f)
    try:
        # Prepare fitter
        constants = pred.load_constants(fit_file_name)
        if rate_only:
            starting_dm2 = dm2ee
        else:
            starting_dm2 = 2.48e-3
        starting_params = pred.FitParams(
            0.15,
            starting_dm2,
        )
        near_ads = None
        n_total_params = len(starting_params.to_list())
        # decide whether to freeze any of the pull parameters
        # (and, if rate-only, also freeze dm2_ee)
        if rate_only:
            frozen_params = [1]
        else:
            frozen_params = []
        if 'all' in pulls:
            pass
        elif len(pulls) == 0:
            frozen_params.extend(list(range(2, n_total_params)))
        else:
            index_map = starting_params.index_map()
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

        # Compute fit
        fit_params = fit.fit_lsq_frozen(starting_params, constants, frozen_params,
                near_ads=near_ads, rate_only=rate_only, avg_near=avg_near)
        # Compute min chi2 and sin2_2theta13 for posterity
        chi2_min = fit.chi_square(constants, fit_params, rate_only=rate_only,
                avg_near=avg_near, near_ads=near_ads, variant='poisson')
        chi2_gof = fit.chi_square(constants, fit_params, rate_only=rate_only,
                avg_near=avg_near, near_ads=near_ads, variant='pearson')
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
    finally:
        # Delete the temporary fit config file
        os.remove(fit_file_name)
    return (label, index, sin2, best_sin2, sin2_error, dm2ee, dm2ee_best, dm2ee_error,
            chi2_min, chi2_gof, rate_only, avg_near)


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
)
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('--label')
    parser.add_argument('--toy-output')
    parser.add_argument('--toy-scripts')
    parser.add_argument('--toy-config')
    parser.add_argument('--find-errors', action='store_true')
    parser.add_argument('--source-category', type=int,
            help='1: no fluctuations, 2: full fluctuations, 3: far stat only, 4: full stat'
            ', 5: near stat only, 6: rel escale, 7: core spectra, 8: reactor power'
            ', 9: no flucts w accs, 10: acc flucts, 11: no flucts w li9'
            ', 12: li9 flucts, 13: no flucts w all bg, 14: all bg flucts'
            ', 15: detection eff, 16: no flucts w fastn, 17: fastn flucts'
            ', 18: no flucts w amc, 19: amc flucts, 20: resolution'
            ', 21: scintillator nonlinearity, 22: all fluctuations (syst & stat)'
            ', 23: all fluctuations except for solar osc, 24: all syst no stat flucts')
    parser.add_argument('--fit-config')
    parser.add_argument('--gen-mc-only', action='store_true')
    parser.add_argument('--dump-mc', action='store_true')
    parser.add_argument('--test-mode', action='store_true')
    parser.add_argument('--rate-only', action='store_true')
    parser.add_argument('--avg-near', action='store_true')
    parser.add_argument(
        "--pulls", nargs='*', choices=pull_choices+('all',), default=[]
    )
    parser.add_argument('--stage', choices=('6ad', '7ad', '8ad'), default='8ad')
    args = parser.parse_args()
    source_categories = {
        1: 'no fluctuations default nGd binning',
        2: 'full fluctuations default nGd binning',
        3: 'far only fluctuations default nGd binning',
        4: 'near+far stat fluctuations default nGd binning',
        5: 'near only fluctuations default nGd binning',
        6: 'relative energy scale default nGd binning',
        7: 'core spectra default nGd binning',
        8: 'reactor power default nGd binning',
        9: 'no fluctuations w acc bg default nGd binning',
        10: 'acc bg fluctuations default nGd binning',
        11: 'no fluctuations w li9 bg default nGd binning',
        12: 'li9 bg fluctuations default nGd binning',
        13: 'no fluctuations w all all bg default nGd binning',
        14: 'all bg fluctuations default nGd binning',
        15: 'efficiency default nGd binning',
        16: 'no fluctuations w fastn bg default nGd binning',
        17: 'fastn bg fluctuations default nGd binning',
        18: 'no fluctuations w amc bg default nGd binning',
        19: 'amc bg fluctuations default nGd binning',
        20: 'resolution default nGd binning',
        21: 'scint nonlinearity default nGd binning',
        22: 'all sys and stat default nGd binning',
        23: 'all sys and stat except solar osc default nGd binning',
        24: 'all syst no stat default nGd binning',
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
        args.gen_mc_only,
        args.dump_mc,
        args.rate_only,
        args.avg_near,
        args.pulls,
        args.stage,
    )
