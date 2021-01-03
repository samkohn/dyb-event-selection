"""Process an entire run from NuWa output through to num coincidences."""

import argparse
import functools
import logging
import multiprocessing
import os
import subprocess
import time

import aggregate_stats
import common
import compute_singles
import create_accidentals
import create_singles
import first_pass
import job_producer
import process
import subtract_accidentals


NUM_MULTIPROCESSING = 64

def time_execution(func):
    @functools.wraps(func)
    def inner_func(*args, **kwargs):
        logging.info('Starting %s', func.__name__)
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        logging.info(
            'Ran %s in %.2f sec',
            func.__name__,
            end_time - start_time,
        )
        return result
    return inner_func

def setup_database(database):
    """Create the database tables to store the analysis results."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.executescript('''
            CREATE TABLE num_coincidences_by_run (
                RunNo INTEGER,
                DetNo INTEGER,
                NumCoincidences INTEGER,
                Label TEXT,
                PRIMARY KEY(RunNo, DetNo, Label)
            );

            CREATE TABLE muon_rates (
                RunNo INTEGER NOT NULL,
                DetNo INTEGER NOT NULL,
                Count INTEGER,
                Livetime_ns INTEGER,
                Rate_Hz REAL,
                Efficiency REAL,
                PRIMARY KEY(RunNo, DetNo)
            );

            CREATE TABLE singles_rates (
                RunNo INTEGER NOT NULL,
                DetNo INTEGER NOT NULL,
                Iteration INTEGER,
                Rate_Hz REAL,
                Rate_Hz_error REAL,
                IsolatedEventCount INTEGER,
                IsolatedEventRate_Hz REAL,
                CorrelatedRate_Hz REAL,
                MultiplicityVetoEfficiency REAL,
                PRIMARY KEY(RunNo, DetNo)
            );

            CREATE TABLE runs (
                RunNo INTEGER PRIMARY KEY,
                Hall INTEGER,
                Start_time INTEGER
            );

            CREATE TABLE accidental_subtraction (
                RunNo INTEGER NOT NULL,
                DetNo INTEGER NOT NULL,
                Label TEXT NOT NULL,
                BaseRate_Hz REAL,
                DistanceTime_DT_Eff REAL,
                AccScaleFactor REAL,
                DTCrossCheck REAL,
                DTCrossCheck_error REAL,
                DistanceCrossCheck REAL,
                DistanceCrossCheck_error REAL,
                Total_Acc_Eff REAL,
                Total_Acc_Eff_err REAL,
                PRIMARY KEY(RunNo, DetNo, Label)
            );
        '''
        )
    return


def setup_directory_structure(raw_output_path, processed_output_path):
    """Create the directory structure for output files."""
    ads = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (3, 4)]
    processed_dir_templates = [
        'processed_ad{det}',
        'hadded_ad{det}',
        'singles_ad{det}',
        'acc_ad{det}',
        'sub_ad{det}',
    ]
    paths = []
    for hall, det in ads:
        paths.append(os.path.join(raw_output_path, f'EH{hall}'))
        paths.append(os.path.join(processed_output_path, f'EH{hall}'))
        for template in processed_dir_templates:
            paths.append(os.path.join(
                processed_output_path,
                f'EH{hall}',
                template.format(det=det),
            ))
    for path in paths:
        os.makedirs(path, exist_ok=True)


def get_site_filenos_for_run(run, run_list_filename):
    """Loop through the run list file and return a tuple of (site, filenos).

    Assumes the input list is sorted with all files from a given run listed
    consecutively.
    """
    run_str = str(run)
    filenos = []
    with open(run_list_filename, 'r') as f:
        found_run = False
        for line in f:
            parts = line.split()
            if parts[0] == run_str:
                found_run = True
                fileno = int(parts[1])
                site = int(parts[2])  # Only need this once but whatever
                filenos.append(fileno)
            elif found_run:
                break
            else:
                pass
    return site, filenos

def get_site_filenos_for_run_range(start_run, end_run, run_list_filename):
    """Loop through the run list file and return a dict of {run: (site, filenos)}.

    Assumes the input list is sorted with all files from a given run listed
    consecutively.
    """
    result = {}
    first_run_str = str(start_run)
    last_run_str = str(end_run)
    found_first_run = False
    current_run = None
    current_site = None
    current_filenos = None
    with open(run_list_filename, 'r') as f:
        for line in f:
            # search for the first run
            parts = line.split()
            if parts[0] == first_run_str:
                found_first_run = True
            else:
                continue
            if found_first_run:
                this_line_run = int(parts[0])
                if current_run != this_line_run:
                    # End the last run
                    if current_run is not None:
                        result[current_run] = (current_site, current_filenos)
                    # Check to see if we're done
                    if this_line_run > end_run:
                        break
                    # Start a new run
                    current_run = this_line_run
                    current_site = int(parts[2])
                    current_filenos = []
                current_filenos.append(int(parts[1]))
        else:  # nobreak, means that the file ended and we didn't save the last run
            # Save the last run
            if current_run is not None:
                result[current_run] = (current_site, current_filenos)
    return result


def _apply_first_pass(run, fileno, output_location):
    """Multiprocessing-friendly call of first_pass."""
    logging.debug('[first_pass] Running on Run %d, file %d', run, fileno)
    num_events = -1
    debug = False
    infile_location = job_producer.get_file_location(run, fileno)
    sentinel_outfile = os.path.join(output_location, f'muons_{run}_{fileno:>04}.root')
    if os.path.isfile(sentinel_outfile):
        logging.debug('[first_pass] Found existing file. Skipping. %s', sentinel_outfile)
    else:
        first_pass.main(
            num_events,
            infile_location,
            output_location,
            (run, fileno),
            debug,
        )
    logging.debug('[first_pass] Finished Run %d, file %d', run, fileno)
    return

@time_execution
def run_first_pass(run, site, filenos, raw_output_path):
    """Convert from NuWa to slimmed-down ROOT files."""
    output_location = os.path.join(raw_output_path, f'EH{site}')
    func_inputs = [(run, fileno, output_location) for fileno in filenos]
    with multiprocessing.Pool(NUM_MULTIPROCESSING) as pool:
        pool.starmap(_apply_first_pass, func_inputs)
    return


def _apply_process(run, site, fileno, input_prefix, processed_output_path):
    """Multiprocessing-friendly call of process (all ADs within a hall together)."""
    logging.debug('[process] Running on Run %d, file %d', run, fileno)
    muons_location = os.path.join(input_prefix, f'muons_{run}_{fileno:>04}.root')
    num_events = -1
    debug = False
    ads = common.dets_for(site, run)
    for ad in ads:
        logging.debug(
            '[process] Running on Run %d, file %d, EH%d-AD%d', run, fileno, site, ad
        )
        events_location = os.path.join(
            input_prefix, f'events_ad{ad}_{run}_{fileno:>04}.root'
        )
        output_location = os.path.join(
            processed_output_path,
            f'EH{site}',
            f'processed_ad{ad}',
            f'out_ad{ad}_{run}_{fileno:>04}.root',
        )
        if os.path.isfile(output_location):
            logging.debug('[process] Found existing file. Skipping. %s', output_location)
        else:
            process.main(
                num_events,
                events_location,
                muons_location,
                output_location,
                (run, fileno),
                ad,
                debug,
            )
    logging.debug('[process] Finished Run %d, file %d', run, fileno)
    return

@time_execution
def run_process(run, site, filenos, raw_output_path, processed_output_path):
    """Run the coincidence search and muon vetos."""
    input_prefix = os.path.join(raw_output_path, f'EH{site}')
    func_inputs = []
    for fileno in filenos:
        func_inputs.append((
            run, site, fileno, input_prefix, processed_output_path
        ))
    with multiprocessing.Pool(NUM_MULTIPROCESSING) as pool:
        pool.starmap(_apply_process, func_inputs)
    return


@time_execution
def run_hadd(run, site, filenos, processed_output_path):
    """Combine the processed output files into 1 per run."""
    path_prefix = os.path.join(processed_output_path, f'EH{site}')
    input_template = os.path.join(
        path_prefix,
        'processed_ad{ad}/out_ad{ad}_' f'{run}' '_{fileno:>04}.root',
    )
    output_template = os.path.join(
        path_prefix,
        'hadded_ad{ad}/out_ad{ad}_' f'{run}.root',
    )
    ads = common.dets_for(site, run)
    for ad in ads:
        infiles = []
        for fileno in filenos:
            infiles.append(input_template.format(ad=ad, fileno=fileno))
        outfile = output_template.format(ad=ad)
        if os.path.isfile(outfile):
            logging.debug('[hadd] Found existing file. Skipping. %s', outfile)
        else:
            command = ['hadd', '-f', '-v', '0', outfile] + infiles
            subprocess.run(command, check=True)
    return


@time_execution
def run_aggregate_stats(run, site, filenos, processed_output_path, database):
    """Aggregate the muon / livetime statistics for each run."""
    path_prefix = os.path.join(processed_output_path, f'EH{site}')
    input_template = os.path.join(
        path_prefix,
        'processed_ad{ad}/out_ad{ad}_' f'{run}' '_{fileno:>04}.json',
    )
    output_template = os.path.join(
        path_prefix,
        'hadded_ad{ad}/out_ad{ad}_' f'{run}.json',
    )
    ads = common.dets_for(site, run)
    for ad in ads:
        infiles = []
        for fileno in filenos:
            infiles.append(input_template.format(ad=ad, fileno=fileno))
        outfile = output_template.format(ad=ad)
        if os.path.isfile(outfile):
            logging.debug('[aggregate_stats] Found existing file. Skipping. %s', outfile)
        else:
            aggregate_stats.main2(run, infiles, site, ad, outfile, database)
    return


@time_execution
def run_create_singles(run, site, processed_output_path):
    """Create a sample of isolated singles (for acc spectrum & DT, not for rate)."""
    path_prefix = os.path.join(processed_output_path, f'EH{site}')
    ads = common.dets_for(site, run)
    ttree_name = 'ad_events'
    for ad in ads:
        infile = os.path.join(
            path_prefix,
            f'hadded_ad{ad}/out_ad{ad}_{run}.root',
        )
        outfile = os.path.join(
            path_prefix,
            f'singles_ad{ad}/singles_ad{ad}_{run}.root',
        )
        if os.path.isfile(outfile):
            logging.debug('[create_singles] Found existing file. Skipping. %s', outfile)
        else:
            create_singles.main(infile, outfile, ttree_name)
    return


@time_execution
def run_compute_singles(run, site, processed_output_path, database):
    """Compute the singles (underlying uncorrelated) rate."""
    path_prefix = os.path.join(processed_output_path, f'EH{site}')
    ads = common.dets_for(site, run)
    update_db = True
    iteration = 0
    extra_cut = '1'
    for ad in ads:
        infile = os.path.join(
            path_prefix,
            f'hadded_ad{ad}/out_ad{ad}_{run}.root',
        )
        with common.get_db(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT
                    COUNT(*)
                FROM
                    singles_rates
                WHERE
                    RunNo = ?
                    AND DetNo = ?
                ''',
                (run, ad)
            )
            num_existing, = cursor.fetchone()
        if num_existing > 0:
            logging.debug(
                '[compute_singles] Found existing rate. Skipping. Run %d, AD %d', run, ad
            )
        else:
            compute_singles.main(
                infile,
                database,
                update_db,
                iteration,
                extra_cut,
            )
    return


@time_execution
def run_create_accidentals(run, site, processed_output_path):
    """Create a synthetic accidentals sample using the sequential pairing."""
    path_prefix = os.path.join(processed_output_path, f'EH{site}')
    ads = common.dets_for(site, run)
    ttree_name = 'singles'
    pairing = 'sequential'
    pairing_note = ''
    update_db = None
    for ad in ads:
        infile = os.path.join(
            path_prefix,
            f'singles_ad{ad}/singles_ad{ad}_{run}.root',
        )
        outfile = os.path.join(
            path_prefix,
            f'acc_ad{ad}/acc_ad{ad}_{run}.root',
        )
        if os.path.isfile(outfile):
            logging.debug(
                '[create_accidentals] Found existing file. Skipping. %s', outfile
            )
        else:
            create_accidentals.main(
                infile,
                outfile,
                ttree_name,
                pairing,
                pairing_note,
                update_db,
            )
    return


@time_execution
def run_subtract_accidentals(run, site, processed_output_path, database):
    path_prefix = os.path.join(processed_output_path, f'EH{site}')
    ads = common.dets_for(site, run)
    override_acc_rate = None
    label = 'test'
    update_db = True
    for ad in ads:
        datafile = os.path.join(
            path_prefix,
            f'hadded_ad{ad}/out_ad{ad}_{run}.root',
        )
        accfile = os.path.join(
            path_prefix,
            f'acc_ad{ad}/acc_ad{ad}_{run}.root',
        )
        outfile = os.path.join(
            path_prefix,
            f'sub_ad{ad}/sub_ad{ad}_{run}.root',
        )
        if os.path.isfile(outfile):
            logging.debug(
                '[subtract_accidentals] Found existing file. Skipping. %s', outfile
            )
        else:
            subtract_accidentals.main(
                outfile,
                datafile,
                accfile,
                database,
                ad,
                override_acc_rate,
                label,
                update_db,
            )
    return


def _tasks_for_whole_run(run, site, filenos, processed_output_path, database):
    """Execute all the tasks that involve the entire run, i.e. not first_pass or
    process."""
    run_hadd(run, site, filenos, processed_output_path)
    run_aggregate_stats(run, site, filenos, processed_output_path, database)
    run_create_singles(run, site, processed_output_path)
    run_compute_singles(run, site, processed_output_path, database)
    run_create_accidentals(run, site, processed_output_path)
    run_subtract_accidentals(run, site, processed_output_path, database)
    return


def many_runs(
    start_run, end_run, run_list_file, raw_output_path, processed_output_path, database,
):
    run_info = get_site_filenos_for_run_range(start_run, end_run, run_list_file)
    logging.debug('Prepping for %d runs', len(run_info))
    logging.debug('First few runs: %s', str(list(run_info.keys())[:5]))
    # Execute each run's first_pass and process in series since they already use
    # multiprocessing pools.
    for run, (site, filenos) in run_info.items():
        run_first_pass(run, site, filenos, raw_output_path)
        run_process(run, site, filenos, raw_output_path, processed_output_path)
    # Execute all runs' remaining tasks in parallel to take advantage of many cores.
    with multiprocessing.Pool(NUM_MULTIPROCESSING) as pool:
        pool.starmap(
            _tasks_for_whole_run,
            [(run, site, filenos, processed_output_path, database)
                for run, (site, filenos) in run_info.items()
            ],
        )
    return


def main(
    run, run_list_file, raw_output_path, processed_output_path, database,
):
    site, filenos = get_site_filenos_for_run(run, run_list_file)
    filenos = filenos[:63]
    run_first_pass(run, site, filenos, raw_output_path)
    run_process(run, site, filenos, raw_output_path, processed_output_path)
    run_hadd(run, site, filenos, processed_output_path)
    run_aggregate_stats(run, site, filenos, processed_output_path, database)
    run_create_singles(run, site, processed_output_path)
    run_compute_singles(run, site, processed_output_path, database)
    run_create_accidentals(run, site, processed_output_path)
    run_subtract_accidentals(run, site, processed_output_path, database)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', type=int)
    parser.add_argument('--end-run', type=int, default=None)
    parser.add_argument('--run-list')
    parser.add_argument('--setup-directories', action='store_true')
    parser.add_argument('--database',
        help='database to store analysis outputs',
    )
    parser.add_argument('--raw-output',
        help='base directory for output of raw "events" and "muons" files',
    )
    parser.add_argument('--processed-output',
        help='base directory for output of processed files',
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    if args.setup_directories:
        setup_directory_structure(args.raw_output, args.processed_output)
    if args.database is not None and not os.path.isfile(args.database):
        setup_database(args.database)
    if args.end_run is not None:
        asdf
    main(
        args.run,
        args.run_list,
        args.raw_output,
        args.processed_output,
        args.database,
    )
