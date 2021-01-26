"""Process an entire run from NuWa output through to num coincidences."""

import argparse
import functools
import glob
import logging
import multiprocessing
import os
import sqlite3
import subprocess
import sys
import time

import tenacity

import aggregate_stats
import common
import compute_acc_total_efficiency
import compute_num_coincidences
import compute_singles
import create_accidentals
import create_accidentals_adtime
import create_singles_adtime
import extract_flasher_vars
import first_pass_adtime
import fit_delayed_adtime
import job_producer
import process_adtime
import subtract_accidentals


NUM_MULTIPROCESSING = 64

NOMINAL_LABEL = "reject_resid_flash_nominal"
ADTIME_LABEL = "reject_resid_flash_adtime"

def time_execution(func):
    """Convention that first arg is always the run number"""
    @functools.wraps(func)
    def inner_func(*args, **kwargs):
        logging.info('Starting %s, Run %d', func.__name__, args[0])
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        logging.info(
            'Ran %s (Run %d) in %.2f sec',
            func.__name__,
            args[0],
            end_time - start_time,
        )
        return result
    return inner_func

def setup_database(database, progress_db):
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

            CREATE TABLE delayed_energy_fits (
                Hall INTEGER NOT NULL,
                DetNo INTEGER NOT NULL,
                Source TEXT NOT NULL,
                Peak REAL,
                Peak_error REAL,
                Resolution REAL,
                Resolution_error REAL,
                ExpScale REAL,
                ExpScale_error REAL,
                PeakFraction REAL,
                PeakFraction_error REAL,
                Normalization REAL,
                Normalization_error REAL,
                ChiSquare REAL,
                NumBins INTEGER,
                NumFitParams INTEGER,
                PRIMARY KEY(Hall, DetNo, Source)
            );

            CREATE TABLE distance_time_cut_efficiency (
                Hall INTEGER NOT NULL,
                DetNo INTEGER NOT NULL,
                Source TEXT NOT NULL,
                Efficiency REAL,
                StatError REAL,
                BinningError REAL,
                PRIMARY KEY(Hall, DetNo, Source)
            );

            CREATE TABLE accidentals_spectrum (
                Label TEXT NOT NULL,
                Hall INTEGER NOT NULL,
                DetNo INTEGER NOT NULL,
                BinningId INTEGER NOT NULL,
                BinIndex INTEGER NOT NULL,
                Spectrum REAL NOT NULL,
                PRIMARY KEY(Label, Hall, DetNo, BinningId, BinIndex)
            );
        '''
        )
    with common.get_db(progress_db) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            CREATE TABLE processing_progress (
                RunNo INTEGER NOT NULL,
                Hall INTEGER NOT NULL,
                DetNo INTEGER NOT NULL,
                FirstPass INTEGER,
                Process INTEGER,
                Hadd INTEGER,
                AggregateStats INTEGER,
                CreateSingles INTEGER,
                ComputeSingles INTEGER,
                CreateAccidentals INTEGER,
                SubtractAccidentals INTEGER,
                AccEfficiency INTEGER,
                NumCoincidences INTEGER,
                PRIMARY KEY(RunNo, DetNo)
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
        'acc_nominal_ad{det}',
        'sub_nominal_ad{det}',
        'acc_using_adtime_ad{det}',
        'sub_using_adtime_ad{det}',
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


@tenacity.retry(
    reraise=True,
    wait=tenacity.wait_random_exponential(max=60),
    retry=tenacity.retry_if_exception_type(sqlite3.Error),
)
def _update_progress_db(database, run, site, ad, column):
    """Update the progress tracker by checking off the given column.

    If ad is None then use all ADs active during that run.
    """
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        if isinstance(column, str):
            # Need to use unsafe string interpolation for column name
            update_string = f'''
                INSERT INTO
                    processing_progress (RunNo, Hall, DetNo, {column})
                VALUES (?, ?, ?, 1)
                ON CONFLICT
                    (RunNo, DetNo)
                DO
                UPDATE
                SET
                    {column} = 1;
            '''
        else:  # other iterable:
            column_list_commas = ', '.join(column)
            column_set_equal_1 = ', '.join([f'{c} = 1' for c in column])
            ones = ', '.join(['1'] * len(column))
            update_string = f'''
                INSERT INTO
                    processing_progress (RunNo, Hall, DetNo, {column_list_commas})
                VALUES (?, ?, ?, {ones})
                ON CONFLICT
                    (RunNo, DetNo)
                DO
                UPDATE
                SET
                    {column_set_equal_1};
            '''
        if ad is None:
            ads = common.dets_for(site, run)
            rows = [(run, site, det) for det in ads]
            cursor.executemany(update_string, rows)
        else:
            cursor.execute(update_string, (run, site, ad))
    if ad is None:
        logging.debug('Updated progress db for Run %d, all ADs, script %s', run, column)
    else:
        logging.debug('Updated progress db for Run %d, AD %d, script %s', run, ad, column)


def _fetch_progress(database):
    """Fetch a dict of progress from the progress tracker,
    keyed by run then script name.
    """
    with common.get_db(database) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('''SELECT * FROM processing_progress''')
        rows = cursor.fetchall()
    progress = {}
    for row in rows:
        progress[row['RunNo'], row['DetNo']] = dict(row)
    return progress


def _should_run(run, site, name, progress):
    """Return True if all ADs have finished the given processing."""
    ads = common.dets_for(site, run)
    should_run = False
    for ad in ads:
        run_progress = progress.get((run, ad), None)
        if run_progress is None:
            should_run = True
        elif run_progress[name] == 0 or run_progress[name] is None:
            should_run = True
    return should_run


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
            if int(parts[0]) >= start_run:
                found_first_run = True
            elif not found_first_run:
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
        else:  # no break, means that the file ended and we didn't save the last run
            # Save the last run
            if current_run is not None:
                result[current_run] = (current_site, current_filenos)
    return result


def _apply_first_pass(run, site, fileno, output_location):
    """Multiprocessing-friendly call of first_pass."""
    logging.debug('[first_pass] Running on Run %d, file %d', run, fileno)
    num_events = -1
    debug = False
    infile_location = job_producer.get_file_location(run, fileno)
    if first_pass_adtime.is_complete(
            run,
            site,
            fileno,
            infile_location,
            output_location
    ):
        logging.debug('[first_pass] Already completed Run %d, file %d', run, fileno)
    else:
        first_pass_adtime.main(
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
    subdir = str(run)[0:3] + '00'  # e.g. 21221 -> 21200
    output_location = os.path.join(raw_output_path, f'EH{site}', subdir)
    func_inputs = [(run, site, fileno, output_location) for fileno in filenos]
    with multiprocessing.Pool(NUM_MULTIPROCESSING) as pool:
        pool.starmap(_apply_first_pass, func_inputs)
    return


def _apply_extract_flashers(run, site, fileno, raw_output_path):
    logging.debug('[extract_flashers] Running on Run %d, file %d', run, fileno)
    subdir = str(run)[0:3] + '00'  # e.g. 21221 -> 21200
    output_location = os.path.join(raw_output_path, f'EH{site}', subdir, 'flashers')
    os.makedirs(output_location, exist_ok=True)
    ads = common.dets_for(site, run)
    events_files = [
        os.path.join(
            raw_output_path,
            f'EH{site}',
            subdir,
            f'events_ad{ad}_{run}_{fileno:>04}.root'
        ) for ad in ads
    ]
    num_events = -1
    debug = False
    infile_location = job_producer.get_file_location(run, fileno)
    extract_flasher_vars.main(
        num_events,
        infile_location,
        output_location,
        (run, fileno, site),
        events_files,
        debug,
    )
    logging.debug('[extract_flashers] Finished Run %d, file %d', run, fileno)
    return

@time_execution
def run_extract_flashers(run, site, filenos, raw_output_path):
    func_inputs = [(run, site, fileno, raw_output_path) for fileno in filenos]
    with multiprocessing.Pool(NUM_MULTIPROCESSING) as pool:
        pool.starmap(_apply_extract_flashers, func_inputs)
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
        if process_adtime.is_complete(events_location, output_location):
            logging.debug('[process] Already completed Run %d, file %d', run, fileno)
        else:
            try:
                process_adtime.main(
                    num_events,
                    events_location,
                    muons_location,
                    output_location,
                    (run, fileno),
                    ad,
                    debug,
                )
            except:
                logging.exception('Error in Run %d, file %d, ad %d, input muons location %s', run,
                    fileno, ad, muons_location
                )
                raise
    logging.debug('[process] Finished Run %d, file %d', run, fileno)
    return

@time_execution
def run_process(run, site, filenos, raw_output_path, processed_output_path):
    """Run the coincidence search and muon vetos."""
    subdir = str(run)[0:3] + '00'
    input_prefix = os.path.join(raw_output_path, f'EH{site}', subdir)
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
    import ROOT
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
        # Check to see if the file already exists
        if os.path.isfile(outfile):
            outfile_opened = ROOT.TFile(outfile, 'READ')
            try:
                if not outfile_opened.IsZombie():
                    # Check to see if the outfile has the right number of events
                    outfile_events = outfile_opened.Get('ad_events').GetEntries()
                    outfile_opened.Close()
                    infile_events = 0
                    for infile in infiles:
                        infile_opened = ROOT.TFile(infile, 'READ')
                        infile_events += infile_opened.Get('ad_events').GetEntries()
                        infile_opened.Close()
                    if infile_events == outfile_events:
                        logging.debug('[hadd] Found existing file. Skipping. %s', outfile)
                        continue
            except ReferenceError:  # If file doesn't have correct TTree
                pass
            finally:
                outfile_opened.Close()
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
        if aggregate_stats.is_complete(run, ad, outfile, database):
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
        if create_singles_adtime.is_complete(infile, outfile):
            logging.debug('[create_singles] Found existing file. Skipping. %s', outfile)
        else:
            create_singles_adtime.main(infile, outfile, ttree_name)
    return


@time_execution
@tenacity.retry(
    reraise=True,
    wait=tenacity.wait_random_exponential(max=60),
    retry=tenacity.retry_if_exception_type(sqlite3.Error),
)
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
        # Ideally should check to see if rate has been computed before.
        # But, this is so fast that I will just re-compute every time.
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
    update_db = None
    for ad in ads:
        infile = os.path.join(
            path_prefix,
            f'singles_ad{ad}/singles_ad{ad}_{run}.root',
        )
        outfile = os.path.join(
            path_prefix,
            f'acc_nominal_ad{ad}/acc_ad{ad}_{run}.root',
        )
        if create_accidentals.is_complete(infile, outfile):
            logging.debug(
                '[create_accidentals] Found existing file. Skipping. %s', outfile
            )
        else:
            logging.debug(
                '[create_accidentals] Running file. %s', outfile
            )
            pairing_note = NOMINAL_LABEL
            create_accidentals.main(
                infile,
                outfile,
                ttree_name,
                pairing,
                pairing_note,
                update_db,
            )
        outfile = os.path.join(
            path_prefix,
            f'acc_using_adtime_ad{ad}/acc_ad{ad}_{run}.root',
        )
        if create_accidentals_adtime.is_complete(infile, outfile):
            logging.debug(
                '[create_accidentals] Found existing file. Skipping. %s', outfile
            )
        else:
            logging.debug(
                '[create_accidentals] Running file. %s', outfile
            )
            pairing_note = ADTIME_LABEL
            create_accidentals_adtime.main(
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
    update_db = True
    for ad in ads:
        datafile = os.path.join(
            path_prefix,
            f'hadded_ad{ad}/out_ad{ad}_{run}.root',
        )
        accfile = os.path.join(
            path_prefix,
            f'acc_nominal_ad{ad}/acc_ad{ad}_{run}.root',
        )
        outfile = os.path.join(
            path_prefix,
            f'sub_nominal_ad{ad}/sub_ad{ad}_{run}.root',
        )
        if os.path.isfile(outfile):
            logging.debug(
                '[subtract_accidentals] Found existing file. Skipping. %s', outfile
            )
        else:
            label = NOMINAL_LABEL
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
        outfile = os.path.join(
            path_prefix,
            f'sub_using_adtime_ad{ad}/sub_ad{ad}_{run}.root',
        )
        accfile = os.path.join(
            path_prefix,
            f'acc_using_adtime_ad{ad}/acc_ad{ad}_{run}.root',
        )
        if os.path.isfile(outfile):
            logging.debug(
                '[subtract_accidentals] Found existing file. Skipping. %s', outfile
            )
        else:
            label = ADTIME_LABEL
            subtract_accidentals.main(
                outfile,
                datafile,
                accfile,
                database,
                ad,
                override_acc_rate,
                label,  # Label containing 'adtime' will use AdTime positions
                update_db,
            )
    return

def run_hadd_sub_files(site, processed_output_path):
    """Combine all accidentals-subtracted files for a given AD."""
    path_prefix = os.path.join(processed_output_path, f'EH{site}')
    input_template = os.path.join(path_prefix, 'sub_{label}_ad{ad}/sub_ad{ad}_?????.root')
    output_template = os.path.join(path_prefix, 'sub_{label}_ad{ad}/sub_ad{ad}.root')
    ads = common.dets_for(site)
    num_cores = str(min(os.cpu_count(), 10))
    for ad in ads:
        for label in ['nominal', 'using_adtime']:
            input_glob = input_template.format(ad=ad, label=label)
            infiles = glob.glob(input_glob)
            outfile = output_template.format(ad=ad, label=label)
            command = ['hadd', '-f', '-v', '0', '-j', num_cores, outfile] + infiles
            subprocess.run(command, check=True)
    return

def run_compute_num_coincidences(processed_output_path, database):
    """Count the number of coincidences for each run.

    Must be run after the delayed energy fits.
    """
    infile_template = os.path.join(
        processed_output_path,
        'EH{site}',
        'hadded_ad{ad}',
        'out_ad{ad}_{run}.root',
    )
    compute_num_coincidences.main(database, database, infile_template, True)
    return


def run_compute_total_acc_efficiency(processed_output_path, database):
    """Compute the accidental sample efficiency."""
    logging.info('beginning nominal acc efficiency')
    infile_template = os.path.join(
        processed_output_path,
        'EH{site}',
        'acc_nominal_ad{ad}',
        'acc_ad{ad}_{run}.root',
    )
    compute_acc_total_efficiency.main(
        database,
        database,
        infile_template,
        NOMINAL_LABEL,
        True,
    )
    logging.info('beginning adtime acc efficiency')
    infile_template = os.path.join(
        processed_output_path,
        'EH{site}',
        'acc_using_adtime_ad{ad}',
        'acc_ad{ad}_{run}.root',
    )
    compute_acc_total_efficiency.main(
        database,
        database,
        infile_template,
        ADTIME_LABEL,
        True,
    )
    return


def run_fit_delayed(site, ad, processed_output_path, database):
    """Fit the delayed energy spectrum and extract the parameters."""
    infile = os.path.join(
        processed_output_path,
        f'EH{site}',
        f'sub_nominal_ad{ad}',
        f'sub_ad{ad}.root',
    )
    outfile = f'fit_nominal_eh{site}_ad{ad}.pdf'
    fit_delayed_adtime.main(infile, outfile, site, ad, NOMINAL_LABEL, database)
    infile = os.path.join(
        processed_output_path,
        f'EH{site}',
        f'sub_using_adtime_ad{ad}',
        f'sub_ad{ad}.root',
    )
    outfile = f'fit_adtime_eh{site}_ad{ad}.pdf'
    fit_delayed_adtime.main(infile, outfile, site, ad, ADTIME_LABEL, database)
    return


def post_processing(processed_output_path, database):
    """Run the scripts that rely on all data being processed through many_runs."""
    sites = (1, 2, 3)
    for site in sites:
        dets = common.dets_for(site)
        run_hadd_sub_files(site, processed_output_path)
        for det in dets:
            run_fit_delayed(site, det, processed_output_path, database)
    run_compute_num_coincidences(processed_output_path, database)
    run_compute_total_acc_efficiency(processed_output_path, database)
    return


def _tasks_for_whole_run(
    run,
    site,
    filenos,
    processed_output_path,
    database,
    progress_db,
    stop_time,
    progress,
):
    """Execute all the tasks that involve the entire run, i.e. not first_pass or
    process."""
    finished = []
    try:
        if _should_run(run, site, 'Hadd', progress):
            run_hadd(run, site, filenos, processed_output_path)
            finished.append('Hadd')
        else:
            logging.debug('[compute_singles] Skipping Run %d based on db progress', run)
        if _should_run(run, site, 'AggregateStats', progress):
            run_aggregate_stats(run, site, filenos, processed_output_path, database)
            finished.append('AggregateStats')
        else:
            logging.debug('[aggregate_stats] Skipping Run %d based on db progress', run)
        if time.time() + 20 * 60 > stop_time:
            return
        if _should_run(run, site, 'CreateSingles', progress):
            run_create_singles(run, site, processed_output_path)
            finished.append('CreateSingles')
        else:
            logging.debug('[create_singles] Skipping Run %d based on db progress', run)
        if time.time() + 20 * 60 > stop_time:
            return
        if _should_run(run, site, 'ComputeSingles', progress):
            run_compute_singles(run, site, processed_output_path, database)
            finished.append('ComputeSingles')
        else:
            logging.debug('[compute_singles] Skipping Run %d based on db progress', run)
        if _should_run(run, site, 'CreateAccidentals', progress):
            run_create_accidentals(run, site, processed_output_path)
            finished.append('CreateAccidentals')
        else:
            logging.debug(
                '[create_accidentals] Skipping Run %d based on db progress', run
            )
        if time.time() > stop_time:
            return
        if _should_run(run, site, 'SubtractAccidentals', progress):
            run_subtract_accidentals(run, site, processed_output_path, database)
            finished.append('SubtractAccidentals')
        else:
            logging.debug(
                '[subtract_accidentals] Skipping Run %d based on db progress', run
            )
        logging.info('Finished full processing for Run %d', run)
    except Exception:
        logging.exception('Exception in Run %d', run)
    finally:
        if len(finished) > 0:
            _update_progress_db(progress_db, run, site, None, finished)
    return


def main(
    start_run,
    end_run,
    run_list_file,
    raw_output_path,
    processed_output_path,
    database,
    progress_db,
    max_runtime_sec,
):
    if max_runtime_sec == -1:
        stop_time = time.time() + 100 * 24 * 60 * 60  # forever
    else:
        stop_time = time.time() + max_runtime_sec
    run_info = get_site_filenos_for_run_range(start_run, end_run, run_list_file)
    progress = _fetch_progress(progress_db)
    logging.info('Prepping for %d runs', len(run_info))
    logging.info('First few runs: %s', str(list(run_info.keys())[:5]))
    # Execute each run's first_pass and process in series since they already use
    # multiprocessing pools.
    for run, (site, filenos) in run_info.items():
        if time.time() > stop_time:
            return
        if _should_run(run, site, 'FirstPass', progress):
            run_first_pass(run, site, filenos, raw_output_path)
            _update_progress_db(progress_db, run, site, None, 'FirstPass')
        else:
            logging.debug('[first_pass] Skipping Run %d based on db progress', run)
        if _should_run(run, site, 'ExtractFlashers', progress):
            run_extract_flashers(run, site, filenos, raw_output_path)
            _update_progress_db(progress_db, run, site, None, 'ExtractFlashers')
        else:
            logging.debug('[extract_flashers] Skipping Run %d based on db progress', run)
        if _should_run(run, site, 'Process', progress):
            run_process(run, site, filenos, raw_output_path, processed_output_path)
            _update_progress_db(progress_db, run, site, None, 'Process')
        else:
            logging.debug('[process] Skipping Run %d based on db progress', run)
    if time.time() > stop_time:
        return
    # Execute all runs' remaining tasks in parallel to take advantage of many cores.
    with multiprocessing.Pool(NUM_MULTIPROCESSING, maxtasksperchild=1) as pool:
        pool.starmap(
            _tasks_for_whole_run,
            [
                (
                    run,
                    site,
                    filenos,
                    processed_output_path,
                    database,
                    progress_db,
                    stop_time,
                    progress
                )
                for run, (site, filenos) in run_info.items()
            ],
            chunksize=2,  # Prioritize load balancing over optimizing overhead
        )
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
    parser.add_argument('--progress-db',
        help='database storing progress. leave blank to use same as --database'
    )
    parser.add_argument('--raw-output',
        help='base directory for output of raw "events" and "muons" files',
    )
    parser.add_argument('--processed-output',
        help='base directory for output of processed files',
    )
    parser.add_argument('--max-runtime-sec', type=int, default=-1)
    args = parser.parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format='%(name)s - %(levelname)s - %(asctime)s - %(message)s',
        stream=sys.stdout,
    )
    if args.progress_db is None:
         progress_db = args.database
     else:
         progress_db = args.progress_db
    if args.setup_directories:
        setup_directory_structure(args.raw_output, args.processed_output)
    if args.database is not None and not os.path.isfile(args.database):
        setup_database(args.database, progress_db)
    if args.end_run is not None:
        main(
            args.run,
            args.end_run,
            args.run_list,
            args.raw_output,
            args.processed_output,
            args.database,
            progress_db,
            args.max_runtime_sec,
        )
    else:
        main(
            args.run,
            args.run,
            args.run_list,
            args.raw_output,
            args.processed_output,
            args.database,
            progress_db,
            args.max_runtime_sec,
        )
