"""Dynamically pick a good set of runs to process."""

import argparse
import itertools
import subprocess
import sys
import time

import common

def get_unfinished_runs(progress_db):
    """Get the list of runs which have not completed SubtractAccidentals."""
    with common.get_db(progress_db) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                DISTINCT RunNo
            FROM
                processing_progress
            WHERE
                SubtractAccidentals IS NULL
                OR SubtractAccidentals = 0
            ORDER BY RunNo
        '''
        )
        results = [x[0] for x in cursor.fetchall()]
    return results


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks

    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"

    Taken from python itertools recipes.
    """
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('num_nodes', type=int,
        help='Number of independent nodes / versions of this script to launch',
    )
    parser.add_argument('--srun', action='store_true')
    parser.add_argument('--run-list', required=True)
    parser.add_argument('--database', required=True,
        help='database to store analysis outputs',
    )
    parser.add_argument('--progress-db',
        help='database storing progress. leave blank to use same as --database'
    )
    parser.add_argument('--raw-output', required=True,
        help='base directory for output of raw "events" and "muons" files',
    )
    parser.add_argument('--processed-output', required=True,
        help='base directory for output of processed files',
    )
    parser.add_argument('--shifter',
        help='Use the given shifter image in srun',
    )
    parser.add_argument('--max-runtime-sec', type=int, default=-1)
    parser.add_argument('--skip-initial-steps', action='store_true')
    parser.add_argument('-d', '--debug', action='store_true')
    args = parser.parse_args()
    if args.progress_db is None:
        progress_db = args.database
    else:
        progress_db = args.progress_db
    incomplete_runs = get_unfinished_runs(progress_db)
    num_per_node = len(incomplete_runs) // args.num_nodes + 1
    # Must list-ify the iterators for use in multiprocessing
    runs_by_node = [list(x) for x in grouper(incomplete_runs, num_per_node)]
    run_range_by_node = [(x[0], x[-1] if x[-1] is not None else max(y for y in x if y is
        not None)) for x in runs_by_node]
    srun_command = 'srun --no-kill --ntasks=1 -N 1 --wait=0'
    shifter_command = f'shifter --image={args.shifter}' if args.shifter else ''
    for first, last in run_range_by_node:
        command = (
            f'python full_process_adtime.py -r {first} --end-run {last} '
            f' --database {args.database} --progress-db {progress_db} '
            f' --run-list {args.run_list} --raw-output {args.raw_output} '
            f' --processed-output {args.processed_output} '
            f' --max-runtime-sec {args.max_runtime_sec} '
            f' {"-d" if args.debug else ""}'
            f' {"--skip-initial-steps" if args.skip_initial_steps else ""}'
        )
        print(command)

        if args.srun:
            subprocess.run(
                f'{srun_command} {shifter_command} {command} 2>&1 >> '
                f'log_{first}_{last}.txt &',
                shell=True
            )
    time.sleep(args.max_runtime_sec)
    print(f'Just elapsed max_runtime_sec = {args.max_runtime_sec}')
