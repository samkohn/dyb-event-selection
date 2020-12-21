import argparse
import re
import os.path
import subprocess

from common import dets_for

def get_run_info(line):
    result = re.match('^(\d+)\s+(\d+)\s+(\d).*$', line)
    run = int(result.group(1))
    fileno = int(result.group(2))
    site = int(result.group(3))
    return (run, fileno, site)

def get_run_info_run_only(line):
    result = re.match('^(\d+)\s+(\d+).*$', line)
    run = int(result.group(1))
    site = int(result.group(2))
    return (run, site)

def get_file_location(run, fileno):
    find_file=os.path.expanduser('~mkramer/projscratch/p17b/code/p17b_find/p17b_find')
    command_result = subprocess.run([find_file, str(run), str(fileno)],
            capture_output=True)
    location = command_result.stdout.decode()
    return location.strip()

def generate_worker_command(line, command_template, run_only, per_ad):
    if run_only:
        run, site = get_run_info_run_only(line)
        args = {'run': run, 'site': site}
    else:
        run, fileno, site = get_run_info(line)
        if 'filepath' in command_template:
            infile = get_file_location(run, fileno)
            args = {'run': run, 'fileno': fileno, 'site': site, 'filepath': infile}
        else:
            args = {'run': run, 'fileno': fileno, 'site': site}
    if per_ad:
        ads = dets_for(site, run)
        commands = []
        for ad in ads:
            args['ad'] = ad
            commands.append(command_template.format(**args))
        return commands
    else:
        command = command_template.format(**args)
        return command

if __name__ == '__main__':
    parser = argparse.ArgumentParser(epilog=
            '''Command template variables:
      - {filepath} = path to input file
      - {run} = run number, as an int
      - {fileno} = sub-run / file number, as an int
      - {site} = site / hall number (1, 2, or 3)

    The command template will be passed to python's string.format function
    so the template can contain more detailed formatting instructions (such
    as zero-padding the fileno using {fileno:>04}).''',
    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('runlist')
    parser.add_argument('ntasks', type=int)
    parser.add_argument('--start-task', type=int, default=1)
    parser.add_argument('command', metavar='COMMAND_TEMPLATE')
    parser.add_argument('--run-only', action='store_true')
    parser.add_argument('--per-ad', action='store_true')
    args = parser.parse_args()
    runlist = args.runlist
    ntasks = args.ntasks
    start_task = args.start_task
    command_template = args.command
    with open(runlist, 'r') as fin:
        for i, line in enumerate(fin):
            if i < start_task - 1:
                continue
            if i >= ntasks + start_task - 1:
                break
            command = generate_worker_command(line, command_template,
                    args.run_only, args.per_ad)
            if args.per_ad:
                for c in command:
                    print(c)
            else:
                print(command)
