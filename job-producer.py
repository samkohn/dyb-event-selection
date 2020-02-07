import argparse
import re
import os.path
import subprocess

def get_run_info(line_number, filename):
    command_result = subprocess.run(['sed', '%dq;d' % line_number,
        filename], capture_output=True)
    line_output = command_result.stdout.decode()
    result = re.match('^(\d+)\s+(\d+)\s+(\d).*$', line_output)
    run = int(result.group(1))
    fileno = int(result.group(2))
    site = int(result.group(3))
    return (run, fileno, site)

def get_file_location(run, fileno):
    find_file=os.path.expanduser('~mkramer/projscratch/p17b/code/p17b_find/p17b_find')
    command_result = subprocess.run([find_file, str(run), str(fileno)],
            capture_output=True)
    location = command_result.stdout.decode()
    return location.strip()

def generate_worker_command(line, runlist, command_template):
    run, fileno, site = get_run_info(line, runlist)
    infile = get_file_location(run, fileno)
    command = command_template.format(filepath=infile, run=run, fileno=fileno,
            site=site)
    return command

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
args = parser.parse_args()
runlist = args.runlist
ntasks = args.ntasks
start_task = args.start_task
command_template = args.command
lines_to_read = range(start_task, ntasks+start_task)
for line in lines_to_read:
    command = generate_worker_command(line, runlist, command_template)
    print(command)
