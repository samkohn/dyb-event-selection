import argparse
import re
import os.path
import subprocess
import zmq

def get_run_info(line_number, filename):
    line_output = subprocess.check_output(['sed', '%dq;d' % line_number,
        filename])
    result = re.match('^(\d+)\s+(\d+)\s+(\d).*$', line_output)
    run = int(result.group(1))
    fileno = int(result.group(2))
    site = int(result.group(3))
    return (run, fileno, site)

def get_file_location(run, fileno):
    find_file=os.path.expanduser('~mkramer/projscratch/p17b/code/p17b_find/p17b_find')
    location = subprocess.check_output([find_file, str(run), str(fileno)])
    return location.strip()

def generate_worker_command(line, runlist, prefix='python '):
    run, fileno, site = get_run_info(line, runlist)
    infile = get_file_location(run, fileno)
    command = prefix + '{script} -i {filepath} -n {nevents} --site {site}'.format(
            filepath=infile, nevents=-1, site=site, script=jobscript)
    return command

parser = argparse.ArgumentParser()
parser.add_argument('runlist')
parser.add_argument('ntasks', type=int)
parser.add_argument('--start-task', type=int, default=1)
parser.add_argument('--jobscript')
parser.add_argument('-o', '--output')
args = parser.parse_args()
runlist = args.runlist
ntasks = args.ntasks
start_task = args.start_task
jobscript = args.jobscript
outfile = args.output
lines_to_read = range(start_task, ntasks+start_task)
with open(outfile, 'w') as f:
    for line in lines_to_read:
        command = generate_worker_command(line, runlist,
                prefix='job.sh ')
        f.write(command + '\n')
