import argparse
import re
import os.path
import subprocess

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

def generate_worker_command(line, runlist, nh, prefix='python '):
    run, fileno, site = get_run_info(line, runlist)
    infile = get_file_location(run, fileno)
    command = prefix + '{script} -i {filepath} -n {nevents} --site {site}{nh}'.format(
            filepath=infile, nevents=-1, site=site, script=jobscript, nh=nh)
    return command

def generate_worker_command_basic(line, runlist, selection_name, outdir,
        prefix='python '):
    run, fileno, site = get_run_info(line, runlist)
    infile = get_file_location(run, fileno)
    command = prefix + ('{script} -i {filepath} -o {outdir}'
        ' --selection {name}').format(script=jobscript, filepath=infile,
                nevents=-1, outdir=outdir, name=selection_name)
    return command

parser = argparse.ArgumentParser()
parser.add_argument('runlist')
parser.add_argument('ntasks', type=int)
parser.add_argument('--start-task', type=int, default=1)
parser.add_argument('--jobscript')
group = parser.add_mutually_exclusive_group()
group.add_argument('--nh', action='store_true')
group.add_argument('--selection')
parser.add_argument('-o', '--output')
parser.add_argument('--outdir', default='./')
args = parser.parse_args()
runlist = args.runlist
ntasks = args.ntasks
start_task = args.start_task
jobscript = args.jobscript
outfile = args.output
nh = ' --nh' if args.nh else ''
selection = args.selection
lines_to_read = range(start_task, ntasks+start_task)
with open(outfile, 'w') as f:
    for line in lines_to_read:
        if len(nh) > 0:
            command = generate_worker_command(line, runlist, nh,
                    prefix='job.sh ')
        else:
            command = generate_worker_command_basic(line, runlist, selection,
                    args.outdir, prefix='job.sh ')
        f.write(command + '\n')
