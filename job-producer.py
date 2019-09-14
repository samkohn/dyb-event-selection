import argparse
import re
import os.path
import subprocess
import zmq

parser = argparse.ArgumentParser()
parser.add_argument('runlist')
parser.add_argument('--task', type=int)
parser.add_argument('--npertask', type=int)
parser.add_argument('--nworkers', type=int)
parser.add_argument('--jobscript')
args = parser.parse_args()
runlist = args.runlist
task_number = args.task
npertask = args.npertask
jobscript = args.jobscript
nworkers = args.nworkers

lines_to_read = range(task_number*npertask+1, (task_number+1)*npertask+1)
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

context = zmq.Context()
producer = context.socket(zmq.REP)
producer.bind('tcp://*:52837')
print('bound')

for line in lines_to_read:
    producer.recv()
    run, fileno, site = get_run_info(line, runlist)
    infile = get_file_location(run, fileno)
    command = 'python {script} -i {filepath} -n {nevents} --site {site}'.format(
            filepath=infile, nevents=-1, site=site, script=jobscript)
    producer.send(command)
    print('sent job command')

for _ in range(nworkers):
    producer.recv()
    producer.send('DONE')
