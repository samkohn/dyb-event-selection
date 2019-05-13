'''
Get the next job to run.

'''
from __future__ import print_function
import re
import sys
import os
import subprocess
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument('--sublist', type=int)
parser.add_argument('--srcdir')
parser.add_argument('--outdir')
args = parser.parse_args()

filename = os.path.join(args.srcdir, 'progress/%d.prog' % args.sublist)
lockname = os.path.join(args.srcdir, 'progress/%d.lock' % args.sublist)
while os.access(lockname, os.F_OK):
    time.sleep(0.1)
subprocess.check_output(['touch', lockname])
with open(filename, 'r') as f:
    first_index = int(f.readline().strip())
    max_index = int(f.readline().strip())
    current_index = int(f.readline().strip())
    next_index = current_index + 1
if current_index > max_index:
    try:
        subprocess.check_output(['rm', lockname])
    except:
        subprocess.check_output(['touch', ('progress/__possible_error_%d__' %
            current_index)])
    finally:
        print('echo done')
        sys.exit(0)

with open(filename, 'w') as f:
    f.write('%d\n%d\n%d\n' % (first_index, max_index, next_index))
try:
    subprocess.check_output(['rm', lockname])
except:
    subprocess.check_output(['touch', ('progress/__possible_error_%d__' %
        current_index)])

output = subprocess.check_output(['sed', '%dq;d' % current_index,
    os.path.join(args.srcdir, 'run_list.txt')])

result = re.match('^(\d+)\s+(\d+)\s+(\d).*$', output)
run = int(result.group(1))
fileno = int(result.group(2))
site = int(result.group(3))

args = {
        'srcdir': args.srcdir,
        'outdir': args.outdir,
        'nevents': -1,
        'run': run,
        'fileno': fileno,
        'site': site,
        'sublist': args.sublist,
        'lineno': current_index,
        }



print('%(srcdir)s/job.sh %(srcdir)s %(outdir)s %(run)d %(fileno)d %(nevents)d '
        '%(sublist)d %(site)d %(lineno)d' % args)
