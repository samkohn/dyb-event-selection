'''
Get the next job to run.

'''
from __future__ import print_function
import re
import sys
import subprocess

with open('progress/0.prog', 'r') as f:
    current_index = int(f.readline().strip())
    next_index = current_index + 1
    max_index = int(f.readline().strip())
if current_index > max_index:
    print('echo done')
    sys.exit(0)

with open('progress/0.prog', 'w') as f:
    f.write('%d\n%d\n' % (next_index, max_index))

output = subprocess.check_output(['sed', '%dq;d' % current_index,
    'run_list.txt'])

result = re.match('^(\d+)\s+(\d+)\s+(\d).*$', output)
run = int(result.group(1))
fileno = int(result.group(2))
site = int(result.group(3))

args = {
        'srcdir': '~/dyb-event-selection',
        'outdir': '/global/cscratch1/sd/skohn/dyb5',
        'nevents': -1,
        'run': run,
        'fileno': fileno,
        'site': site,
        }



print('%(srcdir)s/job.sh %(srcdir)s %(outdir)s %(run)d %(fileno)d %(nevents)d &'
        % args)
