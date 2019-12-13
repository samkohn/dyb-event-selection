from __future__ import print_function
import zmq
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--outdir')
parser.add_argument('--max-tasks', type=int)
args = parser.parse_args()
outdir = args.outdir

os.chdir(outdir)

context = zmq.Context()
worker = context.socket(zmq.REQ)
worker.setsockopt(zmq.RCVTIMEO, 5000)
worker.connect('tcp://0.0.0.0:52837')
print('connected')

ntasks_recv = 0
task_id = 0
while ntasks_recv < args.max_tasks:
    worker.send(b'{}'.format(task_id))
    print('sent command')
    task_id, line = worker.recv_multipart()
    print('received command')
    if line == 'DONE':
        break
    ntasks_recv += 1
    print(line)
    task_id = int(task_id)
    try:
        result = os.system(line)
        if result != 0:
            task_id *= -1
    except:
        task_id *= -1

worker.send(b'{}DONE'.format(task_id))
worker.recv()

print('received %d tasks' % ntasks_recv)
