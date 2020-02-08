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
worker.setsockopt(zmq.RCVTIMEO, 25000)
worker.connect('tcp://0.0.0.0:52837')
print('connected')

ntasks_recv = 0
task_id = 0
while ntasks_recv < args.max_tasks:
    worker.send(str(task_id).encode())
    print('sent command')
    task_id, line = worker.recv_multipart()
    task_id = int(task_id)
    print('received command')
    if line == b'DONE':
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

worker.send(str(task_id).encode() + b'DONE')
worker.recv()

print('received %d tasks' % ntasks_recv)
