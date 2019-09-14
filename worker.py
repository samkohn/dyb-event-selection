import zmq
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--outdir')
args = parser.parse_args()
outdir = args.outdir

os.chdir(outdir)

context = zmq.Context()
worker = context.socket(zmq.REQ)
worker.setsockopt(zmq.RCVTIMEO, 5000)
worker.connect('tcp://localhost:52837')
print('connected')

ntasks_recv = 0
while True:
    worker.send('GIMME')
    print('sent command')
    line = worker.recv()
    print('received command')
    if line == 'DONE':
        break
    ntasks_recv += 1
    print(line)
    os.system(line)

print('received %d tasks' % ntasks_recv)
