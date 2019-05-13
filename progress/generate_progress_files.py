'''
Generate simple files.

'''
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('total', type=int)
parser.add_argument('start', type=int)
parser.add_argument('perfile', type=int)
parser.add_argument('--start-index', type=int, default=0)
parser.add_argument('--suffix', default='')
args = parser.parse_args()

start = args.start
perfile = args.perfile
nfiles = int(args.total/perfile)
start_index = args.start_index

for fileno in range(nfiles):
    filename = '%d.prog%s' % ((fileno + start_index), args.suffix)
    startfile = start + fileno*perfile
    endfile = startfile + perfile - 1
    with open(filename, 'w') as f:
        f.write('%d\n%d\n%d' % (startfile, endfile, startfile))


