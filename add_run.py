import argparse
import json

import common

parser = argparse.ArgumentParser()
parser.add_argument('infile')
args = parser.parse_args()

with open(args.infile, 'r') as f:
    data = json.load(f)
run = data['run']
hall = data['site']
start_time = int(data['start_time']/1e9)

with common.get_db('parameters.db') as conn:
    c = conn.cursor()
    c.execute('INSERT OR REPLACE INTO runs VALUES (?, ?, ?)', (run, hall,
        start_time))
