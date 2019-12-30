import sys
import time
import argparse
import sqlite3

parser = argparse.ArgumentParser(description='Update the status of the '
        'tasks specified by the supplied TaskIDs (over stdin/pipe).')
parser.add_argument('database')
parser.add_argument('newstatus')
args = parser.parse_args()

with sqlite3.Connection(args.database) as conn:
    c = conn.cursor()
    for line in sys.stdin:
        c.execute('''INSERT INTO history
            VALUES (?, ?, ?)''',
            (int(line), args.newstatus, int(time.time())))
