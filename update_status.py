import sys
import time
import argparse
import sqlite3

parser = argparse.ArgumentParser()
parser.add_argument('newstatus')
args = parser.parse_args()

with sqlite3.Connection('example.db') as conn:
    c = conn.cursor()
    for line in sys.stdin:
        c.execute('''INSERT INTO history
            VALUES (?, ?, ?)''',
            (int(line), args.newstatus, int(time.time())))
