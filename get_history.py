import sys
import argparse
import sqlite3
import datetime

with sqlite3.Connection('example.db') as conn:
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    for line in sys.stdin:
        c.execute('''SELECT * FROM history
            WHERE TaskID = ?
            ORDER BY Timestamp''',
            (int(line.strip()),))
        print('-'*40)
        print('TaskID {}'.format(line.strip()))
        for row in c.fetchall():
            print('{:<24}{:<15}'.format(
                datetime.datetime.fromtimestamp(row['Timestamp']).isoformat(' '),
                row['Status']))
