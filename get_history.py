import sys
import argparse
import sqlite3
import datetime

parser = argparse.ArgumentParser(description='Get history for given TaskIDs. '
        'TaskIDs can be provided as command-line arguments or piped in '
        'as output from another command, one TaskID per line.')
parser.add_argument('database')
parser.add_argument('tasks', nargs='*', metavar='task')
args = parser.parse_args()

with sqlite3.Connection(args.database) as conn:
    conn.row_factory = sqlite3.Row
    c = conn.cursor()
    if len(args.tasks) == 0:
        task_ids = sys.stdin
    else:
        task_ids = args.tasks
    for task_id_str in task_ids:
        c.execute('''SELECT * FROM history
            WHERE TaskID = ?
            ORDER BY Timestamp''',
            (int(task_id_str.strip()),))
        print('-'*40)
        print('TaskID {}'.format(task_id_str.strip()))
        for row in c.fetchall():
            print('{:<24}{:<15}'.format(
                datetime.datetime.fromtimestamp(row['Timestamp']).isoformat(' '),
                row['Status']))
