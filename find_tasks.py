from __future__ import print_function
import argparse
import sqlite3

parser = argparse.ArgumentParser(description='Get the TaskID and/or Status '
        'for the requested tasks. If both --taskid and --status are given, '
        'each line of output will contain both the TaskID and Status for '
        'a particular task.')
parser.add_argument('database')
group = parser.add_mutually_exclusive_group()
group.add_argument('--search',
        help='SQLite match string (WHERE Command LIKE <search>)')
group.add_argument('--runfile', type=int, nargs=2,
        help='Run and file numbers for a particular command')
group.add_argument('--run', type=int, default=0,
        help='All the commands with the given run number')
group.add_argument('--all', action='store_true',
        help='All tasks')
parser.add_argument('--taskid', action='store_true',
        help='Print the TaskID for matching tasks, one per line')
parser.add_argument('--status', action='store_true',
        help='Print the Status for matching tasks, one per line')
parser.add_argument('--time', action='store_true',
        help='Print the last updated Timestamp for matching tasks, one per line')
args = parser.parse_args()
if args.search is not None:
    search = args.search
elif args.run != 0:
    search = '%{:07}%'.format(args.run)
elif args.all:
    search = '%'
else:
    search = '%{:07}%{:04}%'.format(*args.runfile)

toprint = []
if args.taskid:
    toprint.append(0)
if args.status:
    toprint.append(1)
if args.time:
    toprint.append(2)

with sqlite3.Connection(args.database) as conn:
    c = conn.cursor()
    if args.taskid and not args.status and not args.time:
        c.execute('SELECT TaskID FROM tasks WHERE Command LIKE ?', (search,))
        for row in c.fetchall():
                print(row[0])
    else:
        c.execute('''SELECT TaskID, Status, Timestamp
            FROM (
                SELECT *,
                ROW_NUMBER() OVER (
                    PARTITION BY TaskID
                    ORDER BY Timestamp DESC) rn
                FROM tasks JOIN history USING (TaskID)
            )
            WHERE rn = 1
            AND Command LIKE ?
            ORDER BY TaskID''',
            (search,))
        for row in c.fetchall():
            print(*[row[index] for index in toprint])
