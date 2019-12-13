from __future__ import print_function
import argparse
import sqlite3

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument('--search')
group.add_argument('--runfile', type=int, nargs=2)
group.add_argument('--run', type=int, default=0)
group.add_argument('--all', action='store_true')
parser.add_argument('--taskid', action='store_true')
parser.add_argument('--status', action='store_true')
args = parser.parse_args()
if args.search is not None:
    search = args.search
elif args.run != 0:
    search = '%{:07}%'.format(args.run)
elif args.all:
    search = '%'
else:
    search = '%{:07}%{:04}%'.format(*args.runfile)

with sqlite3.Connection('example.db') as conn:
    c = conn.cursor()
    if args.taskid and not args.status:
        c.execute('SELECT TaskID FROM tasks WHERE Command LIKE ?', (search,))
        for row in c.fetchall():
                print(row[0])
    else:
        c.execute('''SELECT TaskID, Status
            FROM history AS h JOIN tasks USING (TaskID)
            WHERE Timestamp = (SELECT MAX(Timestamp)
                FROM history as h2
                WHERE h.TaskID = h2.TaskID)
            AND Command LIKE ?
            ORDER BY h.rowid DESC
            LIMIT 1''',
            (search,))
        if args.taskid and args.status:
            for row in c.fetchall():
                print(*row)
        elif args.status:
            for row in c.fetchall():
                print(row[1])
