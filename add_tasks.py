import argparse
import sys
import time
import sqlite3

parser = argparse.ArgumentParser(description='Add tasks from file')
parser.add_argument('database')
parser.add_argument('task_list', default=None, nargs='?')
args = parser.parse_args()

def add_tasks(cursor, tasks):
    cursor.execute('BEGIN EXCLUSIVE TRANSACTION')
    for task in tasks:
        cursor.execute('''INSERT INTO tasks
            (Command) VALUES (?)''',
            (task.strip(),))
        cursor.execute('''INSERT INTO history
            VALUES (last_insert_rowid(),
            "not started", ?)''',
            (int(time.time()),))
    cursor.execute('COMMIT')

with sqlite3.Connection(args.database, isolation_level=None) as conn:
    c = conn.cursor()
    if args.task_list is None:
        add_tasks(c, sys.stdin)
    else:
        with open(args.task_list, 'r') as f:
            add_tasks(c, f)
