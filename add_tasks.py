import argparse
import time
import sqlite3

parser = argparse.ArgumentParser(description='Add tasks from file')
parser.add_argument('database')
parser.add_argument('task_list')
args = parser.parse_args()

with open(args.task_list, 'r') as f:
    with sqlite3.Connection(args.database, isolation_level=None) as conn:
        c = conn.cursor()
        c.execute('BEGIN EXCLUSIVE TRANSACTION')
        for row in f:
            c.execute('''INSERT INTO tasks
                (Command) VALUES (?)''',
                (row.strip(),))
            c.execute('''INSERT INTO history
                VALUES (last_insert_rowid(),
                "not started", ?)''',
                (int(time.time()),))
        c.execute('COMMIT')
