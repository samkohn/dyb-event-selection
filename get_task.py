import argparse
import sqlite3

parser = argparse.ArgumentParser(description='Get task command and info')
parser.add_argument('database')
parser.add_argument('tasks', nargs='*', metavar='task')
args = parser.parse_args()

with sqlite3.Connection(args.database) as conn:
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    if len(args.tasks) == 0:
        task_ids = sys.stdin
    else:
        task_ids = args.tasks
    for task_id_str in task_ids:
        cursor.execute('''SELECT * FROM tasks
            WHERE TaskID = ?''',
            (int(task_id_str.strip()),))
        print('-'*40)
        print('TaskID {}'.format(task_id_str.strip()))
        row = cursor.fetchone()
        print(row['Command'])
