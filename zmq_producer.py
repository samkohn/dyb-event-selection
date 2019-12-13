import zmq
import time
import random
import sqlite3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('timelimit_min', type=float)
args = parser.parse_args()
timelimit_min = args.timelimit_min
timelimit_sec = 60 * args.timelimit_min
start_time = time.time()
finish_time = start_time + timelimit_sec

context = zmq.Context()
producer = context.socket(zmq.REP)
producer.bind('tcp://*:52837')
print('bound')
conn = sqlite3.connect('example.db', isolation_level=None)
conn.row_factory = sqlite3.Row
c = conn.cursor()

def current_time():
    return int(time.time())

def next_task():
    try:
        c.execute('BEGIN EXCLUSIVE TRANSACTION')
        c.execute('''SELECT TaskID, Command, h.Status
            FROM tasks JOIN history AS h USING (TaskID)
            WHERE Timestamp = (SELECT MAX(Timestamp)
                FROM history as h2
                WHERE h2.TaskID = h.TaskID)
                AND h.Status = "not started"
            ORDER BY TaskID ASC
            LIMIT 1''')
        row = c.fetchone()
        if row is None:
            c.execute('COMMIT')
            conn.commit()
            return None
        else:
            c.execute('''INSERT INTO history
                VALUES (?, "started", ?)''',
                (row['TaskID'], current_time()))
            c.execute('COMMIT')
            conn.commit()
            return row
    except sqlite3.Error:
        conn.commit()
        return None

def finish_task(task_id):
    c.execute('''INSERT INTO history
        VALUES (?, "finished", ?)''',
        (task_id, current_time()))
    conn.commit()

def error_task(task_id):
    c.execute('''INSERT INTO history
        VALUES (?, "error", ?)''',
        (task_id, current_time()))
    conn.commit()

task_id = 0
while task_id is not None:
    request = producer.recv()
    worker_done = request[-4:] == b'DONE'
    if worker_done:
        finished_id = int(request[:-4])
    else:
        finished_id = int(request)
    print('received request')
    if finished_id < 0:
        error_task(-1 * finished_id)
    elif finished_id > 0:
        finish_task(finished_id)
    else:
        pass
    if worker_done or time.time() > finish_time:
        producer.send_multipart([b'0', b'DONE'])
    else:
        task = next_task()
        if task is None:
            producer.send_multipart([b'0', b'DONE'])
        else:
            producer.send_multipart([b'%d' % task['TaskID'],
                task['Command'].encode()])
        print('sent job command')
