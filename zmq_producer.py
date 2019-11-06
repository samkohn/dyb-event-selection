import zmq
import random
import sqlite3

context = zmq.Context()
producer = context.socket(zmq.REP)
producer.bind('tcp://*:52837')
print('bound')
conn = sqlite3.connect('example.db')
conn.row_factory = sqlite3.Row
c = conn.cursor()

def next_task():
    lock_id = random.randint(1, 2**32)
    c.execute('''UPDATE tasks
    SET Status = "locked", Lock = ?
    WHERE TaskID = (
    SELECT TaskID FROM tasks
    WHERE Status = "not started"
    ORDER BY TaskID ASC
    LIMIT 1)
    ''', (lock_id,))
    conn.commit()
    if c.rowcount == 1:
        c.execute('''SELECT *
        FROM tasks
        WHERE Lock = ?
        ''', (lock_id,))
        task = c.fetchone()
        c.execute('''UPDATE tasks
        SET Status = "started", Lock = NULL
        WHERE Lock = ?
        ''', (lock_id,))
        conn.commit()
        print('Fetched task')
        print(task)
        return task
    elif c.rowcount == 0:
        return None
    else:
        raise RuntimeError('Lock mechanism failed')

def finish_task(task_id):
    c.execute('''UPDATE tasks
    SET Status = "finished"
    WHERE TaskID = ?
    ''', (task_id,))
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
    if finished_id != 0:
        finish_task(finished_id)
    if worker_done:
        producer.send_multipart([b'0', b'DONE'])
    else:
        task = next_task()
        if task is None:
            producer.send_multipart([b'0', b'DONE'])
        else:
            producer.send_multipart([b'%d' % task['TaskID'],
                task['Command'].encode()])
        print('sent job command')
