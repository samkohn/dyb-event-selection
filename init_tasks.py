import argparse
import sqlite3

def main(db_location):
    with sqlite3.Connection(db_location) as conn:
        c = conn.cursor()
        c.execute('''
            CREATE TABLE tasks (
                TaskID INTEGER PRIMARY KEY ASC AUTOINCREMENT,
                Command VARCHAR(1023)
            )
        ''')
        c.execute('''
            CREATE TABLE history (
                TaskID INTEGER,
                Status VARCHAR(15),
                Timestamp INTEGER
            )
        ''')
    print('Created tables tasks, history at location {}'.format(db_location))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('db_location')
    args = parser.parse_args()
    main(args.db_location)
