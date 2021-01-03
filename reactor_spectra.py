"""Read and store the reactor spectra."""
import argparse

import common

def main(infilename, database, source, data_period):
    with open(infilename, 'r') as f:
        with common.get_db(database) as conn:
            cursor = conn.cursor()
            for line in f:
                items = [float(x) for x in line.split()]
                energy = items[0]
                for core_index, value in enumerate(items[1:]):
                    cursor.execute('''
                        INSERT OR REPLACE INTO reactor_emitted_spectrum
                        VALUES (?, ?, ?, ?, ?)''',
                        (energy, core_index + 1, value, data_period, source)
                    )




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('database')
    parser.add_argument('source')
    parser.add_argument('data_period', choices=('6ad', '8ad', '7ad'))
    args = parser.parse_args()
    print(args)
    main(args.infile, args.database, args.source, args.data_period)
