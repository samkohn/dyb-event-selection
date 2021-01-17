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


def extrapolate_down_one_bin(database, source, data_period):
    """Linear extrapolation based on 0th and 1st bin values."""
    with common.get_db(database) as conn:
        conn.row_factory = common.sqlite3.Row
        cursor = conn.cursor()
        for core in range(1, 7):
            cursor.execute('''
                SELECT
                    Energy, NuPerMeVPerSec
                FROM
                    reactor_emitted_spectrum
                WHERE
                    Source = ?
                    AND DataPeriod = ?
                    AND Core = ?
                ORDER BY Energy
                LIMIT 2
                ''',
                (source, data_period, core)
            )
            bin_0, bin_1 = cursor.fetchall()
            dE = bin_1['Energy'] - bin_0['Energy']
            dN = bin_1['NuPerMeVPerSec'] - bin_0['NuPerMeVPerSec']
            new_E = bin_0['Energy'] - dE
            new_N = bin_0['NuPerMeVPerSec'] - dN
            cursor.execute('''
                INSERT OR REPLACE INTO
                    reactor_emitted_spectrum
                VALUES
                    (?, ?, ?, ?, ?)
                ''',
                (new_E, core, new_N, data_period, source)
            )
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('database')
    parser.add_argument('source')
    parser.add_argument('data_period', choices=('6ad', '8ad', '7ad'))
    args = parser.parse_args()
    print(args)
    main(args.infile, args.database, args.source, args.data_period)
