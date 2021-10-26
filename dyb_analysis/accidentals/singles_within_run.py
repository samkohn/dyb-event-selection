import argparse

import common

def main(infilename, update_db, run, detector):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    singles = infile.Get('singles')
    interval_ns = 3600e9
    interval_s = interval_ns / 1e9
    n_total = singles.GetEntries()
    # first hour
    singles.GetEntry(0)
    t0 = singles.timestamp
    stop_time = t0 + interval_ns
    n_first_hour = 0
    entry = 0
    while singles.timestamp <= stop_time and entry < n_total:
        n_first_hour += 1
        entry += 1
        singles.GetEntry(entry)
    # last hour
    entry = singles.GetEntries() - 1
    singles.GetEntry(entry)
    t0 = singles.timestamp
    stop_time = t0 - interval_ns
    n_last_hour = 0
    while singles.timestamp >= stop_time and entry >= 0:
        n_last_hour += 1
        entry -= 1
        singles.GetEntry(entry)
    if update_db:
        with common.get_db(update_db) as conn:
            cursor = conn.cursor()
            cursor.execute('''INSERT OR REPLACE INTO singles_within_run
                VALUES (?, ?, ?, ?)''',
                (run, detector, n_first_hour, n_last_hour))
    else:
        print(f'First hour: {n_first_hour} ({n_first_hour/interval_s} Hz)')
        print(f'Last hour: {n_last_hour} ({n_last_hour/interval_s} Hz)')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--update-db')
    parser.add_argument('-r', '--run', type=int)
    parser.add_argument('-d', '--detector', type=int)
    args = parser.parse_args()
    main(args.infile, args.update_db, args.run, args.detector)
