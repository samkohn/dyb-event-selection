"""Dump the AmC spectrum to the database."""
import argparse

import numpy as np

import common

def main(infilename, hist_name, database, label, binning_id, binning_db_path):
    import ROOT
    # fetch binning
    with common.get_db(binning_db_path) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                BinEdgeEnergy_keV
            FROM
                reco_binnings
            WHERE
                Id = ?
            ORDER BY BinEdgeIndex
            ''',
            (binning_id,)
        )
        # Must reshape to get 1D array since sqlite3 returns 2D
        bin_edges = np.array(cursor.fetchall(), dtype=float).reshape(-1)/1000
    infile = ROOT.TFile(infilename, 'READ')
    bg_spec = infile.Get(hist_name)
    # Annoyingly, the AmC spectrum is stored as a TF1, so I have to bin it manually
    values = np.zeros((len(bin_edges) - 1,))
    for bin_index, (low_edge, up_edge) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
        values[bin_index] = bg_spec.Integral(low_edge, up_edge)
    total_counts = sum(values)
    values /= total_counts  # normalize
    rows = []
    for bin_index, value in enumerate(values):
        rows.append(
            (label, binning_id, bin_index, value)
        )
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.executemany('''
            INSERT OR REPLACE INTO
                amc_spectrum
            VALUES
                (?, ?, ?, ?)
            ''',
            rows
        )
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('database')
    parser.add_argument('--label')
    parser.add_argument('--binning-id', type=int)
    parser.add_argument('--binning-db-path')
    parser.add_argument('--hist-name')
    args = parser.parse_args()
    main(
        args.infile,
        args.hist_name,
        args.database,
        args.label,
        args.binning_id,
        args.binning_db_path,
    )
