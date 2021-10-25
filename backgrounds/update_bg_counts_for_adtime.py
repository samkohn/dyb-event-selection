"""Use the DT cut efficiency to update the background counts for correlated bgs."""

import argparse

import numpy as np

import common


def main(
    database,
    bg_database,
    adsimple_label,
    adtime_label,
    adsimple_bg_source,
    new_bg_source,
):
    # Get the DT efficiencies
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                Efficiency
            FROM
                distance_time_cut_efficiency
            WHERE
                Source = ?
            ORDER BY
                Hall,
                DetNo
            ''',
            (adsimple_label,)
        )
        adsimple_DT_effs = np.array(cursor.fetchall()).reshape(-1)
        cursor.execute('''
            SELECT
                Efficiency
            FROM
                distance_time_cut_efficiency
            WHERE
                Source = ?
            ORDER BY
                Hall,
                DetNo
            ''',
            (adtime_label,)
        )
        adtime_DT_effs = np.array(cursor.fetchall()).reshape(-1)
    correction_factors = adtime_DT_effs / adsimple_DT_effs
    with common.get_db(bg_database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                Count, Error
            FROM
                bg_counts
            WHERE
                Label = ?
                AND BgName <> "accidental"
            ORDER BY
                BgName,
                Hall,
                DetNo
            ''',
            (adsimple_bg_source,)
        )
        adsimple_bg_counts_errors = np.array(cursor.fetchall())
        cursor.execute('''
            SELECT
                BgName
            FROM
                bg_counts
            WHERE
                Label = ?
                AND BgName <> "accidental"
            ORDER BY
                BgName,
                Hall,
                DetNo
            ''',
            (adsimple_bg_source,)
        )
        # Don't really need np.array here but it's convinient for
        # flattening the list
        adsimple_bg_names = np.array(cursor.fetchall()).reshape(-1)
        num_bgs = len(np.unique(adsimple_bg_names))
        # tile by (num_bgs, 2) so that each type of bg and each error gets corrected
        tiled_corrections = np.tile(correction_factors.reshape((8, 1)), (num_bgs, 2))
        corrected_counts = adsimple_bg_counts_errors * tiled_corrections
        rows = []
        for name, (corr_count, corr_err), (hall, det) in zip(
            adsimple_bg_names, corrected_counts, common.all_ads * num_bgs
        ):
            rows.append((new_bg_source, hall, det, name, corr_count, corr_err))
        print(rows)
        cursor.executemany('''
            INSERT OR REPLACE INTO
                bg_counts
            VALUES
                (?, ?, ?, ?, ?, ?)
            ''',
            rows
        )
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('bg_database')
    parser.add_argument('--adsimple-label')
    parser.add_argument('--adtime-label')
    parser.add_argument('--adsimple-bg-source')
    parser.add_argument('--new-bg-source')
    args = parser.parse_args()
    main(
        args.database,
        args.bg_database,
        args.adsimple_label,
        args.adtime_label,
        args.adsimple_bg_source,
        args.new_bg_source,
    )

