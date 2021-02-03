"""Compute the statistical uncertainty for accidentals."""
import argparse
from collections import defaultdict
import json

import numpy as np

import common
import compute_dataset_summary as summary

ERROR_COEFFICIENT = 1.1

def baserate_uncertainty(database, label, general_label):
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                Hall,
                DetNo,
                RunNo,
                acc.Label,
                singles.Rate_Hz,
                singles.IsolatedEventCount,
                singles.IsolatedEventRate_Hz,
                acc.BaseRate_Hz
            FROM
                singles_rates AS singles
            INNER JOIN
                accidental_subtraction AS acc
            USING (
                RunNo,
                DetNo
            )
            NATURAL JOIN
                runs
            WHERE
                acc.Label = ?
                AND singles.Label = ?
            ''',
            (label, general_label)
        )
        inputs = np.array(cursor.fetchall())
    # Annoying casts because of numpy's homogeneous array types
    halls = inputs[:, 0].astype(int)
    dets = inputs[:, 1].astype(int)
    runs = inputs[:, 2].astype(int)
    labels = inputs[:, 3].astype(str)
    singles_rate = inputs[:, 4].astype(float)
    n_1fold = inputs[:, 5].astype(float)
    r_1fold = inputs[:, 6].astype(float)
    r_ss = inputs[:, 7].astype(float)
    base_rate_error = (
        r_ss * (ERROR_COEFFICIENT * r_1fold / singles_rate + 1) / np.sqrt(n_1fold)
    )
    result_rows = np.vstack((base_rate_error, runs, dets, labels)).T
    print(result_rows.shape)
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.executemany('''
            UPDATE
                accidental_subtraction
            SET
                BaseRate_Hz_error = ?
            WHERE
                RunNo = ?
                AND DetNo = ?
                AND Label = ?
            ''',
            result_rows
        )

def count_errors(database, label, general_label):
    """Return the statistical error on the number of accidentals in each AD.

    Return value is an array sorted from EH1-AD1 to EH3-AD4.
    """
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT  -- easier to compute variance in SQL then take sqrt later
                SUM(  -- add in quadrature
                    BaseRate_Hz * Total_Acc_Eff * Livetime_ns/1e9 *  -- num acc squared
                    BaseRate_Hz * Total_Acc_Eff * Livetime_ns/1e9 *
                    (   -- Add relative errors of R_ss and eps_total in quadrature
                        -- (R_ss = BaseRate_Hz, eps_total = Total_Acc_Eff)
                        BaseRate_Hz_error * BaseRate_Hz_error / (
                            BaseRate_Hz * BaseRate_Hz
                        )
                        + Total_Acc_Eff_err * Total_Acc_Eff_err / (
                            Total_Acc_Eff * Total_Acc_Eff
                        )
                    )
                )
            FROM
                accidental_subtraction AS acc
            INNER JOIN
                muon_rates as mu
            USING (
                RunNo,
                DetNo
            )
            NATURAL JOIN
                runs
            WHERE
                acc.Label = ?
                AND mu.Label = ?
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            ''',
            (label, general_label)
        )
        result = np.sqrt(np.array(cursor.fetchall()).reshape(-1))
    return result

def counts(database, label, general_label):
    """Return an array of predicted accidentals counts."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                SUM(acc.BaseRate_Hz * Total_Acc_Eff * Livetime_ns/1e9)
            FROM
                accidental_subtraction AS acc
            NATURAL JOIN
                runs
            INNER JOIN
                muon_rates as mu
            USING (
                RunNo,
                DetNo
            )
            WHERE
                acc.Label = ?
                AND mu.Label = ?
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            ''',
            (label, general_label)
        )
        counts = np.array(cursor.fetchall()).reshape(-1)
    return counts

def rates(database, label, general_label):
    """Return a 2D array of accidentals rates, corrected for eps_mu and eps_m.

    The resulting rates should be multiplied by DAQ_livetime * eps_mu * eps_m
    to get back the predicted counts.
    """
    acc_counts = counts(database, label, general_label)
    daq_livetime = summary.daq_livetime_days(database, general_label)
    muon_eff = summary.muon_efficiency(database, general_label)
    multiplicity_eff = summary.multiplicity_efficiency(database, general_label)

    corrected_counts = acc_counts / muon_eff / multiplicity_eff
    corrected_rates = corrected_counts / daq_livetime

    return corrected_rates

def rate_errors(database, label, general_label):
    """Return a 2D array of accidentals rate errors, corrected for eps_mu and eps_m.

    The resulting rate errors should be multiplied by DAQ_livetime * eps_mu * eps_m
    to get back the predicted count errors.
    """
    errors = count_errors(database, label, general_label)
    daq_livetime = summary.daq_livetime_days(database, general_label)
    muon_eff = summary.muon_efficiency(database, general_label)
    multiplicity_eff = summary.multiplicity_efficiency(database, general_label)

    corrected_errors = errors / muon_eff / multiplicity_eff
    corrected_rate_errors = corrected_errors / daq_livetime

    return corrected_rate_errors

def dump_accidentals_json(database, label, general_label, filename):
    acc_rates = rates(database, label, general_label)
    acc_rate_errors = rate_errors(database, label, general_label)
    print(acc_rates)
    print(acc_rate_errors)
    acc_dict = defaultdict(dict)
    for (hall, det), rate, error in zip(common.all_ads, acc_rates, acc_rate_errors):
        acc_dict[f'EH{hall}-AD{det}']['rate'] = rate
        acc_dict[f'EH{hall}-AD{det}']['error'] = error
    with open(filename, 'w') as fout:
        json.dump({'accidental': acc_dict}, fout, indent=2)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('--base-rate', action='store_true')
    parser.add_argument('--update-db')
    args = parser.parse_args()
    if args.base_rate:
        baserate_uncertainty(args.database)
