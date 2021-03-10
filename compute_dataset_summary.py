"""Compute total dataset muon and multiplicity efficiencies and livetime."""
import numpy as np

import common

def daq_livetime_s(database, label):
    """Return an array of DAQ livetimes ordered from EH1-AD1 to EH3-AD4."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        # Total DAQ Livetime
        cursor.execute('''
            SELECT
                SUM(Livetime_ns/Efficiency/1e9)
            FROM
                muon_rates
            NATURAL JOIN
                runs
            WHERE
                Label = ?
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            ''',
            (label,)
        )
        livetimes_s = np.array(cursor.fetchall()).reshape(-1)
    return livetimes_s

def daq_livetime_days(database, label):
    """Return an array of DAQ livetimes in days from EH1-AD1 to EH3-AD4."""
    livetime_s = daq_livetime_s(database, label)
    return livetime_s/60/60/24

def unvetoed_livetime_s(database, label):
    """Return an array of unvetoed livetimes ordered from EH1-AD1 to EH3-AD4."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        # Total DAQ Livetime
        cursor.execute('''
            SELECT
                SUM(Livetime_ns/1e9)
            FROM
                muon_rates
            NATURAL JOIN
                runs
            WHERE
                Label = ?
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            ''',
            (label,)
        )
        livetimes_s = np.array(cursor.fetchall()).reshape(-1)
    return livetimes_s

def muon_efficiency(database, label):
    """Return an array of muon efficiencies from EH1-AD1 to EH3-AD4."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                SUM(Livetime_ns)/SUM(Livetime_ns/Efficiency)
            FROM
                muon_rates
            NATURAL JOIN
                runs
            WHERE
                Label = ?
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            ''',
            (label,)
        )
        muon_effs = np.array(cursor.fetchall()).reshape(-1)
    return muon_effs

def muon_total_counts(database, label):
    """Return an array of muon counts from EH1-AD1 to EH3-AD4."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                SUM(`Count`)
            FROM
                muon_rates
            NATURAL JOIN
                runs
            WHERE
                Label = ?
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            ''',
            (label,)
        )
        muon_counts = np.array(cursor.fetchall()).reshape(-1)
    return muon_counts

def muon_rate_Hz(database, label):
    """Return an array of muon rates from EH1-AD1 to EH3-AD4."""
    unvetoed_livetimes = unvetoed_livetime_s(database, label)
    counts = muon_total_counts(database, label)
    return counts / unvetoed_livetimes

def multiplicity_efficiency(database, label):
    """Return an array of multiplicity efficiencies from EH1-AD1 to EH3-AD4."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                SUM(MultiplicityVetoEfficiency * Livetime_ns/Efficiency)/
                SUM(Livetime_ns/Efficiency)
            FROM
                singles_rates
            NATURAL JOIN
                runs
            INNER JOIN
                muon_rates
            USING (
                RunNo,
                DetNo,
                Label
            )
            WHERE
                Label = ?
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            ''',
            (label,)
        )
        mult_effs = np.array(cursor.fetchall()).reshape(-1)
    return mult_effs

def singles_rate_Hz(database, label):
    """Return an array of singles rates from EH1-AD1 to EH3-AD4."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                SUM(singles.Rate_Hz * Livetime_ns/Efficiency)/
                SUM(Livetime_ns/Efficiency)
            FROM
                singles_rates AS singles
            NATURAL JOIN
                runs
            INNER JOIN
                muon_rates
            USING (
                RunNo,
                DetNo,
                Label
            )
            WHERE
                Label = ?
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            ''',
            (label,)
        )
        singles_rates = np.array(cursor.fetchall()).reshape(-1)
    return singles_rates

def coincidences_counts(database, label):
    """Return an array of coincidence counts from EH1-AD1 to EH3-AD4."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                SUM(NumCoincidences)
            FROM
                num_coincidences_by_run
            NATURAL JOIN
                runs
            WHERE
                Label = ?
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            ''',
            (label,)
        )
        counts = np.array(cursor.fetchall()).reshape(-1)
    return counts

def coincidences_rates(database, label, general_label):
    """Return an array of coincidence rates (per day) from EH1-AD1 to EH3-AD4.

    Rates are corrected for muon and multiplicity efficiency.
    """
    counts = coincidences_counts(database, label)
    daq_livetimes = daq_livetime_days(database, general_label)
    mult_effs = multiplicity_efficiency(database, general_label)
    muon_effs = muon_efficiency(database, general_label)
    rates = counts / mult_effs / muon_effs / daq_livetimes
    return rates

