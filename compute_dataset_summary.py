"""Compute total dataset muon and multiplicity efficiencies and livetime."""
import numpy as np

import common

def daq_livetime_s(database):
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
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            '''
        )
        livetimes_s = np.array(cursor.fetchall()).reshape(-1)
    return livetimes_s

def daq_livetime_days(database):
    """Return an array of DAQ livetimes in days from EH1-AD1 to EH3-AD4."""
    livetime_s = daq_livetime_s(database)
    return livetime_s/60/60/24

def muon_efficiency(database):
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
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            '''
        )
        muon_effs = np.array(cursor.fetchall()).reshape(-1)
    return muon_effs

def multiplicity_efficiency(database):
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
                DetNo
            )
            GROUP BY
                Hall,
                DetNo
            ORDER BY
                Hall,
                DetNo
            '''
        )
        mult_effs = np.array(cursor.fetchall()).reshape(-1)
    return mult_effs

