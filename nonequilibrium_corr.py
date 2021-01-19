"""Compute and upload the reactor nonequilibrium correction."""
import argparse
import sqlite3

import numpy as np

def compute_correction(energies):
    """Compute the nonequilibrium correction.

    Re-implements the procedure from Matt's DybBerkFit
    make_combined_spectra_P17B_unblinded.C file.
    """
    energy_bins = np.array([2, 2.5, 3, 3.5, 4])
    bin_width = 0.5
    corrections = np.zeros((820, 4))
    # Each row is an isotope: U235, U238, Pu239, Pu241
    # Each column is an energy bin corresponding to the above bins
    correction_parameters = np.array([
        [5.7, 4.4, 1.5, 0.7, 0.1],
        [0, 0, 0, 0, 0],
        [2.1, 1.7, 0.5, 0, 0],
        [1.9, 1.5, 0.5, 0, 0],
    ]) * 0.01  # Literal values given in percent
    for i, energy in enumerate(energies):
        if energy > 4:
            corrections[i, :] = 0
            continue
        # Need the >= to properly account for the upper bin edge.
        # All other bin edges are properly handled anyways.
        bin_low_index = max(0, np.argmax(energy_bins >= energy) - 1)
        bin_high_index = bin_low_index + 1
        dE_to_low_edge = energy - energy_bins[bin_low_index]
        dE_from_up_edge = energy_bins[bin_high_index] - energy
        corrections_low = correction_parameters[:, bin_low_index]
        corrections_up = correction_parameters[:, bin_high_index]
        correction = (
                (corrections_up * dE_to_low_edge + corrections_low * dE_from_up_edge)
                / bin_width
        )
        corrections[i] = correction

    return 1 + corrections

def main(update_db):
    energies = np.linspace(1.8, 9.99, 820)
    correction_factors = compute_correction(energies)
    if update_db is None:
        print(correction_factors[:10])
    else:
        source = "DybBerkFit make_combined_spectra_P17B_unblinded.C"
        with sqlite3.Connection(update_db) as conn:
            cursor = conn.cursor()
            for energy, (u235, u238, pu239, pu241) in zip(energies, correction_factors):
                cursor.execute('''INSERT OR REPLACE INTO nonequilibrium
                VALUES (?, ?, ?, ?, ?, ?)''',
                (energy, u235, u238, pu239, pu241, source))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--update-db')
    args = parser.parse_args()
    main(args.update_db)

