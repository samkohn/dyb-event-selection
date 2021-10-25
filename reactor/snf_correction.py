"""Compute and upload the spent nuclear fuel (SNF) correction."""
import argparse
import sqlite3

import numpy as np


def compute_correction(energies):
    """Compute the SNF correction for each core.

    Re-implements the procedure from Matt's DybBerkFit
    make_combined_spectra_P17B_unblinded.C file.
    """
    energy_bins = np.linspace(1.75, 4.25, 11)
    bin_width = 0.25
    corrections = np.zeros((820, 6))
    # correction_parameters[core, energy_bin], core = 1 through 6
    correction_parameters = (
        np.array(
            [
                [1.169, 1.317, 1.329, 0.9747, 0.8996, 0.7766, 0.3185,
                    0.2475, 0.03444, 0.0001644, 0.0,
                ],
                [1.186, 1.334, 1.345, 0.9816, 0.9050, 0.7812, 0.3200,
                    0.2487, 0.03460, 0.0001657, 0.0,
                ],
                [1.196, 1.367, 1.412, 1.118, 1.043, 0.9007, 0.3691,
                    0.2870, 0.03996, 0.0002051, 0.0,
                ],
                [1.158, 1.321, 1.362, 1.077, 1.003, 0.8666, 0.3560,
                    0.2769, 0.003857, 0.0002021, 0.0,
                ],
                [1.087, 1.265, 1.342, 1.175, 1.107, 0.9545, 0.3854,
                    0.2998, 0.04178, 0.0002208, 0.0,
                ],
                [0.9552, 1.115, 1.184, 1.009, 0.9525, 0.8192, 0.3258,
                    0.2532, 0.03527, 0.0002063, 0.0,
                ],
            ]
        )
        * 0.01
    )  # Literal values given in percent
    for i, energy in enumerate(energies):
        if energy > energy_bins[-1]:
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
            corrections_up * dE_to_low_edge + corrections_low * dE_from_up_edge
        ) / bin_width
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
            for energy, by_core in zip(energies, correction_factors):
                for i, correction in enumerate(by_core):
                    core = i + 1
                    cursor.execute('''INSERT OR REPLACE INTO snf_correction
                    VALUES (?, ?, ?, ?)''',
                    (energy, core, correction, source))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--update-db')
    args = parser.parse_args()
    main(args.update_db)

