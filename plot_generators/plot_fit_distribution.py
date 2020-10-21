"""Plot a grid of distributions of sinSq2theta13 fitted values."""
import argparse
from dataclasses import dataclass
import sqlite3
import time

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import least_squares


def select_config(arr, true_sin2, true_dm2):
    return arr[(arr[:, 0] == true_sin2) & (arr[:, 2] == true_dm2)]

@dataclass
class FitParams:
    height: float
    mean: float
    stdev: float

gaus = norm.pdf
def gaus_model(x, params):
    return params.height * gaus(x, params.mean, params.stdev)

def gaus_resid(params, xvals, yvals):
    """Expect params = [height, mean, stdev]."""
    params = FitParams(*params)
    model = gaus_model(xvals, params)
    diffs = model - yvals
    errs = np.maximum(np.sqrt(yvals), 1)  # yvals are counts so either 0 or >= 1
    return diffs / errs

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    args = parser.parse_args()

    with sqlite3.Connection(args.database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT TrueSinSqT13, FitSinSqT13, TrueDM2ee, ChiSquare
        FROM fitter_validation_results
        WHERE Category LIKE "%Cori%"
        ORDER BY TrueSinSqT13, TrueDm2ee''')
        results = np.array(cursor.fetchall())

    param_pairs = np.unique(results[:, [0, 2]], axis=0)
    fig, axs = plt.subplots(6, 3, figsize=(12, 14))
    axs_flat = axs.flatten()
    ax_index = 0
    for sin2, dm2 in param_pairs:
        ax = axs_flat[ax_index]
        data = select_config(results, sin2, dm2)
        values, bins, _ = ax.hist(data[:, 1], range=(0.055, 0.105), bins=20)
        xvals = np.diff(bins)/2 + bins[:-1]  # midpoints
        guess_height = max(values)
        guess_stdev = 0.01
        guess_loc = np.mean(data[:, 1])
        guess_params = (guess_height, guess_loc, guess_stdev)
        fitresults = least_squares(gaus_resid, guess_params, args=(xvals, values))
        fitparams = FitParams(*fitresults.x)
        curve_points = np.linspace(bins[0], bins[-1], 100)
        fit_curve = gaus_model(curve_points, fitparams)
        ax.plot(curve_points, fit_curve)
        if sin2 < 0.08:
            horiz_pos = 0.95
            horiz_align = 'right'
        else:
            horiz_pos = 0.05
            horiz_align = 'left'
        ax.text(horiz_pos, 0.95, rf'''Mean: {fitparams.mean:.6f}
True $\sin^{{2}}2\theta_{{13}}$: {sin2:.4f}
True $\Delta m^{{2}}_{{ee}}$: {dm2}
Sigma: {fitparams.stdev:.6f}
''',
        fontsize=10,
        transform=ax.transAxes,
        horizontalalignment=horiz_align,
        verticalalignment='top')
        ax_index += 1
    axs[-1, 1].set_xlabel(r'Fit $\sin^{2}2\theta_{13}$', fontsize=12)
    axs[2, 0].set_ylabel('Number of fake experiments', fontsize=12)
    plt.tight_layout()
    plt.show()
