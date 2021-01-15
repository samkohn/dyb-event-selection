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
    parser.add_argument('--sin2', action='store_true')
    parser.add_argument('--dm2', action='store_true')
    parser.add_argument('--filter')
    parser.add_argument('--rate-only', action='store_true')
    parser.add_argument('--avg-near', action='store_true')
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()

    with sqlite3.Connection(args.database) as conn:
        cursor = conn.cursor()
        if args.dm2:
            fit_dm2_query = ', FitDM2ee'
        else:
            fit_dm2_query = ''
        sql_query_str = f'''SELECT TrueSinSqT13, FitSinSqT13, TrueDM2ee {fit_dm2_query}
        FROM fitter_validation_results
        WHERE Category LIKE ?
        AND IsRateOnly = ?
        AND IsAvgNear = ?
        ORDER BY TrueSinSqT13, TrueDm2ee'''
        cursor.execute(sql_query_str, (args.filter, args.rate_only, args.avg_near))
        results = np.array(cursor.fetchall())

    param_pairs = np.unique(results[:, [0, 2]], axis=0)
    fig_index = 1
    figsize=(18, 12)
    if args.sin2:
        fig, axs = plt.subplots(6, 6, figsize=figsize, sharey=True)
        axs_flat = axs.flatten()
        sin2_fig_index = fig_index
        fig_index += 1
    if args.dm2:
        fig_m2, axs_m2 = plt.subplots(6, 6, figsize=figsize, sharey=True)
        axs_m2_flat = axs_m2.flatten()
        m2_fig_index = fig_index
        fig_index += 1
    summary_plot, summary_ax = plt.subplots(figsize=figsize)
    sin2_fits = []
    if args.sin2:
        sin2_errs = []
        sin2_scatters = []
    else:
        sin2_errs = None
        sin2_scatters = None
    dm2_fits = []
    if args.dm2:
        dm2_errs = []
        dm2_scatters = []
    else:
        dm2_errs = None
        dn2_scatters = None
    ax_index = 0
    for sin2, dm2 in param_pairs:
        data = select_config(results, sin2, dm2)
        if args.sin2:
            ax = axs_flat[ax_index]
            values, bins, _ = ax.hist(data[:, 1], range=(0.055, 0.105), bins=50)
            xvals = np.diff(bins)/2 + bins[:-1]  # midpoints
            guess_height = max(values)
            guess_stdev = 0.003
            guess_loc = np.mean(data[:, 1])
            guess_params = (guess_height, guess_loc, guess_stdev)
            fitresults = least_squares(gaus_resid, guess_params, args=(xvals, values))
            fitparams = FitParams(*fitresults.x)
            sin2_fits.append(fitparams.mean)
            sin2_errs.append(fitparams.stdev/np.sqrt(len(data[:, 1])))
            sin2_scatters.append(fitparams.stdev)
            curve_points = np.linspace(bins[0], bins[-1], 100)
            fit_curve = gaus_model(curve_points, fitparams)
            ax.plot(curve_points, fit_curve)
            if sin2 < 0.08:
                horiz_pos = 0.98
                horiz_align = 'right'
            else:
                horiz_pos = 0.02
                horiz_align = 'left'
            text =  rf'''True $\sin^{{2}}2\theta_{{13}}$: {sin2:.5f}
Fit $\mu$: {fitparams.mean:.5f} $\pm$ {fitparams.stdev/np.sqrt(len(data[:, 1])):.5f}
Fit $\sigma$: {fitparams.stdev:.5f}
$\Delta m^{{2}}_{{ee}}$: {dm2:.2e} eV${{}}^2$
Sample $\mu$: {data[:, 1].mean():.5f}
Sample $\sigma$: {data[:, 1].std():.5f}
'''
            if args.verbose:
                print(text)
            ax.text(horiz_pos, 0.98, text,
            fontsize=10,
            transform=ax.transAxes,
            horizontalalignment=horiz_align,
            verticalalignment='top')
        else:
            dm2_fits.append(dm2)

        if args.dm2:
            ax_m2 = axs_m2_flat[ax_index]
            scaled_data = data[:, 3] * 1000  # units of 1e-3 ev^2
            values_m2, bins_m2, _ = ax_m2.hist(scaled_data, range=(2, 3), bins=50)
            xvals_m2 = np.diff(bins_m2)/2 + bins_m2[:-1]  # midpoints
            guess_height_m2 = max(values_m2)
            guess_stdev_m2 = 0.1
            guess_loc_m2 = np.mean(scaled_data)
            guess_params_m2 = (guess_height_m2, guess_loc_m2, guess_stdev_m2)
            fitresults_m2 = least_squares(gaus_resid, guess_params_m2, args=(xvals_m2,
                values_m2))
            fitparams_m2 = FitParams(*fitresults_m2.x)
            dm2_fits.append(fitparams_m2.mean*1e-3)
            dm2_errs.append(fitparams_m2.stdev*1e-3/np.sqrt(len(data[:, 3])))
            dm2_scatters.append(fitparams_m2.stdev*1e-3)
            curve_points_m2 = np.linspace(bins_m2[0], bins_m2[-1], 100)
            fit_curve_m2 = gaus_model(curve_points_m2, fitparams_m2)
            ax_m2.plot(curve_points_m2, fit_curve_m2)
            if dm2 < 2.5e-3:
                horiz_pos = 0.98
                horiz_align = 'right'
            else:
                horiz_pos = 0.02
                horiz_align = 'left'
            ax_m2.text(horiz_pos, 0.98,
rf'''True $\sin^{{2}}2\theta_{{13}}$: {sin2:.4f}
True $\Delta m^{{2}}_{{ee}}$: {dm2:.3e}
$\mu$: {fitparams_m2.mean*1e-3:.4e} $\pm$ '''
rf'''{fitparams_m2.stdev*1e-3/np.sqrt(len(data[:, 3])):.1e} eV${{}}^2$
$\sigma$: {fitparams_m2.stdev*1e-3:.3e} eV${{}}^2$
''',
            fontsize=10,
            transform=ax_m2.transAxes,
            horizontalalignment=horiz_align,
            verticalalignment='top')
            ax_m2.ticklabel_format(axis='x', scilimits=(0, 0), useMathText=True)
        else:
            dm2_fits.append(dm2)

        ax_index += 1
    hack_label_coords = (-0.09, 0)
    if args.sin2:
        axs[-1, 1].set_xlabel(r'Fit $\sin^{2}2\theta_{13}$', fontsize=12)
        plt.figure(sin2_fig_index)
        plt.tight_layout()
        axs[2, 0].set_ylabel('Number of fake experiments', fontsize=12)
        axs[2, 0].yaxis.set_label_coords(*hack_label_coords)
    if args.dm2:
        axs_m2[-1, 1].set_xlabel(r'Fit $\Delta m^{2}_{ee}$ [$10^{-3}$ eV${}^2$]', fontsize=12)
        plt.figure(m2_fig_index)
        plt.tight_layout()
        axs_m2[2, 0].set_ylabel('Number of fake experiments', fontsize=12)
        axs_m2[2, 0].yaxis.set_label_coords(*hack_label_coords)
    # Summary plot
    legend_labels = []
    if args.sin2:
        summary_ax.hist(sin2_fits/param_pairs[:, 0] - 1, histtype='step', linewidth=3)
        legend_labels.append(r'$\sin^{2}2\theta_{13}$')
    if args.dm2:
        summary_ax.hist(dm2_fits/param_pairs[:, 1] - 1, histtype='step', linewidth=3)
        legend_labels.append(r'$\Delta m^{2}_{ee}$')
    summary_ax.legend(legend_labels, fontsize=12)
    summary_ax.set_xlabel('Relative error from mean of 1000 trials', fontsize=12)
    summary_ax.set_ylabel('Number of grid points (parameters)', fontsize=12)
    #summary_ax.plot(param_pairs[:, 1], param_pairs[:, 0], '.')
    #summary_ax.errorbar(dm2_fits, sin2_fits, xerr=dm2_errs, yerr=sin2_errs, fmt='none',
            #capsize=5)
    #plt.subplots_adjust(hspace=0.1, wspace=0.05)
    plt.show()
