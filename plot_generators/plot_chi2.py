"""Plot the chi2 scan, optionally with the 1, 3, 5-sigma regions marked."""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def get_min_max_bounds(delta_chi2s, sin2_values, num_sigmas):
    threshold_delta = num_sigmas**2
    good_enoughs = np.nonzero(delta_chi2s < threshold_delta)[0]
    # Heuristic for getting good approximations to the actual intersection
    min_sin2 = (
        0.5 * (sin2_values[good_enoughs[0]] + sin2_values[good_enoughs[0] - 1])
    )
    max_sin2 = sin2_values[good_enoughs[-1] + 1]
    return min_sin2, max_sin2

def plot_scan(infilename, mark_sigmas=True):
    with np.load(infilename) as infile:
        theta13_values = infile['theta13_values']
        chi2_values = infile['results']
    sin2_values = np.sin(2 * theta13_values)**2
    chi2_min = min(chi2_values)
    delta_chi2s = chi2_values - chi2_min
    fig, ax = plt.subplots()
    ax.plot(sin2_values, chi2_values, 'k')
    ax.grid(False)
    if mark_sigmas:
        for num_sigmas in [1, 3, 5]:
            min_sin2, max_sin2 = get_min_max_bounds(delta_chi2s, sin2_values, num_sigmas)
            threshold_chi2 = chi2_min + num_sigmas**2
            color = '0.7'
            ax.vlines(
                [min_sin2, max_sin2], ymin=0, ymax=threshold_chi2, color=color
            )
            ax.hlines(threshold_chi2, xmin=min_sin2, xmax=max_sin2, color=color)
            ax.annotate(
                rf'{num_sigmas}$\sigma$',
                xy=((min_sin2 + max_sin2)/2, threshold_chi2),
                textcoords='offset points',
                xytext=(0, mpl.rcParams['font.size']/4),
                horizontalalignment='center',
            )
        ax_limit_sigmas = 7
        ax_ylim = ax_limit_sigmas**2 * 1.1
        ax.set_xlim(get_min_max_bounds(delta_chi2s, sin2_values, 7))
        ax.set_ylim([0, ax_ylim])

    ax.set_xlabel(r'$\sin^22\theta_{13}$')
    ax.set_ylabel(r'$\chi^2$')
    return fig, ax

def load_file(name):
    with np.load(name) as infile:
        theta13_values = infile['theta13_values']
        sin2_values = np.sin(2 * theta13_values)**2
        chi2_values = infile['results']
    return (theta13_values, sin2_values, chi2_values)

def plot_many(
    all_filename, stat_filename, other_filenames, other_labels=None, legend=True
):
    theta13_values_all, sin2_values_all, chi2_values_all = load_file(all_filename)
    theta13_values_stat, sin2_values_stat, chi2_values_stat = load_file(stat_filename)
    other_values = {}
    for filename in other_filenames:
        name = os.path.splitext(os.path.basename(filename))[0]
        other_values[name] = load_file(filename)
    fig, ax = plt.subplots()
    stat_handle, = ax.plot(sin2_values_stat, chi2_values_stat, 'k', label='Stat. only')
    if other_labels is None:
        other_labels = other_values.keys()
    other_handles = []
    linestyles = ['--', ':' ,'-.'] * 5
    for label, (_, sin2, chi2), style in zip(
        other_labels,
        other_values.values(),
        linestyles,
    ):
        other_handle, = ax.plot(sin2, chi2, 'k', linestyle=style, linewidth=2, label=label)
        other_handles.append(other_handle)
    all_handle, = ax.plot(sin2_values_all, chi2_values_all, 'k', label='Full uncertainties')
    if legend:
        ax.legend([all_handle, stat_handle] + other_handles)
    ax.grid()
    return fig, ax

def hardcode_plot_sam_results(base_dir):
    all_filename = os.path.join(base_dir, 'grid_results_nominal.npz')
    stat_filename = os.path.join(base_dir, 'grid_results_stat.npz')
    other_filenames = {
        'Reactor': os.path.join(base_dir, 'grid_results_reactor.npz'),
        'Efficiency': os.path.join(base_dir, 'grid_results_efficiency.npz'),
        'Rel. E scale': os.path.join(base_dir, 'grid_results_rel_escale.npz'),
        'Near hall stats': os.path.join(base_dir, 'grid_results_near_stat.npz'),
        'Accidental bkg.': os.path.join(base_dir, 'grid_results_accidental.npz'),
        'Correlated bkg.': os.path.join(base_dir, 'grid_results_corr_bkg.npz'),
        'Oscillation inputs': os.path.join(base_dir, 'grid_results_osc_params.npz'),
    }
    fig, ax = plot_many(
        all_filename,
        stat_filename,
        other_filenames.values(),
        other_filenames.keys(),
        legend=False,
    )
    ax.set_xlim([0.0, 0.135])
    ax.set_ylim([0, 50])
    ax.set_xlabel(r'$\sin^22\theta_{13}$')
    ax.set_ylabel(r'$\chi^2$')
    text_kwargs = dict(
        y=1,
        transform=ax.transAxes,
        horizontalalignment='left',
        verticalalignment='bottom',
        #rotation=15,
        fontsize=12,
    )
    xytext=(0, 10)
    annotate_kwargs = dict(
        xycoords='data',
        textcoords='offset points',
        fontsize=12,
        arrowprops=dict(
            width=1,
            headwidth=5,
            headlength=5,
            shrink=0.05,
            facecolor='black',
        ),
    )
    ax.annotate(
        xy=(0.094, 50),
        text='Stat. only',
        xytext=xytext,
        **annotate_kwargs,
    )
    ax.annotate(
        xy=(0.0273, 50),
        text='Efficiency',
        xytext=xytext,
        **annotate_kwargs,
    )
    ax.annotate(
        xy=(0.1252, 50),
        text='Corr. bkg.',
        xytext=xytext,
        **annotate_kwargs,
    )
    ax.annotate(
        xy=(0.001, 50),
        text='All',
        xytext=(0, 21),
        **annotate_kwargs,
    )
    ax.annotate(
        xy=(0.0036, 50),
        text='Other',
        xytext=xytext,
        **annotate_kwargs,
    )
    return fig, ax


def main(infilenames):
    if len(infilenames) == 1:
        fig, ax = plot_scan(infilenames[0], mark_sigmas=True)
        plt.show()
    else:
        fig, ax = plot_many(infilenames[0], infilenames[1], infilenames[2:])
        plt.show()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('grid_results_file', nargs='+', help='Numpy .npz format')
    args = parser.parse_args()
    main(args.grid_results_file)
