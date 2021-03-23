"""Plot the fit results including prompt spectrum, L/E, backgrounds, etc."""

import argparse
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import common
import fit
import prediction

def _average_near_hall_predictions(predicted):
    # Average the near-hall predictions
    predicted_avg = {halldet: 0 for halldet in common.far_ads}
    denominator = {halldet: 0 for halldet in common.far_ads}
    for (far_halldet, near_halldet), n_predicted in predicted.items():
        predicted_avg[far_halldet] += n_predicted
        denominator[far_halldet] += 1
    for far_halldet, n_predicted in predicted_avg.items():
        predicted_avg[far_halldet] = n_predicted / denominator[far_halldet]
    return predicted_avg

def _plot_line_hist(ax, bin_edges, vals, endpoints=True, **kwargs):
    if endpoints:
        vals = np.concatenate(([0], vals, [0]))
        bin_edges = np.concatenate(([bin_edges[0]], bin_edges))
    else:
        vals = np.concatenate((vals, [vals[-1]]))
    return ax.plot(bin_edges, vals, drawstyle='steps-post', **kwargs)

def _plot_point_hist(ax, bin_edges, vals, **kwargs):
    """Plot a "step" histogram with no vertical lines connecting the bins.

    Plotting algorithm based on https://stackoverflow.com/a/44961513/1144588.
    A NaN value is inserted between each value.
    Since NaNs are not plotted, the line jumps to the next "real" value.
    """
    bin_centers = bin_edges[:-1] + 0.5 * np.diff(bin_edges)
    X = np.c_[bin_edges[:-1], bin_centers, bin_edges[1:], bin_edges[1:]].flatten()
    Y = np.c_[vals, vals, vals, np.zeros_like(bin_edges[:-1])*np.nan].flatten()
    if 'yerr' in kwargs:
        plain_err = kwargs.pop('yerr')
        err = np.c_[plain_err, plain_err, plain_err, np.zeros_like(bin_edges[:-1])].flatten()
        return ax.errorbar(X, Y, yerr=err, errorevery=(1, 4), **kwargs)
    else:
        return ax.plot(X, Y, **kwargs)

def _spectrum_ratio_plot(
    bin_edges,
    line_hist_dicts,
    errorbar_dicts,
    ratio_ylim,
    line_hist_labels,
    errorbar_labels,
):
    fig, axs_deep = plt.subplots(
        4, 2, gridspec_kw={'height_ratios': [3, 1, 3, 1], 'hspace': 0},
        sharex=True
    )
    axs_flat = axs_deep.T.flatten()
    bin_centers = bin_edges[:-1] + 0.5 * np.diff(bin_edges)
    for i, halldet in enumerate(common.far_ads):
        name = f'EH{halldet[0]}-AD{halldet[1]}'
        ax_index_main = 2 * i
        ax_index_ratio = ax_index_main + 1
        #ax.plot(reco_bin_centers, far_no_osc[halldet], '.')
        ax = axs_flat[ax_index_main]
        for line_hist_dict, label in zip(line_hist_dicts, line_hist_labels):
            _plot_line_hist(ax, bin_edges, line_hist_dict[halldet], label=label)
        for (data, errs), label in zip(errorbar_dicts, errorbar_labels):
            obs = data[halldet]
            err = errs[halldet]
            _plot_point_hist(ax, bin_edges, obs, yerr=err, elinewidth=2, label=label)
            #ax.errorbar(
                #bin_centers, obs, yerr=err, fmt='_', elinewidth=2, markersize=5,
                #label=label
            #)
        ax.set_ylim([0, ax.get_ylim()[1]])
        ax.grid()
        ax.tick_params(axis='both', which='major', labelsize=14)
        #ax.text(0.82, 0.95, name, fontsize=12, transform=ax.transAxes,
            #verticalalignment='top')
        ax.legend(
            fontsize=12,
            loc='upper right',
            title=name,
            title_fontsize=14
        )
        ax = axs_flat[ax_index_ratio]
        denominator = line_hist_dicts[0]
        for line_hist_dict in line_hist_dicts:
            _plot_line_hist(
                ax, bin_edges, line_hist_dict[halldet]/denominator[halldet],
                endpoints=False
            )
        for data, errs in errorbar_dicts:
            obs = data[halldet]/denominator[halldet]
            err = errs[halldet]/denominator[halldet]
            _plot_point_hist(ax, bin_edges, obs, yerr=err, elinewidth=2, label=label)
            #ax.errorbar(bin_centers, obs, yerr=err, fmt='_', elinewidth=2, markersize=5)
        ax.set_ylim(ratio_ylim)
        ax.grid()
        ax.set_xlabel('Prompt energy [MeV]', fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=14)
    return fig

def plot_prompt_spectrum(constants, fit_params):
    far_best_fits_ad_by_ad = prediction.predict_ad_to_ad_obs(constants, fit_params)
    far_best_fits = _average_near_hall_predictions(far_best_fits_ad_by_ad)
    no_osc_params = fit_params.clone()
    no_osc_params.theta13 = 0
    no_osc_params.pull_theta12 = -1  # turn off theta12
    far_no_osc_ad_by_ad = prediction.predict_ad_to_ad_obs(constants, no_osc_params)
    far_no_osc = _average_near_hall_predictions(far_no_osc_ad_by_ad)
    data = constants.observed_candidates
    fig = _spectrum_ratio_plot(
        constants.reco_bins,
        [far_no_osc, far_best_fits],
        [(data, {halldet: np.sqrt(obs) for halldet, obs in data.items()})],
        [0.9, 1.02],
        ['No oscillations', 'Best fit'],
        ['Data'],
    )
    return fig


def plot_subtracted_prompt_spectrum(constants, fit_params):
    far_best_fits_ad_by_ad = prediction.predict_ad_to_ad_obs(constants, fit_params)
    pulled_bg = prediction.bg_with_pulls(constants, fit_params, halls='far')
    far_pred_ibds_ad_by_ad = {}
    total_bg = prediction.ad_dict(0, halls='far')
    for bg_type, bg_dict in pulled_bg.items():
        for halldet, bg_spec in bg_dict.items():
            total_bg[halldet] += bg_spec
    for (far_halldet, near_halldet), n_pred in far_best_fits_ad_by_ad.items():
        far_pred_ibds_ad_by_ad[far_halldet, near_halldet] = (
            n_pred - total_bg[far_halldet]
        )
    far_best_fits = _average_near_hall_predictions(far_pred_ibds_ad_by_ad)
    no_osc_params = fit_params.clone()
    no_osc_params.theta13 = 0
    no_osc_params.pull_theta12 = -1  # turn off theta12
    far_no_osc_ad_by_ad = prediction.predict_ad_to_ad_obs(constants, no_osc_params)
    far_no_osc_ibds_ad_by_ad = {}
    for (far_halldet, near_halldet), n_pred in far_no_osc_ad_by_ad.items():
        far_no_osc_ibds_ad_by_ad[far_halldet, near_halldet] = (
            n_pred - total_bg[far_halldet]
        )
    far_no_osc = _average_near_hall_predictions(far_no_osc_ibds_ad_by_ad)
    data_w_bg = constants.observed_candidates
    data = {}
    for halldet, n_obs in data_w_bg.items():
        if halldet in common.far_ads:
            data[halldet] = n_obs - total_bg[halldet]
    fig = _spectrum_ratio_plot(
        constants.reco_bins,
        [far_no_osc, far_best_fits],
        [(data, {halldet: np.sqrt(obs) for halldet, obs in data_w_bg.items()})],
        [0.9, 1.02],
        ['No oscillations', 'Best fit'],
        ['Data'],
    )
    return fig

def plot_fitted_points(constants, fit_params):
    """Plot the individual terms in the chi-square."""
    far_best_fits_ad_by_ad = prediction.predict_ad_to_ad_obs(constants, fit_params)
    far_best_fits = _average_near_hall_predictions(far_best_fits_ad_by_ad)
    data = constants.observed_candidates
    names = [f'EH{hall}\nAD{det}' for hall, det in prediction.far_ads]
    pretend_AD_bins = np.linspace(-0.5, 3.5, 5, endpoint=True)
    fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3]})
    ax = axs[0]
    _plot_line_hist(
        ax,
        pretend_AD_bins,
        [far_best_fits[halldet].sum() for halldet in prediction.far_ads],
        label='Best fit',
        endpoints=False,
        linestyle='--',
    )
    _plot_point_hist(
        ax,
        pretend_AD_bins,
        [data[halldet].sum() for halldet in prediction.far_ads],
        label='Data',
        yerr=[np.sqrt(data[halldet].sum()) for halldet in prediction.far_ads],
    )
    ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([0, 1, 2, 3]))
    ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(names))

    ax = axs[1]
    pull_contributions = fit.chi_square(
        constants, fit_params, return_array=True, rate_only=True, avg_near=True
    )[4:]
    ax.plot(pull_contributions, '.')
    ax.set_yscale('log')

    return fig



def main(database, fit_id):
    fit_params, constants, fit_info = fit.load_result(database, fit_id)
    #fig = plot_prompt_spectrum(constants, fit_params)
    #plt.show()
    #fig = plot_subtracted_prompt_spectrum(constants, fit_params)
    #plt.show()
    fig = plot_fitted_points(constants, fit_params)
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('fit_id', type=int)
    args = parser.parse_args()
    main(args.database, args.fit_id)
