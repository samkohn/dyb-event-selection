"""Plot the fit results including prompt spectrum, L/E, backgrounds, etc."""

import argparse
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib.pyplot as plt

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

def _plot_line_hist(ax, bin_edges, vals, **kwargs):
    vals = np.concatenate(([0], vals, [0]))
    bin_edges = np.concatenate(([bin_edges[0]], bin_edges))
    ax.plot(bin_edges, vals, drawstyle='steps-post')
    return


def plot_prompt_spectrum(constants, fit_params):
    far_best_fits_ad_by_ad = prediction.predict_ad_to_ad_obs(constants, fit_params)
    far_best_fits = _average_near_hall_predictions(far_best_fits_ad_by_ad)
    no_osc_params = fit_params.clone()
    no_osc_params.theta13 = 0
    no_osc_params.pull_theta12 = -1  # turn off theta12
    far_no_osc_ad_by_ad = prediction.predict_ad_to_ad_obs(constants, no_osc_params)
    far_no_osc = _average_near_hall_predictions(far_no_osc_ad_by_ad)
    data = constants.observed_candidates
    fig, axs_deep = plt.subplots(2, 2)
    axs_flat = axs_deep.flatten()
    reco_bin_centers = constants.reco_bins[:-1] + 0.5 * np.diff(constants.reco_bins)
    for halldet, ax in zip(common.far_ads, axs_flat):
        name = f'EH{halldet[0]}-AD{halldet[1]}'
        #ax.plot(reco_bin_centers, far_no_osc[halldet], '.')
        _plot_line_hist(ax, constants.reco_bins, far_no_osc[halldet])
        _plot_line_hist(ax, constants.reco_bins, far_best_fits[halldet])
        obs = data[halldet]
        err = np.sqrt(obs)
        ax.errorbar(reco_bin_centers, obs, yerr=err, fmt='_', elinewidth=2, markersize=5)
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
    fig, axs_deep = plt.subplots(2, 2)
    axs_flat = axs_deep.flatten()
    reco_bin_centers = constants.reco_bins[:-1] + 0.5 * np.diff(constants.reco_bins)
    for halldet, ax in zip(common.far_ads, axs_flat):
        name = f'EH{halldet[0]}-AD{halldet[1]}'
        #ax.plot(reco_bin_centers, far_no_osc[halldet], '.')
        _plot_line_hist(ax, constants.reco_bins, far_no_osc[halldet])
        _plot_line_hist(ax, constants.reco_bins, far_best_fits[halldet])
        obs = data[halldet]
        err = np.sqrt(data_w_bg[halldet])
        ax.errorbar(reco_bin_centers, obs, yerr=err, fmt='_', elinewidth=2, markersize=5)
    return fig


def main(database, fit_id):
    fit_params, constants, fit_info = fit.load_result(database, fit_id)
    fig = plot_prompt_spectrum(constants, fit_params)
    plt.show()
    fig = plot_subtracted_prompt_spectrum(constants, fit_params)
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('fit_id', type=int)
    args = parser.parse_args()
    main(args.database, args.fit_id)
