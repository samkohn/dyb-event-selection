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

def _plot_shaded_errors_hist(ax, indexes, central_vals, errs, **kwargs):
    heights = 2 * np.array(errs)
    bottoms = np.array(central_vals) - errs
    width = 1
    color = '#D0D0D0'
    return ax.bar(
        indexes, heights, width, bottoms, align='center', color=color, **kwargs
    )

def _plot_point_hist(ax, bin_edges, vals, **kwargs):
    """Plot a "step" histogram with no vertical lines connecting the bins.

    Plotting algorithm based on https://stackoverflow.com/a/44961513/1144588.
    A NaN value is inserted between each value.
    Since NaNs are not plotted, the line jumps to the next "real" value.
    """
    bin_centers = bin_edges[:-1] + 0.5 * np.diff(bin_edges)
    X = np.c_[bin_edges[:-1], bin_centers, bin_edges[1:], bin_edges[1:]].flatten()
    Y = np.c_[vals, vals, vals, np.zeros_like(bin_edges[:-1])*np.nan].flatten()
    if 'marker' in kwargs:
        kwargs['markevery'] = (1, 4)
    if 'yerr' in kwargs:
        plain_err = kwargs.pop('yerr')
        err = np.c_[plain_err, plain_err, plain_err, np.zeros_like(bin_edges[:-1])].flatten()
        return ax.errorbar(X, Y, yerr=err, errorevery=(1, 4), **kwargs)
    else:
        return ax.plot(X, Y, **kwargs)

def _plot_pulls(ax, best_fit_pulls, errors, keys, labels, title):
    num_pulls = len(best_fit_pulls)
    indexes = np.arange(num_pulls)
    pretend_AD_bins = np.arange(num_pulls + 1) - 0.5
    _plot_shaded_errors_hist(
        ax,
        indexes,
        0,
        [errors[key] for key in keys],
        label=r'$1\sigma$ constraint',
    )
    _plot_line_hist(
        ax,
        pretend_AD_bins,
        np.zeros((num_pulls,)),
        label='Expected',
        endpoints=False,
        linestyle='--',
    )
    _plot_point_hist(
        ax,
        pretend_AD_bins,
        [best_fit_pulls[key] for key in keys],
        label='Best fit pull parameters',
        marker='o',
    )
    ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(np.arange(num_pulls)))
    ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(labels))
    if len(labels) > 6:
        ax.tick_params(axis='x', labelsize=8)
    ax.set_title(title, y=0.8)
    return


def _spectrum_ratio_plot(
    bin_edges,
    line_hist_dicts,
    errorbar_dicts,
    ratio_ylim,
    line_hist_labels,
    errorbar_labels,
):
    fig, axs_deep = plt.subplots(
        4, 2, gridspec_kw={'height_ratios': [3, 1, 3, 1], 'hspace': 0, 'wspace': 0},
        sharex=True, sharey='row'
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
        ax.grid()
        ax.tick_params(axis='both', which='major', labelsize=14)
        #ax.text(0.82, 0.95, name, fontsize=12, transform=ax.transAxes,
            #verticalalignment='top')
        if ax_index_main in (0, 2):
            ax.set_ylim([0, ax.get_ylim()[1]*1.1])
            ax.set_ylabel('Number of events', fontsize=16)
        if ax_index_main == 0:
            ax.legend(
                fontsize=12,
                loc='upper right',
                title=name,
                title_fontsize=14,
                frameon=False,
            )
        else:
            ax.legend(
                [],
                title=name,
                title_fontsize=14,
                loc='upper right',
                frameon=False,
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
        if ax_index_main in (0, 2):
            ax.set_ylim(ratio_ylim)
            ax.set_ylabel('Ratio to no osc.', fontsize=12)
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

def plot_data_fit_points(constants, fit_params):
    """Plot the "observation" and "prediction" used in the chi-square."""
    far_best_fits_ad_by_ad = prediction.predict_ad_to_ad_obs(constants, fit_params)
    far_best_fits = _average_near_hall_predictions(far_best_fits_ad_by_ad)
    data = constants.observed_candidates
    far_names = [f'EH{hall}\nAD{det}' for hall, det in prediction.far_ads]
    pretend_AD_bins = np.linspace(-0.5, 3.5, 5, endpoint=True)
    fig, ax = plt.subplots()
    _plot_point_hist(
        ax,
        pretend_AD_bins,
        [far_best_fits[halldet].sum() for halldet in prediction.far_ads],
        label='Best fit',
        linestyle='--',
        linewidth=3,
    )
    _plot_point_hist(
        ax,
        pretend_AD_bins,
        [data[halldet].sum() for halldet in prediction.far_ads],
        label='Data',
        yerr=[np.sqrt(data[halldet].sum()) for halldet in prediction.far_ads],
        linewidth=3,
    )
    ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([0, 1, 2, 3]))
    ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(far_names))
    ax.legend(fontsize=12)
    ax.set_ylabel('Number of events', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.grid(True)
    return fig

def _other_contrib_labels():
    """Generate a list of labels for all chi-square return array entries
    that aren't statistics.
    """
    acc = [f'Acc. EH{hall}-AD{det}' for hall, det in common.all_ads]
    li9 = [f'Li9 EH{hall}' for hall in (1, 2, 3)]
    fastn = [f'Fast neutrons EH{hall}' for hall in (1, 2, 3)]
    amc = ['Am-C']
    radn = [f'Rad. neutrons EH{hall}-AD{det}' for hall, det in common.all_ads]
    reactor = [f'Reactor power {core}' for core in ['D1', 'D2', 'L1', 'L2', 'L3', 'L4']]
    efficiency = [f'Efficiency EH{hall}-AD{det}' for hall, det in common.all_ads]
    rel_escale = [f'Rel. E scale EH{hall}-AD{det}' for hall, det in common.all_ads]
    osc = [r'$\theta_{12}$', r'$\Delta m^2_{21}$', r'$\Delta m^2_{32}$']
    return acc + li9 + fastn + amc + radn + reactor + efficiency + rel_escale + osc

#def _stats_contrib_labels(reco_bins):
    #"""Generate a list of labels for all statistics chi-square return array entries."""
    #far = [f'EH{hall}-AD{det}' for hall, det in prediction.far_ads]
    #near = [f'{reco_bins[i]}-{reco_bins[i+1]} MeV, EH{hall}-AD{det}'
            #for hall, det in prediction.near_ads
            #for i in range(len(reco_bins) - 1)
        #]
    #return far + near

def _far_stats_contrib_labels():
    """Generate a list of labels for the far hall statistics."""
    return [f'Data EH{hall}-AD{det}' for hall, det in common.far_ads]

def _near_stats_contrib_labels():
    """Generate a list of labels for the near hall statistics."""
    return [f'Near stats EH{hall}-AD{det}' for hall, det in common.near_ads]

def plot_chi2_contrib(constants, fit_params):
    """Plot each term's contributions to the chi-square,
    except for the near stat pulls.
    """
    rate_only = True
    avg_near = True
    contribs = fit.chi_square(
        constants,
        fit_params,
        return_array=True,
        rate_only=rate_only,
        avg_near=avg_near,
    )[:-8]  # hack because there are 8 extra 0's due to no alpha-n
    far_stat_contribs = contribs[:4]
    near_stat_contribs_binned = contribs[27 : 27 + 4 * 34]
    near_stat_contribs_total = np.empty((4,))
    for i, halldet in enumerate(common.near_ads):
        near_stat_contribs_total[i] = sum(near_stat_contribs_binned[i*34 : (i+1)*34])
    other_contribs = np.concatenate((contribs[4:27], contribs[27 + 4 * 34 :]))

    fig, ax = plt.subplots(figsize=(6, 7))
    ax.invert_yaxis()
    syst_handle = ax.barh(_other_contrib_labels(), other_contribs, log=True, height=1,
        linewidth=1, edgecolor='w')
    near_handle = ax.barh(_near_stats_contrib_labels(), near_stat_contribs_total,
        log=True, height=1, linewidth=1, edgecolor='w')
    far_handle = ax.barh(_far_stats_contrib_labels(), far_stat_contribs, log=True,
        height=1, linewidth=1, edgecolor='w')
    ax.grid(False)
    ax.grid(True, axis='x', color='w')
    ax.tick_params(axis='y', labelsize=6)
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', right=False, left=False)
    ax.tick_params(axis='x', which='minor', labelbottom=False, length=3, width=1)
    ax.set_ylim([
        len(other_contribs)
        + len(far_stat_contribs)
        + len(near_stat_contribs_total),
        -1,
    ])
    ax.set_xlim([3e-11, 2e-1])
    ax.set_xlabel(r'$\delta\chi^2$', fontsize=9)
    fig.tight_layout()
    ax.legend(
        [syst_handle, near_handle, far_handle],
        ['Systematics', 'Near hall stats', 'Data (Far hall stats)'],
        frameon=False,
        loc='upper right',
        bbox_to_anchor=(1, 0.85),
    )

    #axs[1].invert_yaxis()
    #axs[1].barh(_stats_contrib_labels(constants.reco_bins), stat_contribs, log=True,
        #height=1, linewidth=1, edgecolor='w')
    #axs[1].grid(False)
    #axs[1].grid(True, axis='x', color='w')
    #axs[1].xaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=6))
    #axs[1].xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, numticks=18))
    #axs[1].tick_params(axis='y', labelsize=6)
    #axs[1].tick_params(axis='x', labelsize=8)
    #axs[1].tick_params(axis='y', right=False, left=False)
    #axs[1].tick_params(axis='x', which='minor', labelbottom=False, length=3, width=1)
    #axs[1].set_ylim([-1, len(stat_contribs) + 1])
    #for label in axs[1].xaxis.get_ticklabels()[1::2]:
        #label.set_visible(False)
    fig.savefig('test_pulls.pdf', bbox_inches='tight')
    return fig

def plot_near_stat_pulls(constants, fit_params):
    rate_only = True
    avg_near = True
    contribs = fit.chi_square(
        constants,
        fit_params,
        return_array=True,
        rate_only=rate_only,
        avg_near=avg_near,
    )[:-8]  # hack because there are 8 extra 0's due to no alpha-n
    near_stat_dict = {}
    labels = [f'EH{hall}-AD{det}' for hall, det in common.near_ads]
    fig, ax = plt.subplots()
    for (halldet, near_stat_contrib), label in zip(near_stat_dict.items(), labels):
        _plot_point_hist(ax, constants.reco_bins, near_stat_contrib, label=label)
    ax.set_yscale('log')
    ax.set_xlabel('Prompt energy [MeV]')
    ax.set_ylabel(r'\delta\chi^2')
    ax.legend()
    return fig


def plot_pulls(constants, fit_params):
    """Plot the pull parameter best-fits."""
    far_names = [f'EH{hall}\nAD{det}' for hall, det in prediction.far_ads]
    all_names = [f'EH{hall}\nAD{det}' for hall, det in prediction.all_ads]
    fig, axs_deep = plt.subplots(4, 2)
    axs = axs_deep.flatten()
    axs[2].set_ylabel('Fractional variation', fontsize=14)
    _plot_pulls(
        axs[1],
        fit_params.pull_accidental,
        constants.bg_errors['accidental'],
        prediction.all_ads,
        all_names,
        'Accidentals'
    )
    _plot_pulls(
        axs[2],
        fit_params.pull_reactor,
        {i: constants.reactor_err for i in range(1, 7)},
        range(1, 7),
        [f'D{i}' if i <= 2 else f'L{i-2}' for i in range(1, 7)],
        'Reactor power',
    )
    _plot_pulls(
        axs[3],
        fit_params.pull_efficiency,
        {halldet: constants.efficiency_err for halldet in prediction.all_ads},
        prediction.all_ads,
        all_names,
        'Detection efficiency',
    )
    _plot_pulls(
        axs[4],
        fit_params.pull_rel_escale,
        {halldet: constants.rel_escale_err for halldet in prediction.all_ads},
        prediction.all_ads,
        all_names,
        'Relative energy scale',
    )
    li9_amc_pulls = fit_params.pull_li9.copy()
    li9_amc_pulls['amc'] = fit_params.pull_amc
    li9_amc_errors = {hall: constants.bg_errors['li9'][hall, 1] for hall in [1, 2, 3]}
    li9_amc_errors['amc'] = constants.bg_errors['amc'][1, 1]
    _plot_pulls(
        axs[5],
        li9_amc_pulls,
        li9_amc_errors,
        [1, 2, 3, 'amc'],
        [r'${}^9$Li ' f'EH{i}' for i in [1, 2, 3]] + ['AmC'],
        r'${}^9$Li/${}^8$He and AmC',
    )
    _plot_pulls(
        axs[6],
        fit_params.pull_fast_neutron,
        {hall: constants.bg_errors['fast-neutron'][hall, 1] for hall in [1, 2, 3]},
        [1, 2, 3],
        [f'EH{i}' for i in [1, 2, 3]],
        'Fast neutrons',
    )
    _plot_pulls(
        axs[7],
        {
            'theta12': fit_params.pull_theta12,
            'm2_21': fit_params.pull_m2_21,
            'm2_ee': fit_params.pull_m2_ee,
        },
        {
            'theta12': constants.theta12_err,
            'm2_21': constants.m2_21_err,
            'm2_ee': constants.m2_ee_err,
        },
        ['theta12', 'm2_21', 'm2_ee'],
        [r'$\theta_{12}$', r'$\Delta m^2_{21}$', r'$\Delta m^2_{ee}$'],
        'Input oscillation parameters',
    )



    return fig


def main(database, fit_id):
    fit_params, constants, fit_info = fit.load_result(database, fit_id)
    #fig_prompt = plot_prompt_spectrum(constants, fit_params)
    #fig_prompt_sub = plot_subtracted_prompt_spectrum(constants, fit_params)
    fig_rate_only_data = plot_data_fit_points(constants, fit_params)
    fig_rate_only_data.tight_layout()
    #fig_pulls = plot_pulls(constants, fit_params)
    #fig_contribs = plot_chi2_contrib(constants, fit_params)
    #fig_near_contrib = plot_near_stat_pulls(constants, fit_params)
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('fit_id', type=int)
    args = parser.parse_args()
    main(args.database, args.fit_id)
