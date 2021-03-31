import argparse
from datetime import datetime
import sqlite3

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.ndimage.filters import uniform_filter1d
import numpy as np

from plot_results import _plot_point_hist


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('--label')
    parser.add_argument('--rates-label')
    parser.add_argument('--fig-titles', action='store_true',
        help='Include a title in the figure')
    choices = ('all', 'singles', 'muons', 'delayed-fit',
            'delayed-eff-unc', 'DT-eff', 'acc-DT-eff',
            'singles-within-run')
    parser.add_argument('--types', nargs='+', choices=choices, default='all',
            metavar='PLOT_TYPE', help=f'Choices: {choices}')
    args = parser.parse_args()
    database = args.database
    plot_types = args.types
    plot_all = 'all' in plot_types
    mpl.rcParams['font.size'] = 14
    mpl.rcParams['figure.figsize'] = (5, 5)
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['lines.markersize'] = 5
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['axes.labelsize'] = 14
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['xtick.major.size'] = 6
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['ytick.major.size'] = 6
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rcParams['legend.frameon'] = True
    mpl.rcParams['legend.fontsize'] = 12
    mpl.rcParams['legend.markerscale'] = 2
    wide_fig = (16, 8)
    big_fig = (16, 10)
    colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
    dateconv = np.vectorize(datetime.fromtimestamp)
    def get_dates(timestamp_ns_array):
        return dateconv(timestamp_ns_array/1e9)
    def get_ad_data(hall, det, in_data):
        return in_data[(in_data[:, 1] == hall) & (in_data[:, 2] == det)]
    ads = [(1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (3, 4)]
    near_ads = ads[:4]
    far_ads = ads[4:]
    def name(hall, det):
        return f'EH{hall}-AD{det}'
    ad_names = [name(hall, det) for hall, det in ads]
    ad_names_2line = [name(hall, det).replace('-', '\n') for hall, det
            in ads]
    near_names = [name(hall, det) for hall, det in near_ads]
    far_names = [name(hall, det) for hall, det in far_ads]

    def plot_ad_points(ax, y, **kwargs):
        """Plot the specified y values (and yerr if specified).

        Uses the "point histogram" style where the y-error, if any,
        is a vertical line with no "cap,"
        and the "x error bar" represents the bin width.
        """
        pretend_AD_bins = np.arange(9) - 0.5
        _plot_point_hist(
            ax,
            pretend_AD_bins,
            y,
            **kwargs,
        )
        ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(np.arange(8)))
        ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(ad_names_2line))
        return

    if plot_all or ('singles' in plot_types):
        print('Retrieving and plotting singles data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT
                    RunNo,
                    Hall,
                    DetNo,
                    Start_time,
                    Rate_Hz,
                    Rate_Hz_error
                FROM
                    singles_rates
                INNER JOIN
                    runs
                USING (RunNo)
                WHERE Label = ?
                ORDER BY RunNo
                ''',
                (args.rates_label,)
            )
            data = np.array(cursor.fetchall())
        # All
        fig, ax = plt.subplots()
        for hall, det in ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot(ad_data[:, 0], ad_data[:, 4], '.')
        ax.legend(ad_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('All ADs')
        ax.set_xlabel('Run number')
        ax.set_ylabel('Singles rate [Hz]')
        fig.tight_layout()
        fig.savefig('singles_all_byrun.pdf')

        fig, ax = plt.subplots()
        for hall, det in ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot_date(get_dates(ad_data[:, 3]), ad_data[:, 4], '.')
        ax.legend(ad_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('All ADs')
        ax.set_xlabel('Start time')
        ax.set_ylabel('Singles rate [Hz]')
        fig.tight_layout()
        fig.savefig('singles_all_bydate.pdf')


        # EH1, EH2
        fig, ax = plt.subplots()
        for hall, det in near_ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot(ad_data[:, 0], ad_data[:, 4], '.')
        ax.legend(near_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('Near halls')
        ax.set_xlabel('Run number')
        ax.set_ylabel('Singles rate [Hz]')
        fig.tight_layout()
        fig.savefig('singles_near_byrun.pdf')

        fig, ax = plt.subplots()
        for hall, det in near_ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot_date(get_dates(ad_data[:, 3]), ad_data[:, 4], '.')
        ax.legend(near_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('Near halls')
        ax.set_xlabel('Start time')
        ax.set_ylabel('Singles rate [Hz]')
        fig.tight_layout()
        fig.savefig('singles_near_bydate.pdf')

        #EH3
        fig, ax = plt.subplots()
        for hall, det in far_ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot(ad_data[:, 0], ad_data[:, 4], '.')
        ax.legend(far_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('Far hall')
        ax.set_xlabel('Run number')
        ax.set_ylabel('Singles rate [Hz]')
        fig.tight_layout()
        fig.savefig('singles_far_byrun.pdf')

        fig, ax = plt.subplots()
        for hall, det in far_ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot_date(get_dates(ad_data[:, 3]), ad_data[:, 4], '.')
        ax.legend(far_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('Far hall')
        ax.set_xlabel('Start time')
        ax.set_ylabel('Singles rate [Hz]')
        fig.tight_layout()
        fig.savefig('singles_far_bydate.pdf')

    if plot_all or ('muons' in plot_types):
        print('Retrieving and plotting muon data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT
                    RunNo,
                    Hall,
                    DetNo,
                    Start_time,
                    Efficiency
                FROM
                    muon_rates
                INNER JOIN
                    runs
                USING (RunNo)
                WHERE Label = ?
                ORDER BY RunNo
                ''',
                (args.rates_label,)
            )
            data = np.array(cursor.fetchall())

        #EH1, EH2
        fig, ax = plt.subplots()
        for hall, det in near_ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot(ad_data[:, 0], ad_data[:, 4], '.')
        ax.legend(near_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('Near halls')
        ax.set_xlabel('Run number')
        ax.set_ylabel('Muon efficiency')
        fig.tight_layout()
        fig.savefig('muon_eff_near_byrun.pdf')

        fig, ax = plt.subplots()
        for hall, det in near_ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot_date(get_dates(ad_data[:, 3]), ad_data[:, 4], '.')
        ax.legend(near_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('Near halls')
        ax.set_xlabel('Start time')
        ax.set_ylabel('Muon efficiency')
        fig.tight_layout()
        fig.savefig('muon_eff_near_bydate.pdf')

        #EH3
        fig, ax = plt.subplots()
        for hall, det in far_ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot(ad_data[:, 0], ad_data[:, 4], '.')
        ax.legend(far_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('Far hall')
        ax.set_xlabel('Run number')
        ax.set_ylabel('Muon efficiency')
        fig.tight_layout()
        fig.savefig('muon_eff_far_byrun.pdf')

        fig, ax = plt.subplots()
        for hall, det in far_ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot_date(get_dates(ad_data[:, 3]), ad_data[:, 4], '.')
        ax.legend(far_names, fontsize=12)
        if args.fig_titles:
            ax.set_title('Far hall')
        ax.set_xlabel('Start time')
        ax.set_ylabel('Muon efficiency')
        fig.tight_layout()
        fig.savefig('muon_eff_far_bydate.pdf')

    if plot_all or ('delayed-fit' in plot_types):
        print('Retrieving and plotting delayed energy fit data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT
                    Hall, DetNo, Peak, Peak_error,
                    Resolution, Resolution_error,
                    ExpScale, ExpScale_error,
                    PeakFraction, PeakFraction_error,
                    Normalization, Normalization_error
                FROM delayed_energy_fits
                WHERE
                    Source = ?
                ORDER BY
                    Hall,
                    DetNo
                ''',
                (args.label,)
            )
            data = np.array(cursor.fetchall())

        # Peak locations
        fig = plt.figure()
        spacing = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0)
        ax1 = plt.subplot(spacing[0])
        plot_ad_points(ax1, data[:, 2], yerr=data[:, 3])
        if args.fig_titles:
            ax1.set_title('Delayed energy peak & fit error')
        ax1.set_ylabel('Delayed energy peak value [MeV]')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax2 = plt.subplot(spacing[1], sharex=ax1)
        near_hall_avg = np.average(data[:4, 2], weights=data[:4, 3]**(-2))
        rel_deviations = (data[:, 2] - near_hall_avg)/near_hall_avg
        rel_errors = data[:, 3]/near_hall_avg
        plot_ad_points(ax2, rel_deviations, yerr=rel_errors)
        ax2.set_ylim([-0.007, 0.007])
        fig.tight_layout()
        fig.savefig('delayed_energy_peak.pdf')

        # Resolution
        fig = plt.figure()
        spacing = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0)
        ax1 = plt.subplot(spacing[0])
        plot_ad_points(ax1, data[:, 4], yerr=data[:, 5])
        if args.fig_titles:
            ax1.set_title('Delayed energy "resolution/width" & fit error')
        ax1.set_ylabel('Delayed energy peak width [MeV]')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax2 = plt.subplot(spacing[1], sharex=ax1)
        near_hall_avg = np.average(data[:4, 4], weights=data[:4, 5]**(-2))
        rel_deviations = (data[:, 4] - near_hall_avg)/near_hall_avg
        rel_errors = data[:, 5]/near_hall_avg
        plot_ad_points(ax2, rel_deviations, yerr=rel_errors)
        ax2.set_ylim([-0.025, 0.025])
        fig.tight_layout()
        fig.savefig('delayed_energy_width.pdf')


        # Exponential scale
        fig = plt.figure()
        spacing = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0)
        ax1 = plt.subplot(spacing[0])
        plot_ad_points(ax1, data[:, 6], yerr=data[:, 7])
        if args.fig_titles:
            ax1.set_title('Delayed energy tail slope & fit error')
        ax1.set_ylabel('Delayed energy tail slope [1/MeV]')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax2 = plt.subplot(spacing[1], sharex=ax1)
        near_hall_avg = np.average(data[:4, 6], weights=data[:4, 7]**(-2))
        rel_deviations = (data[:, 6] - near_hall_avg)/near_hall_avg
        rel_errors = data[:, 7]/near_hall_avg
        plot_ad_points(ax2, rel_deviations, yerr=rel_errors)
        fig.tight_layout()
        fig.savefig('delayed_energy_expo_scale.pdf')

        # Peak fraction
        fig = plt.figure()
        spacing = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0)
        ax1 = plt.subplot(spacing[0])
        plot_ad_points(ax1, data[:, 8], yerr=data[:, 9])
        if args.fig_titles:
            ax1.set_title('Delayed energy peak fraction & fit error')
        ax1.set_ylabel('Delayed energy peak fraction')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax2 = plt.subplot(spacing[1], sharex=ax1)
        near_hall_avg = np.average(data[:4, 8], weights=data[:4, 9]**(-2))
        rel_deviations = (data[:, 8] - near_hall_avg)/near_hall_avg
        rel_errors = data[:, 9]/near_hall_avg
        plot_ad_points(ax2, rel_deviations, yerr=rel_errors)
        ax2.set_ylim([-0.035, 0.035])
        fig.tight_layout()
        fig.savefig('delayed_energy_peak_frac.pdf')

        # Upper and lower fit bounds
        fig, axs = plt.subplots(2, 1, sharex=True)
        plot_ad_points(
            axs[0],
            data[:, 2] + 3 * data[:, 4],
            yerr=np.hypot(data[:, 3], data[:, 5]),
        )
        axs[0].set_ylabel('Upper bound [MeV]')
        axs[0].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.005))
        plot_ad_points(
            axs[1],
            data[:, 2] - 3 * data[:, 4],
            yerr=np.hypot(data[:, 3], data[:, 5]),
        )
        axs[1].set_ylabel('Lower bound [MeV]')
        axs[1].yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.005))
        fig.tight_layout()
        fig.savefig('delayed_energy_bounds.pdf')

    if plot_all or ('delayed-eff-unc' in plot_types):
        print('Retrieving and plotting delayed energy eff. uncertainty data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT
                    Hall,
                    DetNo,
                    RelativeDeviation,
                    StatError
                FROM
                    delayed_energy_uncertainty_1
                WHERE
                    Label = ?
                ORDER BY
                    Hall,
                    DetNo
                ''',
                (args.label,)
            )
            data = np.array(cursor.fetchall())

        fig, ax = plt.subplots()
        plot_ad_points(ax, data[:, 2], yerr=data[:, 3])
        if args.fig_titles:
            ax.set_title('Delayed energy efficiency uncertainty')
        ax.set_ylabel('Relative deviation from fit model')
        fig.tight_layout()
        fig.savefig('delayed_energy_uncertainty_method1.pdf')

    if plot_all or ('DT-eff' in plot_types):
        print('Retrieving and plotting distance-time (DT) efficiency data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT
                    Hall,
                    DetNo,
                    Efficiency,
                    StatError
                FROM
                    distance_time_cut_efficiency
                WHERE
                    Source = ?
                ORDER BY
                    Hall,
                    DetNo
                ''',
                (args.label,)
            )
            data = np.array(cursor.fetchall())

        fig, ax = plt.subplots()
        plot_ad_points(ax, data[:, 2], yerr=data[:, 3])
        if args.fig_titles:
            ax.set_title('Distance-time (DT) cut efficiency')
        ax.set_ylabel('Efficiency')
        ax.grid(False, axis='x')
        fig.tight_layout()
        fig.savefig('distance_time_cut_efficiency.pdf')

        #raise RuntimeError("Double check the Source in DB query below!")
        #with sqlite3.Connection(database) as conn:
            #cursor = conn.cursor()
            #cursor.execute('''SELECT Hall, DetNo, Efficiency, StatError
                #FROM distance_time_cut_efficiency
                #WHERE Source = "adtime"
                #ORDER BY Hall, DetNo''')
            #data = np.array(cursor.fetchall())

        #fig, ax = plt.subplots()
        #ax.errorbar(ad_names_2line, data[:, 2], yerr=data[:, 3], fmt='o')
        #ax.set_title('Distance-time (DT) cut efficiency')
        #ax.set_ylabel('Efficiency')
        #ax.grid()
        #fig.tight_layout()
        #fig.savefig('distance_time_cut_efficiency_adtime.pdf')

    if 'acc-DT-eff' in plot_types:  # don't count in "all"
        print('Retrieving and plotting accidentals DT efficiency data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT RunNo, Hall, DetNo, Start_time, Livetime_ns,
                DistanceTime_DT_Eff
                FROM runs NATURAL JOIN accidental_subtraction NATURAL JOIN muon_rates
                ORDER BY RunNo, DetNo''')
            data = np.array(cursor.fetchall())
            cursor.execute('''SELECT RunNo, Hall, DetNo, DistanceTime_DT_Eff,
                DistanceTime_DT_Eff_error
                FROM runs NATURAL JOIN distance_time_eff_study
                WHERE PairingType = "random_N"
                ORDER BY RunNo, DetNo''')
            extra_data = np.array(cursor.fetchall())
            cursor.execute('''SELECT RunNo, Hall, DetNo, Start_time, Livetime_ns, DistanceTime_DT_Eff,
                DistanceTime_DT_Eff_error
                FROM runs NATURAL JOIN distance_time_eff_study NATURAL JOIN muon_rates
                WHERE PairingType = "random_many"
                ORDER BY RunNo, DetNo''')
            data_random_many = np.array(cursor.fetchall())
        ns_per_18h = 60e9 * 60 * 18

        #EH1, EH2
        fig, ax = plt.subplots(figsize=wide_fig)
        for hall, det in near_ads:
            ad_data = get_ad_data(hall, det, data)
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            ax.plot(long_runs[:, 0], long_runs[:, 5], '.')
        ax.legend(near_names, fontsize=12)
        ax.set_title('Near halls, runtime >18h')
        ax.set_xlabel('Run number')
        ax.set_ylabel(r'$\varepsilon_{DT,acc}$')
        fig.tight_layout()
        fig.savefig('acc_DT_eff_near_byrun.pdf')

        fig, ax = plt.subplots(figsize=wide_fig)
        for i, (hall, det) in enumerate(near_ads):
            ad_data = get_ad_data(hall, det, data)
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            average = uniform_filter1d(long_runs[:, 5], size=20, mode='nearest')
            ax.plot(long_runs[10:-10:20, 0], average[10:-10:20], '.')
        ax.legend(near_names, fontsize=12)
        ax.set_title('Near halls, runtime >18h, average over $\pm10$ runs')
        ax.set_xlabel('Run number')
        ax.set_ylabel(r'$\varepsilon_{DT,acc}$')
        fig.tight_layout()
        fig.savefig('acc_DT_eff_avg_near_byrun.pdf')

        fig, ax = plt.subplots(figsize=wide_fig)
        for hall, det in near_ads:
            ad_data = get_ad_data(hall, det, data)
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            ax.plot_date(get_dates(long_runs[:, 3]), long_runs[:, 5], '.')
        ax.legend(near_names, fontsize=12)
        ax.set_title('Near halls, runtime >18h')
        ax.set_xlabel('Start time')
        ax.set_ylabel(r'$\varepsilon_{DT,acc}$')
        fig.tight_layout()
        fig.savefig('acc_DT_eff_near_bydate.pdf')

        #EH3
        fig, ax = plt.subplots(figsize=wide_fig)
        for hall, det in far_ads:
            ad_data = get_ad_data(hall, det, data)
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            ax.plot(long_runs[:, 0], long_runs[:, 5], '.')
        ax.legend(far_names, fontsize=12)
        ax.set_title('Far halls, runtime >18h')
        ax.set_xlabel('Run number')
        ax.set_ylabel(r'$\varepsilon_{DT,acc}$')
        fig.tight_layout()
        fig.savefig('acc_DT_eff_far_byrun.pdf')

        fig, ax = plt.subplots(figsize=wide_fig)
        for hall, det in far_ads:
            ad_data = get_ad_data(hall, det, data)
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            average = uniform_filter1d(long_runs[:, 5], size=20, mode='nearest')
            ax.plot(long_runs[10:-10:20, 0], average[10:-10:20], '.')
        ax.legend(far_names, fontsize=12)
        ax.set_title('Far halls, runtime >18h, average over $\pm10$ runs')
        ax.set_xlabel('Run number')
        ax.set_ylabel(r'$\varepsilon_{DT,acc}$')
        fig.tight_layout()
        fig.savefig('acc_DT_eff_avg_far_byrun.pdf')

        fig, ax = plt.subplots(figsize=wide_fig)
        for hall, det in far_ads:
            ad_data = get_ad_data(hall, det, data)
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            ax.plot_date(get_dates(long_runs[:, 3]), long_runs[:, 5], '.')
        ax.legend(far_names, fontsize=12)
        ax.set_title('Far halls, runtime >18h')
        ax.set_xlabel('Start time')
        ax.set_ylabel(r'$\varepsilon_{DT,acc}$')
        fig.tight_layout()
        fig.savefig('acc_DT_eff_far_bydate.pdf')

    if 'singles-within-run' in plot_types:
        print('Retrieving and plotting accidentals DT efficiency data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT RunNo, Hall, DetNo, Start_time, Livetime_ns,
                NumFirstHour, NumLastHour
                FROM runs NATURAL JOIN muon_rates NATURAL JOIN singles_within_run
                WHERE Livetime_ns > 64800000000000
                ORDER BY RunNo, DetNo''')
            data = np.array(cursor.fetchall())
                    # AND NumFirstHour > 55000 AND NumLastHour > 55000

        fig, axs = plt.subplots(2, 2, figsize=big_fig, sharex=True)
        for ax, name, (hall, det) in zip(axs.flat, near_names, near_ads):
            ad_data = get_ad_data(hall, det, data)
            diffs = (ad_data[:, 6] - ad_data[:, 5])/ad_data[:, 5]
            ax.hist(diffs[np.abs(diffs) < 0.1], histtype='step',
                    linewidth=2)
            ax.legend([name], fontsize=12)
            ax.grid()
        axs[0, 0].set_title('Near halls, runtime > 18h')
        ax.set_xlabel('Fractional increase from first to last hour', fontsize=14)
        ax.set_ylabel('Number of runs')
        fig.tight_layout()
        fig.savefig('singles_within_run_near.pdf')

        fig, axs = plt.subplots(2, 2, figsize=big_fig, sharex=True)
        for ax, name, (hall, det) in zip(axs.flat, far_names, far_ads):
            ad_data = get_ad_data(hall, det, data)
            diffs = (ad_data[:, 6] - ad_data[:, 5])/ad_data[:, 5]
            ax.hist(diffs[np.abs(diffs) < 0.1], histtype='step',
                    linewidth=2)
            ax.legend([name], fontsize=12)
            ax.grid()
        axs[0, 0].set_title('Far halls, runtime > 18h')
        ax.set_xlabel('Fractional increase from first to last hour', fontsize=14)
        ax.set_ylabel('Number of runs')
        fig.tight_layout()
        fig.savefig('singles_within_run_far.pdf')

        fig, axs = plt.subplots(2, 2, figsize=big_fig, sharex=True)
        for ax, name, (hall, det) in zip(axs.flat, near_names, near_ads):
            ad_data = get_ad_data(hall, det, data)
            ax.plot(ad_data[:, 5], ad_data[:, 6], '.')
            ax.legend([name], fontsize=12)
            ax.grid()
        axs[0, 0].set_title('Near halls, runtime > 18h')
        axs[1, 0].set_ylabel('Number of singles in last hour')
        ax.set_xlabel('Number of singles in first hour')
        fig.tight_layout()
        fig.savefig('singles_within_run_near_corr.pdf')

        fig, axs = plt.subplots(2, 2, figsize=big_fig, sharex=True)
        for ax, name, (hall, det) in zip(axs.flat, far_names, far_ads):
            ad_data = get_ad_data(hall, det, data)
            ax.plot(ad_data[:, 5], ad_data[:, 6], '.')
            ax.legend([name], fontsize=12)
            ax.grid()
        axs[0, 0].set_title('Far halls, runtime > 18h')
        axs[1, 0].set_ylabel('Number of singles in last hour')
        ax.set_xlabel('Number of singles in first hour')
        fig.tight_layout()
        fig.savefig('singles_within_run_far_corr.pdf')
