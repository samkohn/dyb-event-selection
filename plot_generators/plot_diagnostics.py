import argparse
from datetime import datetime
import sqlite3

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.ndimage.filters import uniform_filter1d
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    choices = ('all', 'singles', 'muons', 'delayed-fit',
            'delayed-eff-unc', 'DT-eff', 'acc-DT-eff')
    parser.add_argument('--types', nargs='+', choices=choices, default='all',
            metavar='PLOT_TYPE', help=f'Choices: {choices}')
    args = parser.parse_args()
    database = args.database
    plot_types = args.types
    plot_all = 'all' in plot_types
    mpl.rcParams['font.size'] = 18
    mpl.rcParams['figure.figsize'] = (8, 8)
    mpl.rcParams['lines.linewidth'] = 3
    mpl.rcParams['lines.markersize'] = 10
    mpl.rcParams['errorbar.capsize'] = 5
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

    if plot_all or ('singles' in plot_types):
        print('Retrieving and plotting singles data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT RunNo, Hall, DetNo, Start_time, Rate_Hz, Rate_Hz_error
                FROM singles_rates INNER JOIN runs USING (RunNo)
                ORDER BY RunNo''')
            data = np.array(cursor.fetchall())
        # All
        fig, ax = plt.subplots()
        for hall, det in ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot(ad_data[:, 0], ad_data[:, 4], '.')
        ax.legend(ad_names, fontsize=12)
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
        ax.set_title('Far hall')
        ax.set_xlabel('Start time')
        ax.set_ylabel('Singles rate [Hz]')
        fig.tight_layout()
        fig.savefig('singles_far_bydate.pdf')

    if plot_all or ('muons' in plot_types):
        print('Retrieving and plotting muon data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT RunNo, Hall, DetNo, Start_time, Efficiency
                FROM muon_rates INNER JOIN runs USING (RunNo)
                ORDER BY RunNo''')
            data = np.array(cursor.fetchall())

        #EH1, EH2
        fig, ax = plt.subplots()
        for hall, det in near_ads:
            ad_data = get_ad_data(hall, det, data)
            ax.plot(ad_data[:, 0], ad_data[:, 4], '.')
        ax.legend(near_names, fontsize=12)
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
        ax.set_title('Far hall')
        ax.set_xlabel('Start time')
        ax.set_ylabel('Muon efficiency')
        fig.tight_layout()
        fig.savefig('muon_eff_far_bydate.pdf')

    if plot_all or ('delayed-fit' in plot_types):
        print('Retrieving and plotting delayed energy fit data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT Hall, DetNo, Peak, Peak_error,
            Resolution, Resolution_error,
            ExpScale, ExpScale_error,
            PeakFraction, PeakFraction_error,
            Normalization, Normalization_error
                FROM delayed_energy_fits ORDER BY Hall, DetNo''')
            data = np.array(cursor.fetchall())

        # Peak locations
        fig = plt.figure()
        spacing = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax1 = plt.subplot(spacing[0])
        ax1.errorbar(ad_names_2line, data[:, 2], yerr=data[:, 3], fmt='o')
        ax1.set_title('Delayed energy peak & fit error')
        ax1.set_ylabel('Delayed energy peak value [MeV]')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.grid()
        ax2 = plt.subplot(spacing[1], sharex=ax1)
        near_hall_avg = np.average(data[:4, 2], weights=data[:4, 3]**(-2))
        rel_deviations = (data[:, 2] - near_hall_avg)/near_hall_avg
        rel_errors = data[:, 3]/near_hall_avg
        ax2.errorbar(ad_names_2line, rel_deviations, yerr=rel_errors,
                fmt='.')
        ax2.set_ylim([-0.007, 0.007])
        ax2.grid()
        fig.tight_layout()
        fig.savefig('delayed_energy_peak.pdf')

        # Resolution
        fig = plt.figure()
        spacing = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax1 = plt.subplot(spacing[0])
        ax1.errorbar(ad_names_2line, data[:, 4], yerr=data[:, 5], fmt='o')
        ax1.set_title('Delayed energy "resolution/width" & fit error')
        ax1.set_ylabel('Delayed energy peak width [MeV]')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.grid()
        ax2 = plt.subplot(spacing[1], sharex=ax1)
        near_hall_avg = np.average(data[:4, 4], weights=data[:4, 5]**(-2))
        rel_deviations = (data[:, 4] - near_hall_avg)/near_hall_avg
        rel_errors = data[:, 5]/near_hall_avg
        ax2.errorbar(ad_names_2line, rel_deviations, yerr=rel_errors,
                fmt='.')
        ax2.set_ylim([-0.025, 0.025])
        ax2.grid()
        fig.tight_layout()
        fig.savefig('delayed_energy_width.pdf')


        # Exponential scale
        fig = plt.figure()
        spacing = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax1 = plt.subplot(spacing[0])
        ax1.errorbar(ad_names_2line, data[:, 6], yerr=data[:, 7], fmt='o')
        ax1.set_title('Delayed energy tail slope & fit error')
        ax1.set_ylabel('Delayed energy tail slope [1/MeV]')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.grid()
        ax2 = plt.subplot(spacing[1], sharex=ax1)
        near_hall_avg = np.average(data[:4, 6], weights=data[:4, 7]**(-2))
        rel_deviations = (data[:, 6] - near_hall_avg)/near_hall_avg
        rel_errors = data[:, 7]/near_hall_avg
        ax2.errorbar(ad_names_2line, rel_deviations, yerr=rel_errors,
                fmt='.')
        ax2.grid()
        fig.tight_layout()
        fig.savefig('delayed_energy_expo_scale.pdf')

        # Peak fraction
        fig = plt.figure()
        spacing = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        ax1 = plt.subplot(spacing[0])
        ax1.errorbar(ad_names_2line, data[:, 8], yerr=data[:, 9], fmt='o')
        ax1.set_title('Delayed energy peak fraction & fit error')
        ax1.set_ylabel('Delayed energy peak fraction')
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.grid()
        ax2 = plt.subplot(spacing[1], sharex=ax1)
        near_hall_avg = np.average(data[:4, 8], weights=data[:4, 9]**(-2))
        rel_deviations = (data[:, 8] - near_hall_avg)/near_hall_avg
        rel_errors = data[:, 9]/near_hall_avg
        ax2.errorbar(ad_names_2line, rel_deviations, yerr=rel_errors,
                fmt='.')
        ax2.set_ylim([-0.035, 0.035])
        ax2.grid()
        fig.tight_layout()
        fig.savefig('delayed_energy_peak_frac.pdf')

        # Upper and lower fit bounds
        fig, axs = plt.subplots(2, 1, sharex=True, figsize=(9, 6.8))
        axs[0].errorbar(ad_names_2line, data[:, 2] + 3 * data[:, 4],
                yerr=np.hypot(data[:, 3], data[:, 5]), fmt='o')
        axs[0].set_ylabel('Upper bound [MeV]')
        axs[0].grid()
        axs[1].errorbar(ad_names_2line, data[:, 2] - 3 * data[:, 4],
                yerr=np.hypot(data[:, 3], data[:, 5]), fmt='o')
        axs[1].set_ylabel('Lower bound [MeV]')
        axs[1].grid()
        fig.tight_layout()
        fig.savefig('delayed_energy_bounds.pdf')

    if plot_all or ('delayed-eff-unc' in plot_types):
        print('Retrieving and plotting delayed energy eff. uncertainty data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT Hall, DetNo, RelativeDeviation, StatError
                FROM delayed_energy_uncertainty_1 ORDER BY Hall, DetNo''')
            data = np.array(cursor.fetchall())

        fig, ax = plt.subplots(figsize=(9, 5.8))
        ax.errorbar(ad_names_2line, data[:, 2], yerr=data[:, 3], fmt='o')
        ax.set_title('Delayed energy efficiency uncertainty')
        ax.set_ylabel('Relative deviation from fit model')
        ax.grid()
        fig.tight_layout()
        fig.savefig('delayed_energy_uncertainty_method1.pdf')

    if plot_all or ('DT-eff' in plot_types):
        print('Retrieving and plotting distance-time (DT) efficiency data...')
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT Hall, DetNo, Efficiency, StatError
                FROM distance_time_cut_efficiency ORDER BY Hall, DetNo''')
            data = np.array(cursor.fetchall())

        fig, ax = plt.subplots()
        ax.errorbar(ad_names_2line, data[:, 2], yerr=data[:, 3], fmt='o')
        ax.set_title('Distance-time (DT) cut efficiency')
        ax.set_ylabel('Efficiency')
        ax.grid()
        fig.tight_layout()
        fig.savefig('distance_time_cut_efficiency.pdf')

    if plot_all or ('acc-DT-eff' in plot_types):
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
        wide_fig = (16, 8)
        big_fig = (16, 10)

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

        fig, axs = plt.subplots(4, 1, figsize=wide_fig, sharex='col')
        for ax, (hall, det) in zip(axs, near_ads):
            ad_data = get_ad_data(hall, det, data)
            extra_ad_data = get_ad_data(hall, det, extra_data)
            ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            average = uniform_filter1d(long_runs[:, 5], size=20, mode='nearest')
            average_method2 = uniform_filter1d(long_runs[:, 6], size=20, mode='nearest')
            ax.plot_date(get_dates(long_runs[10:-10:20, 3]), average[10:-10:20], '.')
            ax.errorbar(get_dates(long_runs[10:-10:20, 3]), average_method2[10:-10:20],
                    fmt='.', yerr=long_runs[10:-10:20, 7]/np.sqrt(20))
            ax.grid(axis='y')
        for ax, name in zip(axs, near_names):
            ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                    verticalalignment='top')
        axs[0].set_title(r'Near halls, runtime >18h, average over $\pm10$ runs')
        axs[0].legend(['Sequential', 'Random'], fontsize=12)
        axs[3].set_xlabel('Start time')
        axs[3].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
        fig.tight_layout()
        axs[3].yaxis.set_label_coords(-0.1, 2)
        fig.savefig('acc_DT_eff_avg_comparison_near_bydate.pdf')

        fig, axs = plt.subplots(4, 1, figsize=wide_fig, sharex='col')
        for ax, (hall, det) in zip(axs, near_ads):
            ad_data = get_ad_data(hall, det, data)
            extra_ad_data = get_ad_data(hall, det, extra_data)
            ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
            high_stats = get_ad_data(hall, det, data_random_many)
            high_stats_long = high_stats[high_stats[:, 4] > ns_per_18h]
            high_stats_short = high_stats[high_stats[:, 4] <= ns_per_18h]
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            average = uniform_filter1d(long_runs[:, 5], size=50, mode='nearest')
            avg_high_stats = uniform_filter1d(high_stats_long[:, 5], size=50,
                    mode='nearest')
            ax.plot(get_dates(high_stats_long[:, 3]), high_stats_long[:, 5],
                    '.',
                    markersize=4)
            ax.errorbar(get_dates(high_stats_long[25:-25:50, 3]),
                avg_high_stats[25:-25:50], fmt='.',
                yerr=high_stats_long[25:-25:50, 6]/np.sqrt(50))
            ax.errorbar(get_dates(long_runs[25:-25:50, 3]), average[25:-25:50], fmt='.',
                    yerr=long_runs[25:-25:50, 7]/np.sqrt(50))
            ax.grid(axis='y')
        for ax, name in zip(axs, near_names):
            ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                    verticalalignment='top')
        axs[0].set_title(r'Near halls, random with high $N_{pairs}$')
        axs[0].legend(['Random', 'Avg of Random', 'Avg. of Sequential'], fontsize=12)
        axs[3].set_xlabel('Start time')
        axs[3].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
        fig.tight_layout()
        axs[3].yaxis.set_label_coords(-0.1, 2)
        fig.savefig('acc_DT_eff_avg_comparison_highstats_near_bydate.pdf')

        fig = plt.figure()
        for i, (hall, det) in enumerate(near_ads):
            ad_data = get_ad_data(hall, det, data)
            extra_ad_data = get_ad_data(hall, det, extra_data)
            ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            mpl.rcParams['font.size'] = 12
            ax = fig.add_subplot(220 + i + 1)
            dev = (long_runs[:, 6] - long_runs[:, 5])/long_runs[:, 5]
            ax.plot_date(get_dates(long_runs[:, 3]), dev, markersize=2)
            mpl.rcParams['font.size'] = 18
            stdev = np.std(dev)
            ax.text(0.05, 0.95,
                    f'{near_names[i]}\n'
                    fr'$\mu=$({100*np.mean(dev):.2f}$\pm${100*stdev/np.sqrt(len(dev)):.2f})%' '\n'
                    fr'$\sigma=${stdev*100:.2f}%',
                    fontsize=12,
                    transform=ax.transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='w', alpha=0.5))
        fig.tight_layout()
        fig.savefig('acc_DT_eff_pairing_deviation_near.pdf')

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

        fig, axs = plt.subplots(4, 1, figsize=wide_fig, sharex='col')
        for ax, (hall, det) in zip(axs, far_ads):
            ad_data = get_ad_data(hall, det, data)
            extra_ad_data = get_ad_data(hall, det, extra_data)
            ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            average = uniform_filter1d(long_runs[:, 5], size=20, mode='nearest')
            average_method2 = uniform_filter1d(long_runs[:, 6], size=20, mode='nearest')
            ax.plot_date(get_dates(long_runs[10:-10:20, 3]), average[10:-10:20], '.')
            ax.plot_date(get_dates(long_runs[10:-10:20, 3]), average_method2[10:-10:20], '.')
            ax.grid(axis='y')
        for ax, name in zip(axs, far_names):
            ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                    verticalalignment='top')
        axs[0].set_title(r'Far halls, runtime >18h, average over $\pm10$ runs')
        axs[0].legend(['Sequential', 'Random'], fontsize=12)
        axs[3].set_xlabel('Start time')
        axs[3].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
        fig.tight_layout()
        axs[3].yaxis.set_label_coords(-0.1, 2)
        fig.savefig('acc_DT_eff_avg_comparison_far_bydate.pdf')

        fig, axs = plt.subplots(4, 1, figsize=wide_fig, sharex='col')
        for ax, (hall, det) in zip(axs, far_ads):
            ad_data = get_ad_data(hall, det, data)
            extra_ad_data = get_ad_data(hall, det, extra_data)
            ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
            high_stats = get_ad_data(hall, det, data_random_many)
            high_stats_long = high_stats[high_stats[:, 4] > ns_per_18h]
            high_stats_short = high_stats[high_stats[:, 4] <= ns_per_18h]
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            average = uniform_filter1d(long_runs[:, 5], size=50, mode='nearest')
            avg_high_stats = uniform_filter1d(high_stats_long[:, 5], size=50,
                    mode='nearest')
            ax.plot(get_dates(high_stats_long[:, 3]), high_stats_long[:, 5],
                    '.',
                    markersize=4)
            ax.errorbar(get_dates(high_stats_long[25:-25:50, 3]),
                avg_high_stats[25:-25:50], fmt='.',
                yerr=high_stats_long[25:-25:50, 6]/np.sqrt(50))
            ax.errorbar(get_dates(long_runs[25:-25:50, 3]), average[25:-25:50], fmt='.',
                    yerr=long_runs[25:-25:50, 7]/np.sqrt(50))
            ax.grid(axis='y')
        for ax, name in zip(axs, far_names):
            ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                    verticalalignment='top')
        axs[0].set_title(r'Far halls, random with high $N_{pairs}$')
        axs[0].legend(['Random', 'Avg of Random', 'Avg. of Sequential'], fontsize=12)
        axs[3].set_xlabel('Start time')
        axs[3].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
        fig.tight_layout()
        axs[3].yaxis.set_label_coords(-0.1, 2)
        fig.savefig('acc_DT_eff_avg_comparison_highstats_far_bydate.pdf')

        fig = plt.figure()
        for i, (hall, det) in enumerate(far_ads):
            ad_data = get_ad_data(hall, det, data)
            extra_ad_data = get_ad_data(hall, det, extra_data)
            ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            mpl.rcParams['font.size'] = 12
            ax = fig.add_subplot(220 + i + 1)
            dev = (long_runs[:, 6] - long_runs[:, 5])/long_runs[:, 5]
            ax.plot_date(get_dates(long_runs[:, 3]), dev, markersize=2)
            mpl.rcParams['font.size'] = 18
            stdev = np.std(dev)
            ax.text(0.05, 0.95,
                    f'{far_names[i]}\n'
                    fr'$\mu=$({100*np.mean(dev):.2f}$\pm${100*stdev/np.sqrt(len(dev)):.2f})%' '\n'
                    fr'$\sigma=${stdev*100:.2f}%',
                    fontsize=12,
                    transform=ax.transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='w', alpha=0.5))
        fig.tight_layout()
        fig.savefig('acc_DT_eff_pairing_deviation_far.pdf')

        # Singles and acc-DT-eff
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT RunNo, Hall, DetNo, Start_time,
                    singles_rates.Rate_Hz, Rate_Hz_error, Livetime_ns
                FROM (singles_rates INNER JOIN runs USING (RunNo)) INNER JOIN muon_rates
                    USING (RunNo, DetNo)
                ORDER BY RunNo''')
            singles_data = np.array(cursor.fetchall())
        with sqlite3.Connection('/home/skohn/parameters_resid_flashers.db') as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT RunNo, Hall, DetNo, Start_time,
                    singles_rates.Rate_Hz, Rate_Hz_error, Livetime_ns
                FROM (singles_rates INNER JOIN runs USING (RunNo)) INNER JOIN muon_rates
                    USING (RunNo, DetNo)
                ORDER BY RunNo''')
            singles_data_resid = np.array(cursor.fetchall())

        fig, axs = plt.subplots(4, 1, figsize=big_fig, sharex='col')
        axs_singles = []
        for ax, (hall, det) in zip(axs, far_ads):
            ad_data = get_ad_data(hall, det, data)
            extra_ad_data = get_ad_data(hall, det, extra_data)
            ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
            high_stats = get_ad_data(hall, det, data_random_many)
            high_stats_long = high_stats[high_stats[:, 4] > ns_per_18h]
            high_stats_short = high_stats[high_stats[:, 4] <= ns_per_18h]
            long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
            average = uniform_filter1d(long_runs[:, 5], size=50, mode='nearest')
            avg_high_stats = uniform_filter1d(high_stats_long[:, 5], size=50,
                    mode='nearest')
            ax.plot(get_dates(high_stats_long[:, 3]), high_stats_long[:, 5],
                    '.',
                    markersize=4)
            ax.errorbar(get_dates(high_stats_long[25:-25:50, 3]),
                avg_high_stats[25:-25:50], fmt='.',
                yerr=high_stats_long[25:-25:50, 6]/np.sqrt(50))
            ax.errorbar(get_dates(long_runs[25:-25:50, 3]), average[25:-25:50], fmt='.',
                    yerr=long_runs[25:-25:50, 7]/np.sqrt(50))
            ax.grid(axis='y')

            ad_singles_data = get_ad_data(hall, det, singles_data)
            ax_singles = ax.twinx()
            axs_singles.append(ax_singles)
            ad_singles_long = ad_singles_data[ad_singles_data[:, 6] > ns_per_18h]
            ad_singles_long = ad_singles_data[ad_singles_data[:, 3] >
                    1343377097272373062]
            ax_singles.plot_date(get_dates(ad_singles_long[:, 3]), ad_singles_long[:, 4],
                    color=colors[3], markersize=2)
            resid_ad_singles_data = get_ad_data(hall, det, singles_data_resid)
            resid_ad_singles_long = resid_ad_singles_data[resid_ad_singles_data[:, 6] > ns_per_18h]
            resid_ad_singles_long = resid_ad_singles_data[resid_ad_singles_data[:, 3] >
                    1343377097272373062]
            ax_singles.plot_date(get_dates(resid_ad_singles_long[:, 3]), resid_ad_singles_long[:, 4],
                    color=colors[4], markersize=2)

        for ax, name in zip(axs, far_names):
            ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                    verticalalignment='top')
        axs[0].set_title(r'Far halls, random with high $N_{pairs}$')
        axs[3].legend(['Random', 'Avg of Random', 'Avg. of Sequential'], fontsize=12)
        axs[3].set_xlabel('Start time')
        axs[3].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
        axs_singles[3].set_ylabel('Singles Rate [Hz]')
        fig.tight_layout()
        axs[3].yaxis.set_label_coords(-0.1, 2)
        axs_singles[3].yaxis.set_label_coords(1.05, 2)
        axs_singles[3].legend(['Singles', 'Singles w/ resid. flasher cut'], fontsize=12)
        fig.savefig('acc_DT_eff_far_compare_singles_8ad.pdf')
