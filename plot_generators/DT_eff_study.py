import argparse
from datetime import datetime
import sqlite3

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.ndimage.filters import uniform_filter1d
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    args = parser.parse_args()
    database = args.database
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
        cursor.execute('''SELECT RunNo, Hall, DetNo, Start_time, Livetime_ns, DistanceTime_DT_Eff,
            DistanceTime_DT_Eff_error
            FROM runs NATURAL JOIN distance_time_eff_study NATURAL JOIN muon_rates
            WHERE PairingType = "random_many_resid_flasher"
            ORDER BY RunNo, DetNo''')
        data_random_resid_flasher = np.array(cursor.fetchall())
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
    ns_per_18h = 60e9 * 60 * 18
    wide_fig = (16, 8)
    big_fig = (16, 10)

    # fig, axs = plt.subplots(4, 1, figsize=wide_fig, sharex='col')
    # for ax, (hall, det) in zip(axs, near_ads):
        # ad_data = get_ad_data(hall, det, data)
        # extra_ad_data = get_ad_data(hall, det, extra_data)
        # ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
        # long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        # average = uniform_filter1d(long_runs[:, 5], size=20, mode='nearest')
        # average_method2 = uniform_filter1d(long_runs[:, 6], size=20, mode='nearest')
        # ax.plot_date(get_dates(long_runs[10:-10:20, 3]), average[10:-10:20], '.')
        # ax.errorbar(get_dates(long_runs[10:-10:20, 3]), average_method2[10:-10:20],
                # fmt='.', yerr=long_runs[10:-10:20, 7]/np.sqrt(20))
        # ax.grid(axis='y')
    # for ax, name in zip(axs, near_names):
        # ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                # verticalalignment='top')
    # axs[0].set_title(r'Near halls, runtime >18h, average over $\pm10$ runs')
    # axs[0].legend(['Sequential', 'Random'], fontsize=12)
    # axs[3].set_xlabel('Start time')
    # axs[3].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
    # fig.tight_layout()
    # axs[3].yaxis.set_label_coords(-0.1, 2)
    # fig.savefig('acc_DT_eff_avg_comparison_near_bydate.pdf')

    # fig, axs = plt.subplots(4, 1, figsize=wide_fig, sharex='col')
    # for ax, (hall, det) in zip(axs, near_ads):
        # ad_data = get_ad_data(hall, det, data)
        # extra_ad_data = get_ad_data(hall, det, extra_data)
        # ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
        # high_stats = get_ad_data(hall, det, data_random_many)
        # high_stats_long = high_stats[high_stats[:, 4] > ns_per_18h]
        # high_stats_short = high_stats[high_stats[:, 4] <= ns_per_18h]
        # long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        # average = uniform_filter1d(long_runs[:, 5], size=50, mode='nearest')
        # avg_high_stats = uniform_filter1d(high_stats_long[:, 5], size=50,
                # mode='nearest')
        # ax.plot(get_dates(high_stats_long[:, 3]), high_stats_long[:, 5],
                # '.',
                # markersize=4)
        # ax.errorbar(get_dates(high_stats_long[25:-25:50, 3]),
            # avg_high_stats[25:-25:50], fmt='.',
            # yerr=high_stats_long[25:-25:50, 6]/np.sqrt(50))
        # ax.errorbar(get_dates(long_runs[25:-25:50, 3]), average[25:-25:50], fmt='.',
                # yerr=long_runs[25:-25:50, 7]/np.sqrt(50))
        # ax.grid(axis='y')
    # for ax, name in zip(axs, near_names):
        # ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                # verticalalignment='top')
    # axs[0].set_title(r'Near halls, random with high $N_{pairs}$')
    # axs[0].legend(['Random', 'Avg of Random', 'Avg. of Sequential'], fontsize=12)
    # axs[3].set_xlabel('Start time')
    # axs[3].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
    # fig.tight_layout()
    # axs[3].yaxis.set_label_coords(-0.1, 2)
    # fig.savefig('acc_DT_eff_avg_comparison_highstats_near_bydate.pdf')

    # fig = plt.figure()
    # for i, (hall, det) in enumerate(near_ads):
        # ad_data = get_ad_data(hall, det, data)
        # extra_ad_data = get_ad_data(hall, det, extra_data)
        # ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
        # long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        # mpl.rcParams['font.size'] = 12
        # ax = fig.add_subplot(220 + i + 1)
        # dev = (long_runs[:, 6] - long_runs[:, 5])/long_runs[:, 5]
        # ax.plot_date(get_dates(long_runs[:, 3]), dev, markersize=2)
        # mpl.rcParams['font.size'] = 18
        # stdev = np.std(dev)
        # ax.text(0.05, 0.95,
                # f'{near_names[i]}\n'
                # fr'$\mu=$({100*np.mean(dev):.2f}$\pm${100*stdev/np.sqrt(len(dev)):.2f})%' '\n'
                # fr'$\sigma=${stdev*100:.2f}%',
                # fontsize=12,
                # transform=ax.transAxes, verticalalignment='top',
                # bbox=dict(boxstyle='round', facecolor='w', alpha=0.5))
    # fig.tight_layout()
    # fig.savefig('acc_DT_eff_pairing_deviation_near.pdf')

    # fig, axs = plt.subplots(4, 1, figsize=wide_fig, sharex='col')
    # for ax, (hall, det) in zip(axs, far_ads):
        # ad_data = get_ad_data(hall, det, data)
        # extra_ad_data = get_ad_data(hall, det, extra_data)
        # ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
        # long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        # average = uniform_filter1d(long_runs[:, 5], size=20, mode='nearest')
        # average_method2 = uniform_filter1d(long_runs[:, 6], size=20, mode='nearest')
        # ax.plot_date(get_dates(long_runs[10:-10:20, 3]), average[10:-10:20], '.')
        # ax.plot_date(get_dates(long_runs[10:-10:20, 3]), average_method2[10:-10:20], '.')
        # ax.grid(axis='y')
    # for ax, name in zip(axs, far_names):
        # ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                # verticalalignment='top')
    # axs[0].set_title(r'Far halls, runtime >18h, average over $\pm10$ runs')
    # axs[0].legend(['Sequential', 'Random'], fontsize=12)
    # axs[3].set_xlabel('Start time')
    # axs[3].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
    # fig.tight_layout()
    # axs[3].yaxis.set_label_coords(-0.1, 2)
    # fig.savefig('acc_DT_eff_avg_comparison_far_bydate.pdf')

    # fig, axs = plt.subplots(4, 1, figsize=wide_fig, sharex='col')
    # for ax, (hall, det) in zip(axs, far_ads):
        # ad_data = get_ad_data(hall, det, data)
        # extra_ad_data = get_ad_data(hall, det, extra_data)
        # ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
        # high_stats = get_ad_data(hall, det, data_random_many)
        # high_stats_long = high_stats[high_stats[:, 4] > ns_per_18h]
        # high_stats_short = high_stats[high_stats[:, 4] <= ns_per_18h]
        # long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        # average = uniform_filter1d(long_runs[:, 5], size=50, mode='nearest')
        # avg_high_stats = uniform_filter1d(high_stats_long[:, 5], size=50,
                # mode='nearest')
        # ax.plot(get_dates(high_stats_long[:, 3]), high_stats_long[:, 5],
                # '.',
                # markersize=4)
        # ax.errorbar(get_dates(high_stats_long[25:-25:50, 3]),
            # avg_high_stats[25:-25:50], fmt='.',
            # yerr=high_stats_long[25:-25:50, 6]/np.sqrt(50))
        # ax.errorbar(get_dates(long_runs[25:-25:50, 3]), average[25:-25:50], fmt='.',
                # yerr=long_runs[25:-25:50, 7]/np.sqrt(50))
        # ax.grid(axis='y')
    # for ax, name in zip(axs, far_names):
        # ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                # verticalalignment='top')
    # axs[0].set_title(r'Far halls, random with high $N_{pairs}$')
    # axs[0].legend(['Random', 'Avg of Random', 'Avg. of Sequential'], fontsize=12)
    # axs[3].set_xlabel('Start time')
    # axs[3].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
    # fig.tight_layout()
    # axs[3].yaxis.set_label_coords(-0.1, 2)
    # fig.savefig('acc_DT_eff_avg_comparison_highstats_far_bydate.pdf')

    # fig = plt.figure()
    # for i, (hall, det) in enumerate(far_ads):
        # ad_data = get_ad_data(hall, det, data)
        # extra_ad_data = get_ad_data(hall, det, extra_data)
        # ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
        # long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        # mpl.rcParams['font.size'] = 12
        # ax = fig.add_subplot(220 + i + 1)
        # dev = (long_runs[:, 6] - long_runs[:, 5])/long_runs[:, 5]
        # ax.plot_date(get_dates(long_runs[:, 3]), dev, markersize=2)
        # mpl.rcParams['font.size'] = 18
        # stdev = np.std(dev)
        # ax.text(0.05, 0.95,
                # f'{far_names[i]}\n'
                # fr'$\mu=$({100*np.mean(dev):.2f}$\pm${100*stdev/np.sqrt(len(dev)):.2f})%' '\n'
                # fr'$\sigma=${stdev*100:.2f}%',
                # fontsize=12,
                # transform=ax.transAxes, verticalalignment='top',
                # bbox=dict(boxstyle='round', facecolor='w', alpha=0.5))
    # fig.tight_layout()
    # fig.savefig('acc_DT_eff_pairing_deviation_far.pdf')

    fig, axs = plt.subplots(4, 2, figsize=big_fig, sharex='col')
    axs_singles = []
    for i, (name, (hall, det)) in enumerate(zip(near_names, near_ads)):
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
        ax = axs[(i // 2) * 2, i % 2]
        ax.plot(get_dates(high_stats_long[:, 3]), high_stats_long[:, 5],
                '.',
                markersize=4,
                alpha=0.5)
        ax.grid(axis='y')
        ax.set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=18)

        ad_singles_data = get_ad_data(hall, det, singles_data)
        ax_singles = axs[(i // 2) * 2 + 1, i % 2]
        axs_singles.append(ax_singles)
        ad_singles_long = ad_singles_data[ad_singles_data[:, 6] > ns_per_18h]
        ax_singles.plot_date(get_dates(ad_singles_long[:, 3]), ad_singles_long[:, 4],
                markersize=2)
        resid_ad_singles_data = get_ad_data(hall, det, singles_data_resid)
        resid_ad_singles_long = resid_ad_singles_data[resid_ad_singles_data[:, 6] > ns_per_18h]
        ax_singles.plot_date(get_dates(resid_ad_singles_long[:, 3]), resid_ad_singles_long[:, 4],
                markersize=2)
        ax_singles.set_ylabel('Singles Rate [Hz]', fontsize=15)

        ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                verticalalignment='top')
        ax_singles.text(0.02, 0.95, name, fontsize=12, transform=ax_singles.transAxes,
                verticalalignment='top')

    axs = axs.flatten()
    axs[0].set_title(r'Near halls, random with high $N_{pairs}$')
    axs[0].legend([r'$\varepsilon_{DT,acc}$ (Random pairing)'], fontsize=12)
    fig.tight_layout()
    axs_singles[0].legend(['Nominal', 'Resid. flasher cut'], fontsize=12)
    fig.savefig('acc_DT_eff_near_compare_singles.pdf')

    fig, axs = plt.subplots(4, 2, figsize=big_fig, sharex='col')
    axs_singles = []
    for i, (name, (hall, det)) in enumerate(zip(far_names, far_ads)):
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
        ax = axs[(i // 2) * 2, i % 2]
        ax.plot(get_dates(high_stats_long[:, 3]), high_stats_long[:, 5],
                '.',
                markersize=4,
                alpha=0.5)
        ax.grid(axis='y')
        ax.set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=18)

        ad_singles_data = get_ad_data(hall, det, singles_data)
        ax_singles = axs[(i // 2) * 2 + 1, i % 2]
        axs_singles.append(ax_singles)
        ad_singles_long = ad_singles_data[ad_singles_data[:, 6] > ns_per_18h]
        ad_singles_long = ad_singles_data[ad_singles_data[:, 3] >
                1343377097272373062]
        ax_singles.plot_date(get_dates(ad_singles_long[:, 3]), ad_singles_long[:, 4],
                markersize=2)
        resid_ad_singles_data = get_ad_data(hall, det, singles_data_resid)
        resid_ad_singles_long = resid_ad_singles_data[resid_ad_singles_data[:, 6] > ns_per_18h]
        resid_ad_singles_long = resid_ad_singles_data[resid_ad_singles_data[:, 3] >
                1343377097272373062]
        ax_singles.plot_date(get_dates(resid_ad_singles_long[:, 3]), resid_ad_singles_long[:, 4],
                markersize=2)
        ax_singles.set_ylabel('Singles Rate [Hz]', fontsize=15)
        ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                verticalalignment='top')
        ax_singles.text(0.02, 0.95, name, fontsize=12, transform=ax_singles.transAxes,
                verticalalignment='top')

    axs = axs.flatten()
    axs[0].set_title(r'Far halls, random with high $N_{pairs}$')
    axs[0].legend([r'$\varepsilon_{DT,acc}$ (Random pairing)'], fontsize=12)
    fig.tight_layout()
    axs_singles[0].legend(['Nominal', 'Resid. flasher cut'], fontsize=12)
    fig.savefig('acc_DT_eff_far_compare_singles.pdf')

    fig, axs = plt.subplots(2, 2, figsize=big_fig, sharex='col')
    axs_singles = []
    for ax, (hall, det) in zip(axs.flat, near_ads):
        ad_data = get_ad_data(hall, det, data)
        extra_ad_data = get_ad_data(hall, det, extra_data)
        ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
        high_stats = get_ad_data(hall, det, data_random_many)
        high_stats_long = high_stats[high_stats[:, 4] > ns_per_18h]
        high_stats_short = high_stats[high_stats[:, 4] <= ns_per_18h]
        resid_flashers = get_ad_data(hall, det, data_random_resid_flasher)
        resid_flashers_long = resid_flashers[resid_flashers[:, 4] > ns_per_18h]
        long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        average = uniform_filter1d(long_runs[:, 5], size=50, mode='nearest')
        avg_high_stats = uniform_filter1d(high_stats_long[:, 5], size=50,
                mode='nearest')
        ax.plot(get_dates(high_stats_long[:, 3]), high_stats_long[:, 5],
                '.',
                markersize=4)
        ax.plot(get_dates(resid_flashers_long[:, 3]), resid_flashers_long[:, 5],
                '.',
                markersize=4)
        ax.grid(axis='y')

    for ax, name in zip(axs.flat, near_names):
        ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                verticalalignment='top')
    axs = axs.flatten()
    axs[0].set_title(r'Near halls, random with high $N_{pairs}$')
    axs[0].legend(['Nominal (random)', 'Applied cut on residual flashers'], fontsize=15)
    axs[3].set_xlabel('Start time')
    axs[2].set_xlabel('Start time')
    axs[2].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
    fig.tight_layout()
    axs[2].yaxis.set_label_coords(-0.15, 1.18)
    fig.savefig('acc_DT_eff_near_compare_resid_flashers.pdf')

    fig, axs = plt.subplots(2, 2, figsize=big_fig, sharex='col')
    axs_singles = []
    for ax, (hall, det) in zip(axs.flat, far_ads):
        ad_data = get_ad_data(hall, det, data)
        extra_ad_data = get_ad_data(hall, det, extra_data)
        ad_data = np.hstack((ad_data, extra_ad_data[:, 3:]))
        high_stats = get_ad_data(hall, det, data_random_many)
        high_stats_long = high_stats[high_stats[:, 4] > ns_per_18h]
        high_stats_short = high_stats[high_stats[:, 4] <= ns_per_18h]
        resid_flashers = get_ad_data(hall, det, data_random_resid_flasher)
        resid_flashers_long = resid_flashers[resid_flashers[:, 4] > ns_per_18h]
        long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        average = uniform_filter1d(long_runs[:, 5], size=50, mode='nearest')
        avg_high_stats = uniform_filter1d(high_stats_long[:, 5], size=50,
                mode='nearest')
        ax.plot(get_dates(high_stats_long[:, 3]), high_stats_long[:, 5],
                '.',
                markersize=4)
        ax.plot(get_dates(resid_flashers_long[:, 3]), resid_flashers_long[:, 5],
                '.',
                markersize=4)
        ax.grid(axis='y')

    for ax, name in zip(axs.flat, far_names):
        ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                verticalalignment='top')
    axs = axs.flatten()
    axs[0].set_title(r'Far halls, random with high $N_{pairs}$')
    axs[0].legend(['Nominal (random)', 'Applied cut on residual flashers'], fontsize=15)
    axs[3].set_xlabel('Start time')
    axs[2].set_xlabel('Start time')
    axs[2].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
    fig.tight_layout()
    axs[2].yaxis.set_label_coords(-0.23, 1.18)
    fig.savefig('acc_DT_eff_far_compare_resid_flashers.pdf')


    # Better joins between pairing types
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        # Table with runs for which there is both a "random_many"
        # and a "random_many_resid_flasher" computed acc DT efficiency.
        cursor.execute('''SELECT nom.RunNo, Hall, nom.DetNo, Start_time, Livetime_ns,
            nom.DistanceTime_DT_Eff,
            noflash.DistanceTime_DT_Eff,
            fl.Rate_Hz
            FROM (runs NATURAL JOIN muon_rates) AS r
                INNER JOIN
                residual_flasher_rates fl ON r.DetNo = fl.DetNo
                    AND r.RunNo = fl.RunNo
                INNER JOIN (
                    (SELECT RunNo, DetNo, DistanceTime_DT_Eff
                        FROM distance_time_eff_study
                        WHERE PairingType = "random_many") AS nom
                    INNER JOIN
                    (SELECT RunNo, DetNo, DistanceTime_DT_Eff
                       FROM distance_time_eff_study
                       WHERE PairingType = "random_many_resid_flasher") AS noflash
                    ON nom.RunNo = noflash.RunNo
                        AND nom.DetNo = noflash.DetNo
                ) AS dt ON r.RunNo = dt.RunNo AND r.DetNo = dt.DetNo
            ORDER BY nom.RunNo, nom.DetNo''')
        data_paired = np.array(cursor.fetchall())
        cursor.execute('''SELECT runs.RunNo, Hall, fl.DetNo,
                Start_time, Livetime_ns, fl.Rate_Hz
            FROM runs
            NATURAL JOIN muon_rates
            INNER JOIN residual_flasher_rates fl ON muon_rates.DetNo = fl.DetNo
                AND runs.RunNo = fl.RunNO
            ORDER BY runs.RunNo, fl.DetNo''')
        flasher_rates = np.array(cursor.fetchall())

    fig, axs = plt.subplots(4, 2, figsize=big_fig, sharex='col')
    axs_flashers = []
    for i, (name, (hall, det)) in enumerate(zip(near_names, near_ads)):
        ad_data = get_ad_data(hall, det, data_paired)
        long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        ad_flasher_rates = get_ad_data(hall, det, flasher_rates)
        ad_flasher_rates_long = ad_flasher_rates[ad_flasher_rates[:, 4] > ns_per_18h]
        ax = axs[(i // 2) * 2, i % 2]
        ax.plot(get_dates(long_runs[:, 3]), long_runs[:, 6] - long_runs[:, 5],
                '.',
                color=colors[0],
                markersize=4)
        ax.set_ylabel(r'$\Delta \varepsilon_{DT,acc}$', fontsize=20)
        ax_flashers = axs[(i // 2) * 2 + 1, i % 2]
        ax_flashers.plot(get_dates(ad_flasher_rates_long[:, 3]),
                ad_flasher_rates_long[:, 5],
                '.',
                color=colors[1],
                markersize=4)
        axs_flashers.append(ax_flashers)
        ax.grid(axis='both')
        ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                verticalalignment='top')
        ax_flashers.text(0.02, 0.95, name, fontsize=12,
                transform=ax_flashers.transAxes,
                verticalalignment='top')
        ax_flashers.set_ylabel(r'$R_{resid.\ flasher}$ [Hz]', fontsize=15)
        ax_flashers.grid(axis='x')

    axs = axs.flatten()
    axs[0].set_title('Near halls')
    # axs[0].legend(['Nominal (random)', 'Applied cut on residual flashers'], fontsize=15)
    axs[6].set_xlabel('Start time')
    axs[7].set_xlabel('Start time')
    # axs[2].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
    fig.tight_layout()
    # axs[2].yaxis.set_label_coords(-0.23, 1.18)
    # axs_flashers[3].yaxis.set_label_coords(1.05, 1.18)
    fig.savefig('acc_DT_eff_near_resid_flasher_rates.pdf')

    fig, axs = plt.subplots(2, 2, figsize=big_fig, sharex='col')
    axs_flashers = []
    for ax, name, (hall, det) in zip(axs.flat, near_names, near_ads):
        ad_data = get_ad_data(hall, det, data_paired)
        long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        ad_flasher_rates = get_ad_data(hall, det, flasher_rates)
        ad_flasher_rates_long = ad_flasher_rates[ad_flasher_rates[:, 4] > ns_per_18h]
        ax.plot(long_runs[:, 5] - long_runs[:, 6],
                long_runs[:, 7],
                '.',
                color=colors[0],
                markersize=4)
        ax.set_ylabel(r'$R_{resid.\ flasher}$ [Hz]', fontsize=20)
        ax.grid()
        ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                verticalalignment='top')

    axs = axs.flatten()
    axs[0].set_title('Near halls')
    axs[2].set_xlabel(r'$\Delta \varepsilon_{DT,acc}$', fontsize=20)
    axs[3].set_xlabel(r'$\Delta \varepsilon_{DT,acc}$', fontsize=20)
    axs[2].tick_params(axis='x', labelrotation=15)
    axs[3].tick_params(axis='x', labelrotation=15)
    fig.tight_layout()
    fig.savefig('acc_DT_eff_near_resid_flasher_corr.pdf')

    fig, axs = plt.subplots(4, 2, figsize=big_fig, sharex='col')
    axs_flashers = []
    for i, (name, (hall, det)) in enumerate(zip(far_names, far_ads)):
        ad_data = get_ad_data(hall, det, data_paired)
        long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        ad_flasher_rates = get_ad_data(hall, det, flasher_rates)
        ad_flasher_rates_long = ad_flasher_rates[ad_flasher_rates[:, 4] > ns_per_18h]
        ax = axs[(i // 2) * 2, i % 2]
        ax.plot(get_dates(long_runs[:, 3]), long_runs[:, 6] - long_runs[:, 5],
                '.',
                color=colors[0],
                markersize=4)
        ax.set_ylabel(r'$\Delta \varepsilon_{DT,acc}$', fontsize=20)
        ax_flashers = axs[(i // 2) * 2 + 1, i % 2]
        ax_flashers.plot(get_dates(ad_flasher_rates_long[:, 3]),
                ad_flasher_rates_long[:, 5],
                '.',
                color=colors[1],
                markersize=4)
        axs_flashers.append(ax_flashers)
        ax.grid(axis='both')
        ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                verticalalignment='top')
        ax_flashers.text(0.02, 0.95, name, fontsize=12,
                transform=ax_flashers.transAxes,
                verticalalignment='top')
        ax_flashers.set_ylabel(r'$R_{resid.\ flasher}$ [Hz]', fontsize=15)
        ax_flashers.grid(axis='x')

    axs = axs.flatten()
    axs[0].set_title('Far halls')
    # axs[0].legend(['Nominal (random)', 'Applied cut on residual flashers'], fontsize=15)
    axs[6].set_xlabel('Start time')
    axs[7].set_xlabel('Start time')
    # axs[2].set_ylabel(r'$\varepsilon_{DT,acc}$', fontsize=20)
    fig.tight_layout()
    # axs[2].yaxis.set_label_coords(-0.23, 1.18)
    # axs_flashers[3].yaxis.set_label_coords(1.05, 1.18)
    fig.savefig('acc_DT_eff_far_resid_flasher_rates.pdf')

    fig, axs = plt.subplots(2, 2, figsize=big_fig, sharex='col')
    axs_flashers = []
    for ax, name, (hall, det) in zip(axs.flat, far_names, far_ads):
        ad_data = get_ad_data(hall, det, data_paired)
        long_runs = ad_data[ad_data[:, 4] > ns_per_18h]
        ad_flasher_rates = get_ad_data(hall, det, flasher_rates)
        ad_flasher_rates_long = ad_flasher_rates[ad_flasher_rates[:, 4] > ns_per_18h]
        ax.plot(long_runs[:, 5] - long_runs[:, 6],
                long_runs[:, 7],
                '.',
                color=colors[0],
                markersize=4)
        ax.set_ylabel(r'$R_{resid.\ flasher}$ [Hz]', fontsize=20)
        ax.grid()
        ax.text(0.02, 0.95, name, fontsize=12, transform=ax.transAxes,
                verticalalignment='top')

    axs = axs.flatten()
    axs[0].set_title('Far halls')
    axs[2].set_xlabel(r'$\Delta \varepsilon_{DT,acc}$', fontsize=20)
    axs[3].set_xlabel(r'$\Delta \varepsilon_{DT,acc}$', fontsize=20)
    axs[2].tick_params(axis='x', labelrotation=15)
    axs[3].tick_params(axis='x', labelrotation=15)
    fig.tight_layout()
    fig.savefig('acc_DT_eff_far_resid_flasher_corr.pdf')

    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT RunNo, Hall, DetNo, Start_time, Livetime_ns,
            NumFirstHour, NumLastHour
            FROM runs NATURAL JOIN muon_rates NATURAL JOIN singles_within_run
            ORDER BY RunNo, DetNo''')
        data = np.array(cursor.fetchall())
    fig, axs = plt.subplots(2, 2)
