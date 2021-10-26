"""Plot the gain values based on txt files."""

import argparse
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import common

P17B_MAX = np.datetime64('2017-08-31T23:59:59')
SUMMER_12_START = np.datetime64('2012-07-28T23:59:59')
SUMMER_12_END = np.datetime64('2012-10-19T00:00:00')
WINTER_17_START = np.datetime64('2016-12-20T23:59:59')
WINTER_17_END = np.datetime64('2017-01-26T00:00:00')


def main(file_template, output):
    gains = {}
    for hall, det in common.all_ads:
        data_list = []
        filename = file_template.format(site=hall, ad=det)
        with open(filename, 'r') as f:
            marked_summer_12_start = False
            marked_summer_12_end = False
            marked_winter_17_start = False
            marked_winter_17_end = False
            for line in f:
                splits = line.split()
                date = np.datetime64(splits[0] + ' ' + splits[1], 'h')
                if date > P17B_MAX:
                    continue
                if date < SUMMER_12_END and (
                    (hall, det) == (2, 2) or (hall, det) == (3, 4)
                ):
                    continue
                if date > SUMMER_12_START and date < SUMMER_12_END:
                    if not marked_summer_12_start:
                        marked_summer_12_start = True
                        data_list.append((date, np.nan, 1))
                    else:
                        pass
                    continue
                if date > SUMMER_12_END and not marked_summer_12_end:
                    data_list.append((date, np.nan, 1))
                    marked_summer_12_end = True
                if date > WINTER_17_START and date < WINTER_17_END:
                    if not marked_winter_17_start:
                        marked_winter_17_start = True
                        data_list.append((date, np.nan, 1))
                    else:
                        pass
                    continue
                if date > WINTER_17_END and not marked_winter_17_end:
                    data_list.append((date, np.nan, 1))
                    marked_winter_17_end = True
                if date > WINTER_17_START and (hall, det) == (1, 1):
                    break
                data_list.append((
                    date,
                    float(splits[2]),
                    float(splits[3]),
                ))
        data = np.array(
            data_list,
            dtype=[
                ('time', 'datetime64[h]'),
                ('gain', float),
                ('rms', float),
            ]
        )
        gains[hall, det] = data

    fig, ax = plt.subplots()
    for (hall, det), data in gains.items():
        ax.plot(data['time'], data['gain'], label=f'EH{hall}-AD{det}')
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xlabel('Date')
    ax.set_ylabel('Gain [ADC counts / PE]')
    plt.autoscale()
    plt.show()
    fig.savefig('gains.pdf', bbox_inches='tight')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_template', help='using {site} and {ad} for EH{site}-AD{ad}')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    main(args.file_template, args.output)
