"""Plot the light yield values based on txt files."""

import argparse
from collections import defaultdict
import datetime
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

UTC_TO_CST = datetime.timedelta(hours=8)
def timestamp_to_datetime64(timestamp):
    utc_datetime = datetime.datetime.utcfromtimestamp(timestamp)
    cst_datetime = utc_datetime + UTC_TO_CST
    return np.datetime64(cst_datetime.isoformat())


def main(file_template, output):
    light_yields = {}
    for hall in [1, 2, 3]:
        data_lists = defaultdict(list)
        filename = file_template.format(site=hall)
        with open(filename, 'r') as f:
            marked_summer_12_start = defaultdict(bool)
            marked_summer_12_end = defaultdict(bool)
            marked_winter_17_start = defaultdict(bool)
            marked_winter_17_end = defaultdict(bool)
            for line in f:
                splits = line.split()
                det = int(splits[0])
                date = timestamp_to_datetime64(int(splits[1]))
                if date > P17B_MAX:
                    continue
                if date < SUMMER_12_END and (
                    (hall, det) == (2, 2) or (hall, det) == (3, 4)
                ):
                    continue
                if date > SUMMER_12_START and date < SUMMER_12_END:
                    if not marked_summer_12_start[det]:
                        marked_summer_12_start[det] = True
                        data_lists[det].append((date, np.nan, 1))
                    else:
                        pass
                    continue
                if date > SUMMER_12_END and not marked_summer_12_end[det]:
                    data_lists[det].append((date, np.nan, 1))
                    marked_summer_12_end[det] = True
                if date > WINTER_17_START and date < WINTER_17_END:
                    if not marked_winter_17_start[det]:
                        marked_winter_17_start[det] = True
                        data_lists[det].append((date, np.nan, 1))
                    else:
                        pass
                    continue
                if date > WINTER_17_END and not marked_winter_17_end[det]:
                    data_lists[det].append((date, np.nan, 1))
                    marked_winter_17_end[det] = True
                if date > WINTER_17_START and (hall, det) == (1, 1):
                    continue
                if float(splits[3]) == 0:
                    continue
                data_lists[det].append((
                    date,
                    float(splits[2]),
                    float(splits[3]),
                ))
            data = {
                (hall, det): np.array(
                    det_data,
                    dtype=[
                        ('time', 'datetime64[h]'),
                        ('light_yield', float),
                        ('error', float),
                    ]
                )
                for det, det_data in data_lists.items()
            }
        light_yields.update(data)

    fig, ax = plt.subplots()
    for (hall, det), data in light_yields.items():
        ax.plot(data['time'], data['light_yield'], label=f'EH{hall}-AD{det}')
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1, 1))
    ax.set_xlabel('Date')
    ax.set_ylabel('Light yield [PE / MeV]')
    plt.autoscale()
    plt.show()
    fig.savefig('light_yield.pdf', bbox_inches='tight')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_template', help='using {site} and {ad} for EH{site}-AD{ad}')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    main(args.file_template, args.output)
