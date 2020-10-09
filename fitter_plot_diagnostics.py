import argparse
from datetime import datetime
import sqlite3

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.ndimage.filters import uniform_filter1d
import numpy as np

import prediction as pred

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config')
    choices = ('all',)
    parser.add_argument('--types', nargs='+', choices=choices, default='all',
            metavar='PLOT_TYPE', help=f'Choices: {choices}')
    args = parser.parse_args()
    constants = pred.load_constants(args.config)
    plot_types = args.types
    plot_all = 'all' in plot_types
    mpl.rcParams['font.size'] = 18
    mpl.rcParams['figure.figsize'] = (8, 8)
    mpl.rcParams['lines.linewidth'] = 3
    mpl.rcParams['lines.markersize'] = 10
    mpl.rcParams['errorbar.capsize'] = 5
    wide_fig = (16, 8)
    big_fig = (16, 10)
    colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
    dateconv_s = np.vectorize(datetime.fromtimestamp)
    def dateconv_ns(timestamp_ns_array):
        return dateconv_s(timestamp_ns_array/1e9)
    def get_ad_data(hall, det, in_data):
        return in_data[(in_data[:, 1] == hall) & (in_data[:, 2] == det)]
    ads = pred.all_ads
    near_ads = pred.near_ads
    far_ads = pred.far_ads
    def name(hall, det):
        return f'EH{hall}-AD{det}'
    ad_names = [name(hall, det) for hall, det in ads]
    ad_names_2line = [name(hall, det).replace('-', '\n') for hall, det
            in ads]
    near_names = [name(hall, det) for hall, det in near_ads]
    far_names = [name(hall, det) for hall, det in far_ads]
