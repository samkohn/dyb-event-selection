"""Dump LBNL ToyMC prompt E histos to the parameters DB."""

import argparse
import json
import sqlite3

from prediction import all_ads

def main(infilename, entry_number, update_db, source, binning_type):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    host_ttree = infile.Get('tr')
    host_ttree.GetEntry(entry_number)
    all_hists = {}
    for i, halldet in enumerate(all_ads):
        for stage in (2,):
            name = f'h_stage{stage}_ad{i+1}'
            all_hists[halldet] = getattr(host_ttree, name)
    num_ibds_default = {}
    bins_default = {}
    for halldet, hist in all_hists.items():
        nbins = hist.GetNbinsX()
        num_ibds_halldet = []
        bins_halldet = []
        for bin_index in range(1, nbins+1):
            value = hist.GetBinContent(bin_index)
            num_ibds_halldet.append(value)
            bin_low_edge = int(hist.GetBinLowEdge(bin_index) * 1000)
            bins_halldet.append(bin_low_edge)
        bins_halldet.append(int(hist.GetBinLowEdge(nbins + 1) * 1000))  # last bin upper edge
        num_ibds_default[halldet] = num_ibds_halldet
        bins_default[halldet] = bins_halldet
    if binning_type == 'default':
        num_ibds = num_ibds_default
        bins = bins_default
    elif binning_type == 'nH nominal':
        num_ibds = {}
        # 1.5, 1.6, 1.8, 2.0, 2.2, ..., 7.8, 8.0, 12 MeV
        bins_halldet = [1500] + list(range(1600, 8001, 200)) + [12000]  # keV
        bins = {}
        for halldet, num_ibds_binned in num_ibds_default.items():
            num_ibds_halldet = []
            split_bin = num_ibds_binned[3]  # 1.4, 1.6 MeV
            next_bin = num_ibds_binned[4]  # 1.6, 1.8 MeV
            # Linear interpolation: need # from 1.5 to 1.6 MeV
            val_at_1_5 = split_bin
            val_at_1_6 = (split_bin + next_bin)/2
            # Trapezoid volume: 1/2 (b1 + b2) * width
            trap_volume = 0.5 * (val_at_1_5 + val_at_1_6) * 0.5
            num_ibds_halldet.append(trap_volume)
            for val in num_ibds_binned[4:]:
                num_ibds_halldet.append(val)
            num_ibds[halldet] = num_ibds_halldet
            bins[halldet] = bins_halldet
    elif binning_type == 'nH modified 1':
        num_ibds = {}
        bins_halldet = list(range(1600, 8001, 200)) + [12000]  # keV
        bins = {}
        for halldet, num_ibds_binned in num_ibds_default.items():
            num_ibds_halldet = []
            for val in num_ibds_binned[4:]:
                num_ibds_halldet.append(val)
            num_ibds[halldet] = num_ibds_halldet
            bins[halldet] = bins_halldet
    elif binning_type == 'rate-only':
        # Create rate-only binning
        num_ibds_rateonly = {}
        bins_rateonly = {}
        for halldet, num_ibds_binned in num_ibds_default.items():
            num = 0
            split_bin = num_ibds_binned[3]  # 1.4, 1.6 MeV
            next_bin = num_ibds_binned[4]  # 1.6, 1.8 MeV
            # Linear interpolation: need # from 1.5 to 1.6 MeV
            val_at_1_5 = split_bin
            val_at_1_6 = (split_bin + next_bin)/2
            # Trapezoid volume: 1/2 (b1 + b2) * width
            trap_volume = 0.5 * (val_at_1_5 + val_at_1_6) * 0.5
            num += trap_volume
            for val in num_ibds_binned[4:]:
                num += val
            num_ibds_rateonly[halldet] = [num]
            bins_rateonly[halldet] = [1500, 12000]
            num_ibds = num_ibds_rateonly
            bins = bins_rateonly
    if update_db is None:
        print('Not updating the db')
    else:
        with sqlite3.Connection(update_db) as conn:
            cursor = conn.cursor()
            for (hall, det), num_ibds_binned in num_ibds.items():
                binning = bins[halldet]
                cursor.execute('''INSERT OR REPLACE INTO num_coincidences
                    VALUES (?, ?, ?, ?, ?)''',
                    (hall, det, json.dumps(num_ibds_binned), json.dumps(binning), source))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('toy_output')
    parser.add_argument('-n', '--entry-number', type=int, default=0)
    parser.add_argument('--update-db')
    parser.add_argument('--source', default='LBNL toymc')
    parser.add_argument('--binning', default='default', choices=('default',
        'nH nominal', 'rate-only', 'nH modified 1'))
    args = parser.parse_args()

    main(args.toy_output, args.entry_number, args.update_db, args.source, args.binning)
