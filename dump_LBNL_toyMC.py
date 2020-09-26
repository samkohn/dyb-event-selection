"""Dump LBNL ToyMC prompt E histos to the parameters DB."""

import argparse
import json
import sqlite3

from prediction import all_ads

def main(infilename, entry_number, update_db, source):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    host_ttree = infile.Get('tr')
    host_ttree.GetEntry(entry_number)
    all_hists = {}
    for i, halldet in enumerate(all_ads):
        for stage in (2,):
            name = f'h_stage{stage}_ad{i+1}'
            all_hists[halldet] = getattr(host_ttree, name)
    num_ibds = {}
    bins = {}
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
        num_ibds[halldet] = num_ibds_halldet
        bins[halldet] = bins_halldet
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
    args = parser.parse_args()

    main(args.toy_output, args.entry_number, args.update_db, args.source)

