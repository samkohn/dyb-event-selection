"""Dump LBNL ToyMC prompt E histos to the parameters DB."""

import argparse
import json
import sys

import common
from prediction import all_ads, near_ads, far_ads

def dump_to_db(database, rows):
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.executemany('''INSERT OR REPLACE INTO num_coincidences
            VALUES (?, ?, ?, ?, ?)''', rows)
            #(hall, det, json.dumps(num_ibds_binned), json.dumps(binning), source))

def create_db(database_name):
    """Create a new database from scratch with the standard tables."""
    # for BinningId = 0, aka nominal nGd bins (0.7, 1, 1.2, 1.4, ..., 7.8, 8, 12) MeV
    bin_edges_keV = [700]
    for edge in range(1000, 8001, 200):
        bin_edges_keV.append(edge)
    bin_edges_keV.append(12000)
    rows = []
    for i, edge in enumerate(bin_edges_keV):
        rows.append((0, "nGd default (LBNL)", i, edge))
    with common.get_db(database_name) as conn:
        cursor = conn.cursor()
        cursor.execute("""
            CREATE TABLE
                reco_binnings (
                    Id INTEGER NOT NULL,
                    Name TEXT,
                    BinEdgeIndex INTEGER NOT NULL,
                    BinEdgeEnergy_keV INTEGER,
                    PRIMARY KEY(Id, BinEdgeIndex)
                )
            """
        )
        cursor.executemany("""
            INSERT OR REPLACE
            INTO
                reco_binnings
            VALUES
                (?, ?, ?, ?)
            """,
            rows
        )
        cursor.execute("""
            CREATE TABLE
                num_coincidences (
                    Hall INTEGER NOT NULL,
                    DetNo INTEGER NOT NULL,
                    BinningId INTEGER NOT NULL,
                    NumCoincidences_binned TEXT,
                    Source TEXT NOT NULL,
                    PRIMARY KEY(Hall, DetNo, BinningId, Source),
                    FOREIGN KEY(BinningId)
                        REFERENCES reco_binnings(Id)
                )
            """
        )
    return


def multiple(infilename, entries, update_db, sources, binning_type, nominal_near,
        nominal_far, stage_label):
    rows = []
    for entry_number, source in zip(entries, sources):
        entry_rows = single(
            infilename, entry_number, None, source, binning_type, nominal_near,
            nominal_far, stage_label
        )
        rows.extend(entry_rows)
    if update_db is None:
        return rows
    else:
        dump_to_db(update_db, rows)
        return rows


def single(infilename, entry_number, update_db, source, binning_type, nominal_near,
        nominal_far, stage_label):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    host_ttree = infile.Get('tr')
    host_ttree.GetEntry(entry_number)
    stages = {'6ad': 1, '8ad': 2, '7ad': 3}
    stage = stages[stage_label]
    all_hists = {}
    for i, halldet in enumerate(all_ads):
        if (
            (nominal_near and halldet in near_ads)
            or (nominal_far and halldet in far_ads)
        ):
            name = f'h_nominal_stage{stage}_ad{i+1}'
            all_hists[halldet] = infile.Get(name)
        else:
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
        binning_id = 0
    elif binning_type == 'default minus lowest':
        binning_id = 1
        num_ibds = {}
        bins = bins_default
        for halldet, ibds_halldet in num_ibds_default.items():
            new_ibds_halldet = ibds_halldet.copy()
            new_ibds_halldet[0] = 0
            new_ibds_halldet[1] = 0
            new_ibds_halldet[2] = 0
            new_ibds_halldet[3] = 0
            num_ibds[halldet] = new_ibds_halldet
    elif binning_type == 'nH nominal':
        binning_id = 2
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
        binning_id = 3
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
        binning_id = 4
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
    rows = []
    for (hall, det), num_ibds_binned in num_ibds.items():
        binning = bins[halldet]
        row = (hall, det, binning_id, json.dumps(num_ibds_binned), source)
        rows.append(row)
    if update_db is None:
        #print('Not updating the db')
        return rows
    else:
        dump_to_db(update_db, rows)
        return rows
        #with sqlite3.Connection(update_db) as conn:
            #cursor = conn.cursor()
            #for (hall, det), num_ibds_binned in num_ibds.items():
                #binning = bins[halldet]
                #cursor.execute('''INSERT OR REPLACE INTO num_coincidences
                    #VALUES (?, ?, ?, ?, ?)''',
                    #(hall, det, json.dumps(num_ibds_binned), json.dumps(binning), source))

main = single

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('toy_output')
    parser.add_argument('-n', '--entry-number', type=int, default=0)
    parser.add_argument('--update-db')
    parser.add_argument('--source', default='LBNL toymc')
    parser.add_argument('--binning', default='default', choices=('default',
        'default minus lowest',
        'nH nominal', 'rate-only', 'nH modified 1'))
    parser.add_argument('--nominal-near', action='store_true')
    parser.add_argument('--nominal-far', action='store_true')
    parser.add_argument('--stage', choices=('6ad', '7ad', '8ad'))
    parser.add_argument('--create-db', action='store_true')
    args = parser.parse_args()

    if args.create_db:
        print(f"Creating new db at location \"{args.toy_output}\"")
        create_db(args.toy_output)
        sys.exit(0)

    rows = main(args.toy_output, args.entry_number, args.update_db, args.source, args.binning,
            args.nominal_near, args.nominal_far, args.stage)
    if not args.update_db:
        print(rows)
