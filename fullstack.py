'''
Translate a Daya Bay recon.*.root file into a simpler .root file.

'''
from array import array
import argparse

from ROOT import TTree, TFile, TChain

import translate
import process

def main(filenames, nevents):
    if len(filenames) == 0:
        filenames = [
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0001.root",
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0002.root",
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0003.root",
                ]

    outfile = TFile('out3.root', 'RECREATE')
    outdata, outdata_buf = translate.create_data_TTree(outfile)
    computed, computed_buf = process.create_computed_TTree(outfile)
    calibStats, adSimple = translate.initialize_indata(filenames)
    computed_helper = process.ProcessHelper()

    n_entries = calibStats.GetEntries() if nevents == -1 else min(nevents,
            calibStats.GetEntries())
    if nevents == -1 and adSimple.GetEntries() != n_entries:
        print('Discrepant number of entries')
        return

    for entry_number in range(n_entries):
        calibStats.LoadTree(entry_number)
        calibStats.GetEntry(entry_number)
        adSimple.LoadTree(entry_number)
        adSimple.GetEntry(entry_number)

        translate.copy(outdata_buf, calibStats, adSimple)
        indata_list = []
        indata_list.append(outdata_buf.timeStamp[0])
        indata_list.append(outdata_buf.detector[0])
        indata_list.append(outdata_buf.site[0])
        indata_list.append(outdata_buf.nHit[0])
        indata_list.append(outdata_buf.charge[0])
        indata_list.append(outdata_buf.fMax[0])
        indata_list.append(outdata_buf.fQuad[0])
        indata_list.append(outdata_buf.fPSD_t1[0])
        indata_list.append(outdata_buf.fPSD_t2[0])
        indata_list.append(outdata_buf.f2inch_maxQ[0])
        indata_list.append(outdata_buf.energy[0])

        outdata.Fill()

        process.one_iteration(entry_number, indata_list, computed,
                computed_buf, computed_helper)
    # After the event loop is finished, fill the remaining events from
    # the event_cache into the output TTree
    for cached_event in computed_helper.event_cache:
        translate.assign_value(cached_event.dt_next_WSMuon, -1)
        translate.assign_value(cached_event.tag_WSMuonVeto, 2)
        translate.assign_value(cached_event.dt_next_DelayedLike, -1)
        cached_event.copyTo(computed_buf)
        computed.Fill()

    outfile.Write()
    outfile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file-list', default=None)
    parser.add_argument('-n', '--num-events', default=-1, type=int)
    args = parser.parse_args()
    infile_list = args.file_list
    infiles = []
    if infile_list is not None:
        with open(infile_list, 'r') as f:
            for line in f:
                if len(line) > 5:
                    infiles.append(line[:-1])
    main(infiles, args.num_events)
