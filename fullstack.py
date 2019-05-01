'''
Translate a Daya Bay recon.*.root file into a simpler .root file.

'''
from array import array
import argparse
import logging

from ROOT import TTree, TFile, TChain

import translate
import process
import rate_calculations

def main(filenames, nevents, start_event):
    if len(filenames) == 0:
        filenames = [
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0001.root",
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0002.root",
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0003.root",
                ]

    outfile = TFile('out.root', 'RECREATE')
    outdata, outdata_buf = translate.create_data_TTree(outfile)
    computed, computed_buf = process.create_computed_TTree(outfile)
    calibStats, adSimple = translate.initialize_indata(filenames)
    computed_helper = process.ProcessHelper()
    rate_helper = rate_calculations.RateHelper()

    end_event = (calibStats.GetEntries() if nevents == -1 else
            min(nevents+start_event, calibStats.GetEntries()))
    if nevents == -1 and adSimple.GetEntries() != end_event:
        print('Discrepant number of entries')
        return

    callback = rate_calculations.callback_adapter(rate_helper, start_event,
            end_event)
    for entry_number in range(start_event, end_event):
        logging.debug('Event %d', entry_number)
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
                computed_buf, computed_helper, callback)
    # After the event loop is finished, fill the remaining events from
    # the event_cache into the output TTree
    for cached_event in computed_helper.event_cache:
        translate.assign_value(cached_event.dt_next_WSMuon, -1)
        translate.assign_value(cached_event.tag_WSMuonVeto, 2)
        translate.assign_value(cached_event.dt_next_DelayedLike, -1)
        callback(cached_event)
        cached_event.copyTo(computed_buf)
        computed.Fill()

    outfile.Write()
    outfile.Close()
    rate_calculations.print_results(rate_helper)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file-list', default=None)
    parser.add_argument('-n', '--num-events', default=-1, type=int)
    parser.add_argument('-s', '--start-event', default=0, type=int)
    parser.add_argument('-d', '--debug', action='store_true')
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    infile_list = args.file_list
    infiles = []
    if infile_list is not None:
        with open(infile_list, 'r') as f:
            for line in f:
                if len(line) > 5:
                    infiles.append(line[:-1])
    main(infiles, args.num_events, args.start_event)
