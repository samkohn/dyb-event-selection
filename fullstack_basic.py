'''
Translate a Daya Bay recon.*.root file into a simpler .root file.

'''
from array import array
import argparse
import logging
import json
import re
import subprocess

def extract_run_fileno(filename):
    run = int(re.search('\d{7}', filename).group(0))
    fileno = int(re.search('_\d{4}', filename).group(0)[1:])
    return (run, fileno)

def main(filename, nevents, start_event, site):
    from ROOT import TTree, TFile, TChain

    import translate
    import process_basic as process
    if filename is None:
        filename = "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0001.root"
        site = 1

    run, fileno = extract_run_fileno(filename)

    outfile = TFile('basic_%d_%05d.root' % (run, fileno), 'RECREATE',
            'DYB Run %d file %d, git %s' % (run, fileno,
                translate.git_describe()))
    #outdata, outdata_buf = translate.create_data_TTree(outfile)
    computed, computed_buf = process.create_computed_TTree('slimmed', outfile,
            'Slimmed down dataset with muons, high-energy events and '
            'coincidences (git: %s)')
    calibStats, adSimple = translate.initialize_indata([filename])
    computed_helper = process.ProcessHelper()
    computed_helper.run = run
    computed_helper.fileno = fileno
    end_event = (calibStats.GetEntries() if nevents == -1 else
            min(nevents+start_event, calibStats.GetEntries()))
    if nevents == -1 and adSimple.GetEntries() != end_event:
        print('Discrepant number of entries')
        return

    for entry_number in range(start_event, end_event):
        logging.debug('Event %d', entry_number)
        calibStats.LoadTree(entry_number)
        calibStats.GetEntry(entry_number)
        adSimple.LoadTree(entry_number)
        adSimple.GetEntry(entry_number)

        translate.copy(computed_buf, calibStats, adSimple, run, fileno)

        process.one_iteration(entry_number, computed,
                computed_buf, computed_helper)
    # After the event loop is finished, fill the remaining events from
    # the event_cache into the output TTree
    process.finish_emptying_cache(computed, computed_buf,
            computed_helper.event_cache)
    outfile.Write()
    outfile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', default=None)
    parser.add_argument('-n', '--num-events', default=-1, type=int)
    parser.add_argument('-s', '--start-event', default=0, type=int)
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('--site', type=int)
    parser.add_argument('--lineno', type=int, default=0)
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    try:
        main(args.infile, args.num_events, args.start_event, args.site)
    except Exception as e:
        logging.exception(e)
        subprocess.check_output(['touch',
            '/global/homes/s/skohn/dyb-event-selection-production/'
            'progress/__possible_error_%d__' % args.lineno])
