'''
Translate a Daya Bay recon.*.root file into a simpler .root file.

'''
import argparse
import os
import subprocess
from root_util import *

def git_describe(directory=None):
    if directory is None:
        directory = os.path.dirname(__file__)
    current_dir = os.getcwd()
    try:
        os.chdir(directory)
    finally:
        output = subprocess.check_output(['git', 'describe',
            '--always', '--long', '--dirty']).strip()
        os.chdir(current_dir)
        return output.decode()



def main(filenames, nevents):
    from ROOT import TFile
    if len(filenames) == 0:
        filenames = [
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0001.root",
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0002.root",
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0003.root",
                ]


    calibStats, adSimple = initialize_indata(filenames)
    outfile = TFile('out.root', 'RECREATE')
    outdata, buf = create_data_TTree(outfile)
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
        copy(buf, calibStats, adSimple, 0, 0)
        outdata.Fill()

    outfile.Write()
    outfile.Close()

def initialize_indata_onefile(tfile_object, rec_name):
    calibStats = tfile_object.Get('/Event/Data/CalibStats')
    adSimple = tfile_object.Get(f'/Event/Rec/{rec_name}')

    activeBranches = {
            calibStats: [
                'triggerNumber',
                'context.mTimeStamp.mSec',
                'context.mTimeStamp.mNanoSec',
                'context.mDetId',
                'nHit',
                'NominalCharge',
                'Quadrant',
                'MaxQ',
                'time_PSD',
                'time_PSD1',
                'MaxQ_2inchPMT',
            ],
            adSimple: [
                'triggerType',
                'energy',
                'x', 'y', 'z',
            ]
    }
    for ttree, branch_names in activeBranches.items():
        ttree.SetBranchStatus('*', 0)
        for branch_name in branch_names:
            ttree.SetBranchStatus(branch_name, 1)
    return calibStats, adSimple

def initialize_indata(filenames):
    from ROOT import TChain
    calibStats = TChain('/Event/Data/CalibStats')
    adSimple = TChain('/Event/Rec/AdSimple')
    for filename in filenames:
        calibStats.Add(filename)
        adSimple.Add(filename)

    activeBranches = {
            calibStats: [
                'triggerNumber',
                'context*',
                'nHit',
                'NominalCharge',
                'Quadrant',
                'MaxQ',
                'time_PSD',
                'time_PSD1',
                'MaxQ_2inchPMT',
            ],
            adSimple: [
                'context*',
                'triggerType',
                'energy',
                'x', 'y', 'z',
            ]
    }
    for chain, branch_names in activeBranches.items():
        chain.SetBranchStatus('*', 0)
        for branch_name in branch_names:
            chain.SetBranchStatus(branch_name, 1)
    return calibStats, adSimple

def create_data_TTree(host_file):
    from ROOT import TTree
    host_file.cd()
    git_description = git_describe()
    outdata = TTree('raw_data', 'Daya Bay Data by Sam Kohn (git: %s)' %
            git_description)

    buf = TreeBuffer()
    buf.run = unsigned_int_value()
    buf.fileno = unsigned_int_value()
    buf.triggerNumber = int_value()
    buf.timestamp_seconds = int_value()
    buf.timestamp_nanoseconds = int_value()
    buf.timestamp = long_value()
    buf.detector = int_value()
    buf.site = int_value()
    buf.triggerType = unsigned_int_value()
    buf.nHit = int_value()
    buf.charge = float_value()
    buf.fQuad = float_value()
    buf.fMax = float_value()
    buf.fPSD_t1 = float_value()
    buf.fPSD_t2 = float_value()
    buf.f2inch_maxQ = float_value()
    buf.energy = float_value()
    buf.x = float_value()
    buf.y = float_value()
    buf.z = float_value()

    outdata.Branch('run', buf.run, 'run/i')
    outdata.Branch('fileno', buf.fileno, 'fileno/i')
    outdata.Branch('triggerNumber', buf.triggerNumber, 'triggerNumber/I')
    outdata.Branch('timestamp_seconds', buf.timestamp_seconds,
            'timestamp_seconds/I')
    outdata.Branch('timestamp_nanoseconds', buf.timestamp_nanoseconds,
            'timestamp_nanoseconds/I')
    outdata.Branch('timestamp', buf.timestamp, 'timestamp/L')
    outdata.Branch('detector', buf.detector, 'detector/I')
    outdata.Branch('site', buf.site, 'site/I')
    outdata.Branch('triggerType', buf.triggerType, 'triggerType/i')
    outdata.Branch('nHit', buf.nHit, 'nHit/I')
    outdata.Branch('charge', buf.charge, 'charge/F')
    outdata.Branch('fQuad', buf.fQuad, 'fQuad/F')
    outdata.Branch('fMax', buf.fMax, 'fMax/F')
    outdata.Branch('fPSD_t1', buf.fPSD_t1, 'fPSD_t1/F')
    outdata.Branch('fPSD_t2', buf.fPSD_t2, 'fPSD_t2/F')
    outdata.Branch('f2inch_maxQ', buf.f2inch_maxQ, 'f2inch_maxQ/F')
    outdata.Branch('energy', buf.energy, 'energy/F')
    outdata.Branch('x', buf.x, 'x/F')
    outdata.Branch('y', buf.y, 'y/F')
    outdata.Branch('z', buf.z, 'z/F')

    return outdata, buf



def copy(buf, calibStats, adSimple, run, fileno, site):
    assign_value(buf.run, run)
    assign_value(buf.fileno, fileno)
    assign_value(buf.triggerNumber, fetch_value(calibStats,
        'triggerNumber', int))
    assign_value(buf.timestamp_seconds, fetch_value(calibStats,
        'context.mTimeStamp.mSec', int))
    assign_value(buf.timestamp_nanoseconds, fetch_value(calibStats,
        'context.mTimeStamp.mNanoSec', int))
    assign_value(buf.timestamp, buf.timestamp_seconds[0]*(10**9) +
            buf.timestamp_nanoseconds[0])
    assign_value(buf.detector, fetch_value(calibStats,
        'context.mDetId', int))
    assign_value(buf.site, site)
    assign_value(buf.triggerType, fetch_value(adSimple,
        'triggerType', int))
    assign_value(buf.nHit, fetch_value(calibStats, 'nHit', int))
    assign_value(buf.charge, fetch_value(calibStats,
        'NominalCharge', float))
    assign_value(buf.fQuad, fetch_value(calibStats, 'Quadrant', float))
    assign_value(buf.fMax, fetch_value(calibStats, 'MaxQ', float))
    assign_value(buf.fPSD_t1, fetch_value(calibStats, 'time_PSD', float))
    assign_value(buf.fPSD_t2, fetch_value(calibStats, 'time_PSD1', float))
    assign_value(buf.f2inch_maxQ, fetch_value(calibStats,
        'MaxQ_2inchPMT', float))
    assign_value(buf.energy, fetch_value(adSimple, 'energy', float))
    assign_value(buf.x, fetch_value(adSimple, 'x', float))
    assign_value(buf.y, fetch_value(adSimple, 'y', float))
    assign_value(buf.z, fetch_value(adSimple, 'z', float))

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
