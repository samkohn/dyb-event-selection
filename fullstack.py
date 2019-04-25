'''
Translate a Daya Bay recon.*.root file into a simpler .root file.

'''
from array import array
import argparse

from ROOT import TTree, TFile, TChain

class TreeBuffer(object):
    def copyTo(self, other):
        '''
        Copy the values in this TreeBuffer to an existing TreeBuffer set
        up with the same array attributes.

        '''
        for attr, value in self.__dict__.items():
            otherArray = getattr(other, attr)
            for i, entry in enumerate(value):
                otherArray[i] = entry

    def clone(self):
        '''
        Create a new independent TreeBuffer with the same array
        attributes and values.

        '''
        new = TreeBuffer()
        for attr, value in self.__dict__.items():
            setattr(new, attr, value[:])
        return new

    def clone_type(self):
        '''
        Create a new TreeBuffer with the same array attributes
        containing values of 0.

        '''
        new = TreeBuffer()
        for attr, value in self.__dict__.items():
            setattr(new, attr, array(value.typecode,
                [0]*len(value)))
        return new

def int_value(length=1):
    return array('i', [0]*length)
def long_value(length=1):
    return array('l', [0]*length)
def unsigned_int_value(length=1):
    return array('I', [0]*length)
def float_value(length=1):
    return array('f', [0]*length)

def fetch_value(ttree, branch_name, type_cast=None):
    try:
        new_value = type_cast(getattr(ttree, branch_name))
    except TypeError:
        new_value = fetch_annoying_value(ttree, branch_name, type_cast)
    return new_value

def assign_value(buf_value, new_value, index=0):
    buf_value[index] = new_value

def fetch_annoying_value(ttree, branch_name, type_cast):
    return type_cast(ttree.GetBranch(branch_name).GetValue(0, 0))

def main(filenames, nevents):
    if len(filenames) == 0:
        filenames = [
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0001.root",
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0002.root",
                "/project/projectdirs/dayabay/data/exp/dayabay/2015/p15a/Neutrino/0126/recon.Neutrino.0050958.Physics.EH1-Merged.P15A-P._0003.root",
                ]

    calibStats = TChain('/Event/Data/CalibStats')
    adSimple = TChain('/Event/Rec/AdSimple')
    for filename in filenames:
        calibStats.Add(filename)
        adSimple.Add(filename)

    outfile = TFile('out.root', 'RECREATE')
    outdata = TTree('data', 'Daya Bay Data by Sam Kohn')

    buf = TreeBuffer()
    buf.triggerNumber = int_value()
    buf.timeStamp_seconds = int_value()
    buf.timeStamp_nanoseconds = int_value()
    buf.timeStamp = long_value()
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


    outdata.Branch('triggerNumber', buf.triggerNumber, 'triggerNumber/I')
    outdata.Branch('timeStamp_seconds', buf.timeStamp_seconds,
            'timeStamp_seconds/I')
    outdata.Branch('timestamp_nanoseconds', buf.timeStamp_nanoseconds,
            'timeStamp_nanoseconds/I')
    outdata.Branch('timeStamp', buf.timeStamp, 'timeStamp/L')
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

        assign_value(buf.triggerNumber, fetch_value(calibStats,
            'triggerNumber', int))
        assign_value(buf.timeStamp_seconds, fetch_value(calibStats,
            'context.mTimeStamp.mSec', int))
        assign_value(buf.timeStamp_nanoseconds, fetch_value(calibStats,
            'context.mTimeStamp.mNanoSec', int))
        assign_value(buf.timeStamp, buf.timeStamp_seconds[0]*(10**9) +
                buf.timeStamp_nanoseconds[0])
        assign_value(buf.detector, fetch_value(calibStats,
            'context.mDetId', int))
        assign_value(buf.site, fetch_value(adSimple, 'context.mSite',
            int))
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
        outdata.Fill()

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