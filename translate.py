'''
Translate a Daya Bay recon.*.root file into a simpler .root file.

'''
from array import array
from ROOT import TTree, TFile

class TreeBuffer():
    pass

def int_value():
    return array('i', [0])
def unsigned_int_value():
    return array('I', [0])
def float_value():
    return array('f', [0])

def assign_value(buf_value, ttree, branch_name, type_cast=None, index=0):
    try:
        new_value = getattr(ttree, branch_name)
        buf_value[index] = new_value
    except TypeError:
        new_value = fetch_annoying_value(ttree, branch_name, type_cast)
        buf_value[index] = new_value

def fetch_annoying_value(ttree, branch_name, type_cast):
    return type_cast(ttree.GetBranch(branch_name).GetValue(0, 0))

def main():
    filename = "/project/projectdirs/dayabay/data/exp/dayabay/2015/p14b/Neutrino/0107/recon.Neutrino.0048415.Physics.EH1-Merged.P14B-P._0278.root"
    infile = TFile(filename)
    calibStats = infile.Get('/Event/Data/CalibStats')
    adSimple = infile.Get('/Event/Rec/AdSimple')

    outfile = TFile('out.root', 'RECREATE')
    outdata = TTree('data', 'Daya Bay Data by Sam Kohn')

    buf = TreeBuffer()
    buf.triggerNumber = int_value()
    buf.timeStamp_seconds = int_value()
    buf.timeStamp_nanoseconds = int_value()
    buf.detector = int_value()
    buf.site = int_value()
    buf.triggerType = unsigned_int_value()
    buf.nHit = int_value()
    buf.fQuad = float_value()
    buf.fMax = float_value()
    buf.fPSD_t1 = float_value()
    buf.fPSD_t2 = float_value()
    buf.f2inch_maxQ = float_value()
    buf.energy = float_value()
    buf.x = float_value()
    buf.y = float_value()
    buf.z = float_value()



    outdata.Branch('triggerNumber', buf.triggerNumber, 'triggerNumber/I')
    outdata.Branch('timeStamp_seconds', buf.timeStamp_seconds,
            'timeStamp_seconds/I')
    outdata.Branch('timestamp_nanoseconds', buf.timeStamp_nanoseconds,
            'timeStamp_nanoseconds/I')
    outdata.Branch('detector', buf.detector, 'detector/I')
    outdata.Branch('site', buf.site, 'site/I')
    outdata.Branch('triggerType', buf.triggerType, 'triggerType/i')
    outdata.Branch('nHit', buf.nHit, 'nHit/I')
    outdata.Branch('fQuad', buf.fQuad, 'fQuad/F')
    outdata.Branch('fMax', buf.fMax, 'fMax/F')
    outdata.Branch('fPSD_t1', buf.fPSD_t1, 'fPSD_t1/F')
    outdata.Branch('fPSD_t2', buf.fPSD_t2, 'fPSD_t2/F')
    outdata.Branch('f2inch_maxQ', buf.f2inch_maxQ, 'f2inch_maxQ/F')
    outdata.Branch('energy', buf.energy, 'energy/F')
    outdata.Branch('x', buf.x, 'x/F')
    outdata.Branch('y', buf.y, 'y/F')
    outdata.Branch('z', buf.z, 'z/F')

    n_entries = calibStats.GetEntries()
    if adSimple.GetEntries() != n_entries:
        print('Discrepant number of entries')
        return

    for entry_number in range(n_entries):
        calibStats.LoadTree(entry_number)
        calibStats.GetEntry(entry_number)
        adSimple.LoadTree(entry_number)
        adSimple.GetEntry(entry_number)

        assign_value(buf.triggerNumber, calibStats, 'triggerNumber')
        assign_value(buf.timeStamp_seconds, calibStats, 'context.mTimeStamp.mSec', int)
        assign_value(buf.timeStamp_nanoseconds, calibStats, 'context.mTimeStamp.mNanoSec', int)
        assign_value(buf.detector, calibStats, 'context.mDetId', int)
        assign_value(buf.site, adSimple, 'context.mSite', int)
        assign_value(buf.triggerType, adSimple, 'triggerType', int)
        assign_value(buf.nHit, calibStats, 'nHit')
        assign_value(buf.fQuad, calibStats, 'Quadrant')
        assign_value(buf.fMax, calibStats, 'MaxQ')
        assign_value(buf.fPSD_t1, calibStats, 'time_PSD')
        assign_value(buf.fPSD_t2, calibStats, 'time_PSD1')
        assign_value(buf.f2inch_maxQ, calibStats, 'MaxQ_2inchPMT')
        assign_value(buf.energy, adSimple, 'energy', float)
        assign_value(buf.x, adSimple, 'x', float)
        assign_value(buf.y, adSimple, 'y', float)
        assign_value(buf.z, adSimple, 'z', float)
        outdata.Fill()

    outfile.Write()
    outfile.Close()
    infile.Close()

if __name__ == '__main__':
    main()
