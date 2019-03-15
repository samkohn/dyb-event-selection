'''
Translate a Daya Bay recon.*.root file into a simpler .root file.

'''
from array import array
from ROOT import TTree, TFile

class TreeBuffer():
    pass

def main():
    filename = "/project/projectdirs/dayabay/data/exp/dayabay/2015/p14b/Neutrino/0107/recon.Neutrino.0048415.Physics.EH1-Merged.P14B-P._0278.root"
    infile = TFile(filename)
    calibStats = infile.Get('/Event/Data/CalibStats')
    adSimple = infile.Get('/Event/Rec/AdSimple')

    outfile = TFile('out.root', 'RECREATE')
    outdata = TTree('data', 'Daya Bay Data by Sam Kohn')

    buf = TreeBuffer()
    buf.triggerNumber = array('L', [0])

    outdata.Branch('triggerNumber', buf.triggerNumber, 'triggerNumber/l')

    n_entries = calibStats.GetEntries()
    if adSimple.GetEntries() != n_entries:
        print('Discrepant number of entries')
        return

    entry_number = 0
    calibStats.LoadTree(entry_number)
    calibStats.GetEntry(entry_number)
    adSimple.LoadTree(entry_number)

    buf.triggerNumber[0] = calibStats.triggerNumber

    outdata.Fill()
    outfile.Write()
    outfile.Close()
    infile.Close()

if __name__ == '__main__':
    main()
