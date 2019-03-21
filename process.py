'''
Process the basic data by computing various cuts, quantities and tags.

'''
from __future__ import print_function
from ROOT import TFile, TTree
from flashers import fID, fPSD, isFlasher
from translate import (TreeBuffer, float_value, assign_value,
        fetch_value, unsigned_int_value)

def main():

    filename = 'out.root'
    infile = TFile(filename, 'UPDATE')
    indata = infile.Get('data')
    outdata = TTree('computed', 'Computed quantities by Sam Kohn')

    buf = TreeBuffer()
    buf.fID = float_value()
    buf.fPSD = float_value()
    buf.tag_flasher = unsigned_int_value()

    outdata.Branch('fID', buf.fID, 'fID/F')
    outdata.Branch('fPSD', buf.fPSD, 'fPSD/F')
    outdata.Branch('tag_flasher', buf.tag_flasher, 'tag_flasher/i')

    for event_number in xrange(indata.GetEntries()):
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)

        fMax = fetch_value(indata, 'fMax', float)
        fQuad = fetch_value(indata, 'fQuad', float)
        fPSD_t1 = fetch_value(indata, 'fPSD_t1', float)
        fPSD_t2 = fetch_value(indata, 'fPSD_t2', float)

        try:
            event_fID = fID(fMax, fQuad)
            assign_value(buf.fID, event_fID)
            event_fPSD = fPSD(fPSD_t1, fPSD_t2)
            assign_value(buf.fPSD, event_fPSD)
            event_isFlasher = isFlasher(event_fID, event_fPSD)
            assign_value(buf.tag_flasher, event_isFlasher)
        except:
            print(fMax, fQuad, fPSD_t1, fPSD_t2)
            raise

        outdata.Fill()

    infile.Write()
    infile.Close()

if __name__ == '__main__':
    main()
