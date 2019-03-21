'''
Process the basic data by computing various cuts, quantities and tags.

'''
from __future__ import print_function
from ROOT import TFile, TTree
from flashers import fID, fPSD, isFlasher
from muons import isWSMuon, isADMuon
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
    buf.tag_WSMuon = unsigned_int_value()
    buf.tag_ADMuon = unsigned_int_value()
    buf.tag_ShowerMuon = unsigned_int_value()

    outdata.Branch('fID', buf.fID, 'fID/F')
    outdata.Branch('fPSD', buf.fPSD, 'fPSD/F')
    outdata.Branch('tag_flasher', buf.tag_flasher, 'tag_flasher/i')
    outdata.Branch('tag_WSMuon', buf.tag_WSMuon, 'tag_WSMuon/i')
    outdata.Branch('tag_ADMuon', buf.tag_ADMuon, 'tag_ADMuon/i')
    outdata.Branch('tag_ShowerMuon', buf.tag_ShowerMuon,
            'tag_ShowerMuon/i')

    for event_number in xrange(indata.GetEntries()):
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)

        detector = fetch_value(indata, 'detector', int)
        nHit = fetch_value(indata, 'nHit', int)
        charge = fetch_value(indata, 'charge', float)
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
            event_isWSMuon = isWSMuon(detector, nHit)
            assign_value(buf.tag_WSMuon, event_isWSMuon)
            event_isADMuon = isADMuon(charge)
            assign_value(buf.tag_ADMuon, event_isADMuon)
            event_isShowerMuon = isShowerMuon(charge)
            assign_value(buf.tag_ShowerMuon, event_isShowerMuon)
        except:
            print(fMax, fQuad, fPSD_t1, fPSD_t2)
            raise

        outdata.Fill()

    infile.Write()
    infile.Close()

if __name__ == '__main__':
    main()
