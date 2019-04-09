'''
Compute the livetime inefficiency due to muons.

'''
from __future__ import print_function

from ROOT import TFile, TTree

from translate import fetch_value
import muons

def main():
    filename = 'out.root'
    infile = TFile(filename, 'UPDATE')
    indata = infile.Get('data')
    incomputed = infile.Get('computed')

    total_nonvetoed_livetime = {1:0, 2:0}  # by AD
    start_time = 0
    end_time = 0

    for event_number in indata.GetEntries():
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)
        incomputed.LoadTree(event_number)
        indata.GetEntry(event_number)

        timestamp = fetch_value(indata, 'timeStamp', int)
        detector = fetch_value(inddata, 'detector', int)
        isWSMuon = fetch_value(incomputed, 'tag_WSMuon', bool)
        isADMuon = fetch_value(incomputed, 'tag_ADMuon', bool)
        isShowerMuon = fetch_value(incomputed, 'tag_ShowerMuon', bool)
        dt_last_WSMuon = fetch_value(incomputed, 'dt_previous_WSMuon', int)
        dt_last_ADMuon = fetch_value(incomputed, 'dt_previous_ADMuon', int)
        dt_last_ShowerMuon = fetch_value(incomputed,
                'dt_previous_ShowerMuon', int)

        if event_number == 0:
            start_time = timestamp
        if event_number == indata.GetEntries() - 1:
            end_time = timestamp

        if isWSMuon or isADMuon or isShowerMuon:
            dt_past_WSVeto = max(0, dt_last_WSMuon -
                    muons._WSMUON_VETO_LAST_NS)
            dt_past_ADVeto = max(0, dt_last_ADMuon -
                    muons._ADMUON_VETO_LAST_NS)
            dt_past_ShowerVeto = max(0, dt_last_ShowerMuon -
                    muons._SHOWER_MUON_VETO_LAST_NS)

            new_livetime = min(dt_past_WSVeto, dt_past_ADVeto,
                    dt_past_ShowerVeto)

            if isWSMuon:
                new_livetime -= muons._WSMUON_VETO_NEXT_NS
                new_livetime = max(0, new_livetime)
            total_nonvetoed_livetime[detector] += new_livetime
    print('total DAQ livetime:')
    daq_livetime = end_time - start_time
    print(daq_livetime)
    print('total nonvetoed livetime:')
    print(total_nonvetoed_livetime)
    efficiency = {n: nonvetoed/daq_livetime for n, nonvetoed in
            total_nonvetoed_livetime.items()}

if __name__ == '__main__':
    main()

