import argparse
import math
import random

import delayeds
from root_util import assign_value
from first_pass_adtime import create_event_TTree

def is_single(ttree_event, energy):
    'Applies the single event criteria to the loaded entry in this TTree.'
    return (
            energy < 12 and energy > 1.5
            and ttree_event.multiplicity == 1
            and ttree_event.dt_cluster_to_prev_ADevent >
                delayeds._NH_THU_MAX_TIME)

def main(infilename, outfile, ttree_name):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    outfile = ROOT.TFile(outfile, 'RECREATE')
    singles, singles_buf = create_event_TTree(outfile)
    singles.SetTitle('Singles sample')
    singles.SetName('singles')
    outfile.cd()
    computed = infile.Get(ttree_name)
    entries = computed.GetEntries()
    index = 0
    while index < entries:
        # Increment until we get a good single
        computed.GetEntry(index)
        energy = computed.energy[0]
        if is_single(computed, energy):
            assign_value(singles_buf.run, computed.run)
            assign_value(singles_buf.fileno, computed.fileno)
            assign_value(singles_buf.site, computed.site)
            assign_value(singles_buf.loopIndex, computed.loopIndex[0])
            assign_value(singles_buf.timestamp, computed.timestamp[0])
            assign_value(singles_buf.timestamp_seconds, computed.timestamp_seconds[0])
            assign_value(singles_buf.timestamp_nanoseconds,
                    computed.timestamp_nanoseconds[0])
            assign_value(singles_buf.detector, computed.detector[0])
            assign_value(singles_buf.triggerNumber, computed.triggerNumber[0])
            assign_value(singles_buf.triggerType, computed.triggerType[0])
            assign_value(singles_buf.nHit, computed.nHit[0])
            assign_value(singles_buf.charge, computed.charge[0])
            assign_value(singles_buf.fQuad, computed.fQuad[0])
            assign_value(singles_buf.fMax, computed.fMax[0])
            assign_value(singles_buf.fPSD_t1, computed.fPSD_t1[0])
            assign_value(singles_buf.fPSD_t2, computed.fPSD_t2[0])
            assign_value(singles_buf.f2inch_maxQ, computed.f2inch_maxQ[0])
            assign_value(singles_buf.energy, computed.energy[0])
            assign_value(singles_buf.energy_AdTime, computed.energy_AdTime[0])
            assign_value(singles_buf.x, computed.x[0])
            assign_value(singles_buf.y, computed.y[0])
            assign_value(singles_buf.z, computed.z[0])
            assign_value(singles_buf.x_AdTime, computed.x_AdTime[0])
            assign_value(singles_buf.y_AdTime, computed.y_AdTime[0])
            assign_value(singles_buf.z_AdTime, computed.z_AdTime[0])
            singles.Fill()
        index += 1
    outfile.Write()
    outfile.Close()
    infile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('--ttree-name', default='ad_events')
    args = parser.parse_args()

    main(args.infile, args.outfile, args.ttree_name)
