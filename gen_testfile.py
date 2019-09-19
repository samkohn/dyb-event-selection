'''
Generate a test ROOT file with user-defined events.

'''

import csv

from translate import (TreeBuffer, float_value, assign_value,
        fetch_value, int_value, unsigned_int_value, long_value, git_describe)

def prep_CalibStats(host_file):
    '''
    Return a tuple of (calibStats, buffer) containing the TTree object and the
    buffer used to fill its TBranches.

    '''
    buf = TreeBuffer()
    buf.triggerNumber = int_value()
    buf.detector = int_value()
    buf.timestamp_seconds = int_value()
    buf.timestamp_nanoseconds = int_value()
    buf.nHit = int_value()
    buf.charge = float_value()
    buf.fQuad = float_value()
    buf.fMax = float_value()
    buf.fPSD_t1 = float_value()
    buf.fPSD_t2 = float_value()
    buf.f2inch_maxQ = float_value()

    host_file.cd()
    event_subdir = host_file.Get('Event')
    if not bool(event_subdir):
        event_subdir = host_file.mkdir('Event')
    event_subdir.cd()
    data_subdir = event_subdir.Get('Data')
    if not bool(data_subdir):
        data_subdir = event_subdir.mkdir('Data')
    data_subdir.cd()
    calibStats = TTree('CalibStats', 'Tree at /Event/Data/CalibStats holding '
            'Data_CalibStats')
    calibStats.Branch('triggerNumber', buf.triggerNumber, 'triggerNumber/I')
    calibStats.Branch('context.mTimeStamp.mSec', buf.timestamp_seconds,
            'context.mTimeStamp.mSec/I')
    calibStats.Branch('context.mTimeStamp.mNanoSec', buf.timestamp_nanoseconds,
            'context.mTimeStamp.mNanoSec/I')
    calibStats.Branch('context.mDetId', buf.detector, 'context.mDetId/I')
    calibStats.Branch('nHit', buf.nHit, 'nHit/I')
    calibStats.Branch('NominalCharge', buf.charge, 'NominalCharge/F')
    calibStats.Branch('Quadrant', buf.fQuad, 'Quadrant/F')
    calibStats.Branch('MaxQ', buf.fMax, 'MaxQ/F')
    calibStats.Branch('time_PSD', buf.fPSD_t1, 'time_PSD/F')
    calibStats.Branch('time_PSD1', buf.fPSD_t2, 'time_PSD1/F')
    calibStats.Branch('MaxQ_2inchPMT', buf.f2inch_maxQ, 'MaxQ_2inchPMT/F')
    return calibStats, buf

def prep_AdSimple(host_file):
    '''
    Return a tuple of (adSimple, buffer) containing the TTree object and the
    buffer used to fill its TBranches.

    '''
    buf = TreeBuffer()
    buf.triggerType = unsigned_int_value()
    buf.site = int_value()
    buf.energy = float_value()
    buf.x = float_value()
    buf.y = float_value()
    buf.z = float_value()

    host_file.cd()
    event_subdir = host_file.Get('Event')
    if not bool(event_subdir):
        event_subdir = host_file.mkdir('Event')
    event_subdir.cd()
    rec_subdir = event_subdir.Get('Rec')
    if not bool(rec_subdir):
        rec_subdir = event_subdir.mkdir('Rec')
    rec_subdir.cd()
    adSimple = TTree('AdSimple', 'Tree at /Event/Rec/AdSimple holding '
            'Rec_AdSimple')
    adSimple.Branch('context.mSite', buf.site, 'context.mSite/I')
    adSimple.Branch('triggerType', buf.triggerType, 'triggerType/i')
    adSimple.Branch('energy', buf.energy, 'energy/F')
    adSimple.Branch('x', buf.x, 'x/F')
    adSimple.Branch('y', buf.y, 'y/F')
    adSimple.Branch('z', buf.z, 'z/F')
    return adSimple, buf

def assign_flasher(cs_buf, ads_buf):
    '''
    Assign the given TreeBuffers to the prototypical flasher event.

    '''
    assign_value(cs_buf.nHit, 100)
    assign_value(cs_buf.charge, 80)
    assign_value(cs_buf.fQuad, 1.234)
    assign_value(cs_buf.fMax, 2.345)
    assign_value(cs_buf.fPSD_t1, 0.1234)
    assign_value(cs_buf.fPSD_t2, 0.2345)
    assign_value(cs_buf.f2inch_maxQ, 105)
    assign_value(ads_buf.energy, 0.9)
    assign_value(ads_buf.x, -0.1)
    assign_value(ads_buf.y, 1)
    assign_value(ads_buf.z, 1.45)

def assign_prompt_like(cs_buf, ads_buf, nh):
    '''
    Assign the given TreeBuffers to the prototypical prompt-like event.

    '''
    assign_value(cs_buf.nHit, 100)
    assign_value(cs_buf.charge, 300)
    assign_value(cs_buf.fQuad, 0.1)
    assign_value(cs_buf.fMax, 0.2)
    assign_value(cs_buf.fPSD_t1, 0.9)
    assign_value(cs_buf.fPSD_t2, 0.95)
    assign_value(cs_buf.f2inch_maxQ, 1)
    assign_value(ads_buf.energy, 3.1)
    assign_value(ads_buf.x, 0.5)
    assign_value(ads_buf.y, 1.3)
    assign_value(ads_buf.z, -0.25)

def assign_delayed_like(cs_buf, ads_buf, nh):
    '''
    Assign the given TreeBuffers to the prototypical delayed-like event.

    '''
    assign_value(cs_buf.nHit, 120)
    assign_value(cs_buf.charge, 1360)
    assign_value(cs_buf.fQuad, 0.1)
    assign_value(cs_buf.fMax, 0.2)
    assign_value(cs_buf.fPSD_t1, 0.9)
    assign_value(cs_buf.fPSD_t2, 0.95)
    assign_value(cs_buf.f2inch_maxQ, 1)
    if nh:
        energy = 2.2
    else:
        energy = 8.1
    assign_value(ads_buf.energy, energy)
    assign_value(ads_buf.x, 1.5)
    assign_value(ads_buf.y, -0.3)
    assign_value(ads_buf.z, 0.65)

def assign_ADMuon(cs_buf, ads_buf):
    '''
    Assign the given TreeBuffers to the prototypical AD Muon event.

    '''
    assign_value(cs_buf.nHit, 180)
    assign_value(cs_buf.charge, 5000)
    assign_value(cs_buf.fQuad, 0.123)
    assign_value(cs_buf.fMax, 0.210)
    assign_value(cs_buf.fPSD_t1, 0.912)
    assign_value(cs_buf.fPSD_t2, 0.956)
    assign_value(cs_buf.f2inch_maxQ, 2)
    assign_value(ads_buf.energy, 30)
    assign_value(ads_buf.x, 0.1)
    assign_value(ads_buf.y, 0.2)
    assign_value(ads_buf.z, 0.3)

def assign_WSMuon(cs_buf, ads_buf):
    '''
    Assign the given TreeBuffers to the prototypical WS Muon event.

    '''
    assign_value(cs_buf.nHit, 20)
    assign_value(cs_buf.charge, 1000)
    assign_value(cs_buf.fQuad, 0)
    assign_value(cs_buf.fMax, 0)
    assign_value(cs_buf.fPSD_t1, 1)
    assign_value(cs_buf.fPSD_t2, 1)
    assign_value(cs_buf.f2inch_maxQ, 0)
    assign_value(ads_buf.energy, 0)
    assign_value(ads_buf.x, 0)
    assign_value(ads_buf.y, 0)
    assign_value(ads_buf.z, 0)

def fill_event(calibStats, cs_buf, adSimple, ads_buf,
        timestamp, triggerNumber, event_type, detector, nh):
    assign_value(cs_buf.triggerNumber, triggerNumber)
    timestamp_seconds = timestamp // int(1e9)
    timestamp_nanoseconds = timestamp % int(1e9)
    assign_value(cs_buf.timestamp_seconds, timestamp_seconds)
    assign_value(cs_buf.timestamp_nanoseconds, timestamp_nanoseconds)
    assign_value(cs_buf.detector, detector)
    assign_value(ads_buf.triggerType, 0x10001100)
    assign_value(ads_buf.site, 1)
    if event_type == 'flasher':
        assign_flasher(cs_buf, ads_buf)
    elif event_type == 'promptlike':
        assign_prompt_like(cs_buf, ads_buf, nh)
    elif event_type == 'delayedlike':
        assign_delayed_like(cs_buf, ads_buf, nh)
    elif event_type == 'admuon':
        assign_ADMuon(cs_buf, ads_buf)
    elif event_type == 'wsmuon':
        assign_WSMuon(cs_buf, ads_buf)
    calibStats.Fill()
    adSimple.Fill()

def main(infilename, outfilename, nh):
    outfile = TFile(outfilename, 'RECREATE')
    calibStats, cs_buf = prep_CalibStats(outfile)
    adSimple, ads_buf = prep_AdSimple(outfile)
    with open(infilename, 'r') as eventspec:
        event_reader = csv.DictReader(eventspec)
        for i, row in enumerate(event_reader):
            fill_event(calibStats, cs_buf, adSimple, ads_buf,
                    int(row['timestamp']), i, row['type'].strip(),
                    int(row['detector']), nh)
    outfile.Write()
    outfile.Close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default='mockevents.csv')
    parser.add_argument('-o', '--output', default='test.root')
    parser.add_argument('--nh', action='store_true')
    args = parser.parse_args()
    from ROOT import TFile, TTree
    main(args.input, args.output, args.nh)
