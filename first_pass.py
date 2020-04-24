import argparse
import logging
from collections import deque
import os
from functools import lru_cache

from common import *
import process
import translate
import flashers
from translate import initialize_indata_onefile as initialize
from root_util import *

class RawFileAdapter():
    ATTR_LOOKUP = {
            'triggerNumber': ('triggerNumber', int),
            'timestamp_seconds': ('context.mTimeStamp.mSec', int),
            'timestamp_nanoseconds': ('context.mTimeStamp.mNanoSec',
                int),
            'detector': ('context.mDetId', int),
            'nHit': ('nHit', int),
            'charge': ('NominalCharge', float),
            'fQuad': ('Quadrant', float),
            'fMax': ('MaxQ', float),
            'fPSD_t1': ('time_PSD', float),
            'fPSD_t2': ('time_PSD1', float),
            'f2inch_maxQ': ('MaxQ_2inchPMT', float),
            'triggerType': ('triggerType', int),
            'energy': ('energy', float),
            'x': ('x', float),
            'y': ('y', float),
            'z': ('z', float),
    }
    def __init__(self, ttree_w_friend, run, fileno):
        from ROOT import TFile
        self.ttree = ttree_w_friend
        self.run = run
        self.fileno = fileno
        self.ttree.SetBranchStatus('execNumber', 1)
        # Hack to extract site number
        old_status = self.ttree.GetBranchStatus('context.mSite')
        if int(old_status) == 0:
            self.ttree.SetBranchStatus('context.mSite', 1)
        self.ttree.GetEntry(0)
        self.site = {1:1, 2:2, 4:3}[fetch_value(self.ttree, 'context.mSite', int)]
        self.ttree.SetBranchStatus('context.mSite', int(old_status))
        # End hack

    def GetEntry(self, index):
        self.ttree.GetEntry(index)
        self.__getattr__.cache_clear()

    def GetEntries(self):
        return self.ttree.GetEntries()

    def _timestamp(self):
        sec = self.timestamp_seconds
        nano = self.timestamp_nanoseconds
        return sec*1000000000 + nano

    @lru_cache(maxsize=32)
    def __getattr__(self, name):
        if name == 'timestamp':
            result = self._timestamp()
            return result
        attr_name = self.ATTR_LOOKUP.get(name, None)
        if attr_name is None:
            raise AttributeError('No attribute "{}"'.format(name))
        return fetch_value(self.ttree, attr_name[0], attr_name[1])

    def SetBranchStatus(self, *args):
        self.ttree.SetBranchStatus(*args)


def main_loop(events, indata, muon_ttree, event_ttrees, ads, debug):
    loopIndex = 0
    while loopIndex < events:
        indata.GetEntry(loopIndex)
        if isMuon(indata):
            fill_muon_TTree(muon_ttree, indata, loopIndex)
        if isADEvent(indata, ads):
            ttree, fill_buf = event_ttrees[indata.detector]
            load_adevent_buf(fill_buf, indata, loopIndex)
            ttree.Fill()
        loopIndex += 1
    return


def isADEvent(indata, ads):
    return (indata.detector in ads
            and indata.energy > 0.7
            and int(flashers.isFlasher_nH(flashers.fID(indata.fMax, indata.fQuad),
                None, indata.f2inch_maxQ, indata.detector)) == 0
            and (indata.triggerType & 0x1100) > 0)


def isMuon(indata):
    is_WP = (indata.detector in WP_DETECTORS
            and indata.nHit > 11)
    is_AD = (indata.detector in AD_DETECTORS
            and indata.energy > 18)
    good_trigger = (indata.triggerType & 0x1100) > 0
    return (is_WP or is_AD) and good_trigger

def fill_muon_TTree(ttree, indata, loopIndex):
    ttree, buf = ttree
    load_basic_TTree(buf, indata, loopIndex)
    ttree.Fill()
    return

def load_basic_TTree(buf, indata, loopIndex):
    assign_value(buf.run, indata.run)
    assign_value(buf.fileno, indata.fileno)
    assign_value(buf.site, indata.site)
    assign_value(buf.detector, indata.detector)
    assign_value(buf.loopIndex, loopIndex)
    assign_value(buf.timestamp, indata.timestamp)
    assign_value(buf.timestamp_seconds, indata.timestamp_seconds)
    assign_value(buf.timestamp_nanoseconds, indata.timestamp_nanoseconds)
    assign_value(buf.triggerNumber, indata.triggerNumber)
    assign_value(buf.triggerType, indata.triggerType)
    assign_value(buf.nHit, indata.nHit)
    assign_value(buf.energy, indata.energy)
    return

def initialize_basic_TTree(ttree, buf):
    buf.run = unsigned_int_value()
    buf.fileno = unsigned_int_value()
    buf.site = unsigned_int_value()
    buf.detector = unsigned_int_value()
    buf.loopIndex = unsigned_int_value()
    buf.timestamp = long_value()
    buf.timestamp_seconds = int_value()
    buf.timestamp_nanoseconds = int_value()
    buf.triggerNumber = int_value()
    buf.triggerType = unsigned_int_value()
    buf.nHit = int_value()
    buf.charge = float_value()
    buf.energy = float_value()

    def branch(name, typecode):
        ttree.Branch(name, getattr(buf, name), '{}/{}'.format(name,
            typecode))
        return

    branch('run', 'i')
    branch('fileno', 'i')
    branch('site', 'i')
    branch('detector', 'i')
    branch('loopIndex', 'i')
    branch('timestamp', 'L')
    branch('timestamp_seconds', 'I')
    branch('timestamp_nanoseconds', 'I')
    branch('triggerNumber', 'I')
    branch('triggerType', 'I')
    branch('nHit', 'I')
    branch('charge', 'F')
    branch('energy', 'F')
    return

def create_muon_TTree(host_file):
    from ROOT import TTree
    host_file.cd()
    title = 'Muon-like events by Sam Kohn (git: {})'.format(
            translate.git_describe())
    out = TTree('muons', title)
    buf = TreeBuffer()
    initialize_basic_TTree(out, buf)
    return out, buf

def load_adevent_buf(buf, indata, loopIndex):
    load_basic_TTree(buf, indata, loopIndex)
    assign_value(buf.fQuad, indata.fQuad)
    assign_value(buf.fMax, indata.fMax)
    assign_value(buf.fPSD_t1, indata.fPSD_t1)
    assign_value(buf.fPSD_t2, indata.fPSD_t2)
    assign_value(buf.f2inch_maxQ, indata.f2inch_maxQ)
    assign_value(buf.x, indata.x)
    assign_value(buf.y, indata.y)
    assign_value(buf.z, indata.z)
    assign_value(buf.fID, flashers.fID(indata.fMax, indata.fQuad))
    assign_value(buf.fPSD, flashers.fPSD(indata.fPSD_t1, indata.fPSD_t2))

def create_singles_TTree(host_file):
    from ROOT import TTree
    host_file.cd()
    title = 'Single events by Sam Kohn (git: {})'.format(
            translate.git_describe())
    out = TTree('singles', title)
    buf = TreeBuffer()
    initialize_basic_TTree(out, buf)
    buf.fQuad = float_value()
    buf.fMax = float_value()
    buf.fPSD_t1 = float_value()
    buf.fPSD_t2 = float_value()
    buf.f2inch_maxQ = float_value()
    buf.x = float_value()
    buf.y = float_value()
    buf.z = float_value()
    buf.fID = float_value()
    buf.fPSD = float_value()
    buf.num_nearby_events = unsigned_int_value()
    buf.nearby_dt = int_value(10)
    buf.nearby_energy = float_value(10)

    def branch(name, typecode):
        out.Branch(name, getattr(buf, name), '{}/{}'.format(name,
            typecode))
        return

    branch('fQuad', 'F')
    branch('fMax', 'F')
    branch('fPSD_t1', 'F')
    branch('fPSD_t2', 'F')
    branch('f2inch_maxQ', 'F')
    branch('x', 'F')
    branch('y', 'F')
    branch('z', 'F')
    branch('fID', 'F')
    branch('fPSD', 'F')
    branch('num_nearby_events', 'I')
    out.Branch('nearby_dt', buf.nearby_dt, 'nearby_dt[num_nearby_events]/I')
    out.Branch('nearby_energy', buf.nearby_energy,
            'nearby_energy[num_nearby_events]/F')
    return out, buf

def create_event_TTree(host_file):
    from ROOT import TTree
    host_file.cd()
    title = 'AD events by Sam Kohn (git: {})'.format(
            translate.git_describe())
    out = TTree('events', title)
    buf = TreeBuffer()
    initialize_basic_TTree(out, buf)
    buf.fQuad = float_value()
    buf.fMax = float_value()
    buf.fPSD_t1 = float_value()
    buf.fPSD_t2 = float_value()
    buf.f2inch_maxQ = float_value()
    buf.x = float_value()
    buf.y = float_value()
    buf.z = float_value()
    buf.fID = float_value()
    buf.fPSD = float_value()

    def branch(name, typecode):
        out.Branch(name, getattr(buf, name), '{}/{}'.format(name,
            typecode))
        return

    branch('fQuad', 'F')
    branch('fMax', 'F')
    branch('fPSD_t1', 'F')
    branch('fPSD_t2', 'F')
    branch('f2inch_maxQ', 'F')
    branch('x', 'F')
    branch('y', 'F')
    branch('z', 'F')
    branch('fID', 'F')
    branch('fPSD', 'F')
    return out, buf

def create_outfiles(out_location, run, fileno, ads):
    from ROOT import TFile
    muon_name = 'muons_{}_{:>04}.root'.format(run, fileno)
    events_name = 'events_ad{}_{}_{:>04}.root'.format('{}', run, fileno)
    muonFile = TFile(os.path.join(out_location, muon_name), 'RECREATE')
    eventsFiles = {ad: TFile(os.path.join(out_location, events_name.format(ad)),
        'RECREATE') for ad in ads}
    return {'muon': muonFile, 'events': eventsFiles}

def main(events, infile, out_location, run_and_file, debug):
    from ROOT import TFile
    logging.debug(events)
    logging.debug(infile)
    logging.debug(out_location)
    logging.debug(run_and_file)
    logging.debug(debug)
    run, fileno = run_and_file
    infile = TFile(infile, 'READ')
    calibStats, adSimple = initialize(infile)
    calibStats.AddFriend(adSimple)
    indata = RawFileAdapter(calibStats, run, fileno)
    ads = dets_for(indata.site, run)
    outfiles = create_outfiles(out_location, run, fileno, ads)
    muon_ttree = create_muon_TTree(outfiles['muon'])
    event_ttrees = {ad: create_event_TTree(f) for ad, f in
            outfiles['events'].items()}

    if events == -1:
        events = indata.GetEntries()
    main_loop(events, indata, muon_ttree, event_ttrees, ads, debug)
    outfiles['muon'].Write()
    outfiles['muon'].Close()
    for x in outfiles['events'].values():
        x.Write()
        x.Close()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file name')
    parser.add_argument('-o', '--output', help='Output location')
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--events', type=int, default=-1)
    parser.add_argument('-r', '--runfile', type=int, nargs=2,
        help='<run> <fileno>')
    parser.add_argument('--default', action='store_true')
    args = parser.parse_args()
    if args.input is None:
        args.input = '/global/projecta/projectdirs/dayabay/data/dropbox/p17b/lz4.skim.3/recon.Neutrino.0058043.Physics.EH1-Merged.P17B-P._0001.root'
    if args.default:
        args.runfile = [58043, 1]
        args.output = '.'
        args.debug = True
        args.events = 1000
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.events, args.input, args.output, args.runfile, args.debug)
