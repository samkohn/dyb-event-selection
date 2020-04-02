import argparse
import logging
from collections import deque
import os

from common import *
import process
import translate
import flashers
from translate import (initialize_indata_onefile as initialize, TreeBuffer,
        float_value, int_value, unsigned_int_value, long_value, assign_value)

def main_loop(events, indata, muon_ttree, singles_ttrees, coinc_ttrees, ads, debug):
    loopIndex = 0
    singles_caches = {ad: deque() for ad in ads}
    coinc_caches = {ad: (None, False) for ad in ads}
    singles_buf_copy = singles_ttrees[0][1].clone()
    coinc_buf_copy = coinc_ttrees[0][1].clone()
    while loopIndex < events:
        indata.GetEntry(loopIndex)
        if isMuon(indata):
            fill_muon_TTree(muon_ttree, indata, loopIndex)
        if isADEvent(indata):
            detector = indata.detector
            # Process singles
            add_to_singles_cache(singles_caches[detector], indata,
                    loopIndex, singles_buf_copy)
            ready_to_fills = process_singles(singles_caches)
            fill_multiple_TTrees(singles_ttrees, ready_to_fills)
            # Process coincidences
            if isCoincidenceCandidate(indata):
                ready_to_fills = update_coinc_caches(coinc_caches, indata, loopIndex,
                        coinc_buf_copy)
                fill_multiple_TTrees(coinc_ttrees, ready_to_fills)
        loopIndex += 1
    return

def update_coinc_caches(caches, indata, loopIndex, template_buf):
    buf = template_buf.clone()
    load_coinc_buf(buf, indata, loopIndex)
    cache_buf, cache_is_coinc = caches[indata.detector]
    if cache_buf is None:
        caches[indata.detector] = (buf, False)
        return []
    dt = indata.timestamp - cache_buf.timestamp[0]
    if dt < 600e3:
        cache_is_coinc = True
        current_is_coinc = True
    else:
        current_is_coinc = False
    caches[indata.detector] = (buf, current_is_coinc)
    if cache_is_coinc:
        return [cache_buf]
    else:
        return []

def isCoincidenceCandidate(indata):
    return indata.energy < 14

def isADEvent(indata):
    return (indata.detector in AD_DETECTORS
            and indata.energy > 0.7
            and int(flashers.isFlasher_nH(flashers.fID(indata.fMax, indata.fQuad),
                None, indata.f2inch_maxQ, indata.detector)) == 0)

class CacheItem:
    def __init__(self):
        self.buf = None
        # The use of sets here allows the same neighbor's data to be added
        # multiple times without saving the duplicates. This is useful if the
        # cache is traversed multiple times before finally deciding if an item
        # is good or bad.
        self.pre = set()
        self.post = set()
        self.is_good_neighbor = True
        self.is_good_single = True
        return

    def __repr__(self):
        return "<CacheItem [E={:.03}] pre={} post={} single={} neighbor={}>".format(
                self.buf.energy[0], self.pre, self.post, self.is_good_single,
                self.is_good_neighbor)

def add_to_singles_cache(cache, indata, loopIndex, template_buf):
    buf = template_buf.clone()
    load_singles_buf_noneighbors(buf, indata, loopIndex)
    item = CacheItem()
    item.buf = buf
    item.is_good_neighbor = isGoodSingleNeighbor(buf.energy[0])
    item.is_good_single = isGoodSingle(buf.energy[0])
    cache.append(item)
    return

def isGoodSingleNeighbor(energy):
    return energy < 2

def isGoodSingle(energy):
    return energy < 14

def process_singles(caches):
    good = []
    for cache in caches.values():
        revisit_current = False
        while len(cache) > 0 and not revisit_current:
            current = cache.popleft()
            t0 = current.buf.timestamp[0]
            dt = 0
            # This loop serves 2 purposes: 1. do future events cause
            # problems for the current event, and 2. does this event cause problems for
            # any future events.
            for neighbor in cache:
                dt = neighbor.buf.timestamp[0] - t0
                if dt < 600e3:
                    # Do future events cause problems for this one?
                    if neighbor.is_good_neighbor:
                        current.post.add((neighbor.buf.energy[0], dt))
                    else:
                        current.is_good_single = False
                    # Does this event cause problems for future ones?
                    if current.is_good_neighbor:
                        neighbor.pre.add((current.buf.energy[0], dt))
                    else:
                        neighbor.is_good_single = False
                else:
                    # We've come far enough that new events have no bearing on the
                    # current one. If the current event is still a new single, then
                    # save it.
                    if current.is_good_single:
                        load_singles_buf_neighbors(current.buf,
                            current.pre, current.post)
                        good.append(current.buf)
                    break
            revisit_current = dt < 600e3
            if revisit_current:
                cache.appendleft(current)
    return good


def fill_multiple_TTrees(ttrees, bufs):
    for buf in bufs:
        ad_index = buf.detector[0] - 1
        ttree, fill_buf = ttrees[ad_index]
        buf.copyTo(fill_buf)
        ttree.Fill()

def isMuon(indata):
    is_WP = (indata.detector in WP_DETECTORS
            and indata.nHit > 11)
    is_AD = (indata.detector in AD_DETECTORS
            and indata.energy > 18)
    return is_WP or is_AD

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

def load_coinc_buf(buf, indata, loopIndex):
    load_singles_buf_noneighbors(buf, indata, loopIndex)

def load_singles_buf_noneighbors(buf, indata, loopIndex):
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

def load_singles_buf_neighbors(buf, pre_iter, post_iter):
    """Load only the neighboring events into the buffer.

    Assumes pre and post are iterables of (energy, dt) tuples, not necessarily
    sorted, where dt is always positive. (dt will be deliberately negated for
    pre events.)
    """
    pre = sorted(pre_iter, key=lambda x: x[1], reverse=True)
    post = sorted(post_iter, key=lambda x: x[1])
    npre = len(pre)
    mult = npre + len(post)
    assign_value(buf.num_nearby_events, mult)
    for i, (energy, dt) in enumerate(pre):
        assign_value(buf.nearby_dt, -dt, i)
        assign_value(buf.nearby_energy, energy, i)
    for i, (energy, dt) in enumerate(post):
        assign_value(buf.nearby_dt, dt, i + npre)
        assign_value(buf.nearby_energy, energy, i + npre)
    return

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

def create_coinc_TTree(host_file):
    from ROOT import TTree
    host_file.cd()
    title = 'Coincidence events by Sam Kohn (git: {})'.format(
            translate.git_describe())
    out = TTree('coincs', title)
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
    singles_name = 'singles_ad{}_{}_{:>04}.root'.format('{}', run, fileno)
    coinc_name = 'coinc_ad{}_{}_{:>04}.root'.format('{}', run, fileno)
    muonFile = TFile(os.path.join(out_location, muon_name), 'RECREATE')
    singlesFiles = [TFile(os.path.join(out_location, singles_name.format(ad)),
        'RECREATE') for ad in ads]
    coincFiles = [TFile(os.path.join(out_location, coinc_name.format(ad)),
        'RECREATE') for ad in ads]
    return {'muon': muonFile, 'singles': singlesFiles, 'coincs': coincFiles}

def main(events, infile, out_location, run_and_file, debug):
    from ROOT import TFile
    logging.debug(events)
    logging.debug(infile)
    logging.debug(out_location)
    logging.debug(run_and_file)
    logging.debug(debug)
    run, fileno = run_and_file
    ads = process.get_ads(run)
    infile = TFile(infile, 'READ')
    calibStats, adSimple = initialize(infile)
    calibStats.AddFriend(adSimple)
    indata = process.RawFileAdapter(calibStats, run, fileno)
    outfiles = create_outfiles(out_location, run, fileno, ads)
    muon_ttree = create_muon_TTree(outfiles['muon'])
    singles_ttrees = [create_singles_TTree(f) for f in outfiles['singles']]
    coinc_ttrees = [create_coinc_TTree(f) for f in outfiles['coincs']]

    if events == -1:
        events = indata.GetEntries()
    main_loop(events, indata, muon_ttree, singles_ttrees, coinc_ttrees, ads, debug)
    outfiles['muon'].Write()
    outfiles['muon'].Close()
    [x.Write() for x in outfiles['singles']]
    [x.Close() for x in outfiles['singles']]
    [x.Write() for x in outfiles['coincs']]
    [x.Close() for x in outfiles['coincs']]
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
