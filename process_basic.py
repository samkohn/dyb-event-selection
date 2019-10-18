'''
Extract from the bare data stream the following:

    - WP events with NHIT > 12
    - AD events with energy > 12 MeV
    - AD events that survive the Ellipse and 2-inch flasher cuts
    - AD events with the following coincidence criteria:
        - energy > 0.5 MeV, and
        - within 1ms of another AD event (whose energy is > 0.5 MeV)

'''
from __future__ import print_function

from collections import deque
import time
import argparse
import logging
import math

from ROOT import TFile, TTree
from common import *
from flashers import fID, fPSD
import flashers
import muons
import adevent
from translate import (TreeBuffer, float_value, assign_value,
        fetch_value, int_value, unsigned_int_value, long_value, git_describe)

def done_with_cache(buf):
    '''
    Test whether an event should be removed from the event cache, based
    on whether certain values have been assigned to it.

    '''
    # An event is done if it is a flasher, or if not, then if the next
    # WSMuon and the next delayed-like events have been processed. (If
    # the event is not an AD event then the delayed-like event does not
    # have to have been processed.)
    found_next_ADevent = buf.dt_next_ADevent[0] != -1
    detector = buf.detector[0]
    return found_next_ADevent or detector not in AD_DETECTORS

def create_computed_TTree(name, host_file, title=None):
    git_description = git_describe()
    if title is None:
        title = ('Computed quantities by Sam Kohn (git: %s)' %
            git_description)
    else:
        try:
            title = title % git_description
        except:
            pass
    host_file.cd()
    outdata = TTree(name, title)
    # Initialize the "buffer" used to fill values into the TTree
    fill_buf = TreeBuffer()
    fill_buf.noTree_loopIndex = int_value()
    fill_buf.timestamp = long_value()
    fill_buf.timestamp_seconds = int_value()
    fill_buf.timestamp_nanoseconds = int_value()
    fill_buf.detector = int_value()
    fill_buf.site = int_value()
    fill_buf.run = unsigned_int_value()
    fill_buf.fileno = unsigned_int_value()
    fill_buf.triggerNumber = int_value()
    fill_buf.detector = int_value()
    fill_buf.site = int_value()
    fill_buf.triggerType = unsigned_int_value()
    fill_buf.nHit = int_value()
    fill_buf.charge = float_value()
    fill_buf.fQuad = float_value()
    fill_buf.fMax = float_value()
    fill_buf.fPSD_t1 = float_value()
    fill_buf.fPSD_t2 = float_value()
    fill_buf.f2inch_maxQ = float_value()
    fill_buf.energy = float_value()
    fill_buf.x = float_value()
    fill_buf.y = float_value()
    fill_buf.z = float_value()
    fill_buf.fID = float_value()
    fill_buf.fPSD = float_value()
    fill_buf.tag_flasher = unsigned_int_value()
    fill_buf.tag_ADevent = unsigned_int_value()
    fill_buf.tag_WPevent = unsigned_int_value()
    fill_buf.tag_coincidence = unsigned_int_value()
    fill_buf.dt_last_ADevent = long_value()
    fill_buf.dt_next_ADevent = long_value()

    # Initialize the new TTree so that each TBranch reads from the
    # appropriate TreeBuffer attribute
    outdata.Branch('run', fill_buf.run, 'run/i')
    outdata.Branch('fileno', fill_buf.fileno, 'fileno/i')
    outdata.Branch('triggerNumber', fill_buf.triggerNumber, 'triggerNumber/I')
    outdata.Branch('timestamp_seconds', fill_buf.timestamp_seconds,
            'timestamp_seconds/I')
    outdata.Branch('timestamp_nanoseconds', fill_buf.timestamp_nanoseconds,
            'timestamp_nanoseconds/I')
    outdata.Branch('timestamp', fill_buf.timestamp, 'timestamp/L')
    outdata.Branch('detector', fill_buf.detector, 'detector/I')
    outdata.Branch('site', fill_buf.site, 'site/I')
    outdata.Branch('triggerType', fill_buf.triggerType, 'triggerType/i')
    outdata.Branch('nHit', fill_buf.nHit, 'nHit/I')
    outdata.Branch('charge', fill_buf.charge, 'charge/F')
    outdata.Branch('fQuad', fill_buf.fQuad, 'fQuad/F')
    outdata.Branch('fMax', fill_buf.fMax, 'fMax/F')
    outdata.Branch('fPSD_t1', fill_buf.fPSD_t1, 'fPSD_t1/F')
    outdata.Branch('fPSD_t2', fill_buf.fPSD_t2, 'fPSD_t2/F')
    outdata.Branch('f2inch_maxQ', fill_buf.f2inch_maxQ, 'f2inch_maxQ/F')
    outdata.Branch('energy', fill_buf.energy, 'energy/F')
    outdata.Branch('x', fill_buf.x, 'x/F')
    outdata.Branch('y', fill_buf.y, 'y/F')
    outdata.Branch('z', fill_buf.z, 'z/F')
    outdata.Branch('fID', fill_buf.fID, 'fID/F')
    outdata.Branch('fPSD', fill_buf.fPSD, 'fPSD/F')
    outdata.Branch('tag_flasher', fill_buf.tag_flasher, 'tag_flasher/i')
    outdata.Branch('tag_ADevent', fill_buf.tag_ADevent, 'tag_ADevent/i')
    outdata.Branch('tag_WPevent', fill_buf.tag_WPevent, 'tag_WPevent/i')
    outdata.Branch('tag_coincidence', fill_buf.tag_coincidence, 'tag_coincidence/i')
    outdata.Branch('dt_last_ADevent', fill_buf.dt_last_ADevent,
            'dt_last_ADevent/L')
    outdata.Branch('dt_next_ADevent', fill_buf.dt_next_ADevent,
            'dt_next_ADevent/L')
    return outdata, fill_buf

class ProcessHelper(object):
    def __init__(self):
        self.run = 0
        self.fileno = 0
        # Initialize the machinery for dt_next_* and dt_previous_*
        # attributes. The event_cache stores the TreeBuffer for each event
        # until the next WSMuon and next DelayedLike have been found, thus
        # preserving the event order between the input and output TTrees.
        self.event_cache = deque()
        self.first_index_without_next_ADevent = 0
        self.go_through_whole_cache = True
        self.last_ADevent_time = {n:-1 for n in range(9)}

def main(entries, debug):

    filename = 'out.root'
    infile = TFile(filename, 'UPDATE')
    indata = infile.Get('raw_data')
    outdata, fill_buf = create_computed_TTree('data', infile)

    helper = ProcessHelper()
    if entries == -1:
        entries = xrange(indata.GetEntries())
    else:
        entries = xrange(entries)
    if debug:
        entries = xrange(10000)
    def callback(event_cache):
        logging.debug(event_cache.noTree_loopIndex[0])
    for event_number in entries:
        # Load the current event in the input TTree
        indata.LoadTree(event_number)
        indata.GetEntry(event_number)
        fetch_indata(indata, fill_buf)

        one_iteration(event_number, outdata, fill_buf, helper, callback)
    finish_emptying_cache(outdata, fill_buf, helper.event_cache, callback)
    outdata.Write()
    infile.Close()

def finish_emptying_cache(outdata, fill_buf, cache, callback=lambda e:None):
    for cached_event in cache:
        e = cached_event
        coincidence = adevent.hasCoincidence(
                e.tag_ADevent[0],
                e.dt_last_ADevent[0],
                e.dt_next_ADevent[0],
                e.tag_flasher[0])
        high_energy = adevent.hasHighEnergy(e.tag_ADevent[0], e.energy[0])
        logging.debug('cached event: %s', cached_event)
        logging.debug('hasCoincidence: %s', coincidence)
        logging.debug('hasHighEnergy: %s', high_energy)
        logging.debug('isWPevent: %s', cached_event.tag_WPevent[0])
        assign_value(cached_event.tag_coincidence, coincidence)
        callback(cached_event)
        if coincidence or high_energy or cached_event.tag_WPevent[0]:
            cached_event.copyTo(fill_buf)
            outdata.Fill()


def fetch_indata(indata, buf):
    # Fetch the necessary values from the input TTree
    assign_value(buf.timestamp, fetch_value(indata, 'timestamp', int))
    assign_value(buf.timestamp_seconds, fetch_value(indata, 'timestamp_seconds',
        int))
    assign_value(buf.timestamp_nanoseconds, fetch_value(indata,
            'timestamp_nanoseconds', int))
    assign_value(buf.triggerType, fetch_value(indata, 'triggerType', int))
    assign_value(buf.detector, fetch_value(indata, 'detector', int))
    assign_value(buf.site, fetch_value(indata, 'site', int))
    assign_value(buf.nHit, fetch_value(indata, 'nHit', int))
    assign_value(buf.charge, fetch_value(indata, 'charge', float))
    assign_value(buf.fMax, fetch_value(indata, 'fMax', float))
    assign_value(buf.fQuad, fetch_value(indata, 'fQuad', float))
    assign_value(buf.fPSD_t1, fetch_value(indata, 'fPSD_t1', float))
    assign_value(buf.fPSD_t2, fetch_value(indata, 'fPSD_t2', float))
    assign_value(buf.f2inch_maxQ, fetch_value(indata, 'f2inch_maxQ', float))
    assign_value(buf.energy, fetch_value(indata, 'energy', float))
    assign_value(buf.x, fetch_value(indata, 'x', float))
    assign_value(buf.y, fetch_value(indata, 'y', float))
    assign_value(buf.z, fetch_value(indata, 'z', float))
    return buf

def one_iteration(event_number, outdata, fill_buf, helper, callback=lambda e:None):
    timestamp = fill_buf.timestamp[0]
    timestamp_seconds = fill_buf.timestamp_seconds[0]
    timestamp_nanoseconds = fill_buf.timestamp_nanoseconds[0]
    triggerType = fill_buf.triggerType[0]
    detector = fill_buf.detector[0]
    site = fill_buf.site[0]
    nHit = fill_buf.nHit[0]
    charge = fill_buf.charge[0]
    fMax = fill_buf.fMax[0]
    fQuad = fill_buf.fQuad[0]
    fPSD_t1 = fill_buf.fPSD_t1[0]
    fPSD_t2 = fill_buf.fPSD_t2[0]
    f2inch_maxQ = fill_buf.f2inch_maxQ[0]
    energy = fill_buf.energy[0]
    x = fill_buf.x[0]
    y = fill_buf.y[0]
    z = fill_buf.z[0]

    isFlasher = flashers.isFlasher
    isWSMuon = muons.isWSMuon
    isADEvent = adevent.isADEvent

    logging.debug('event cache size: %d', len(helper.event_cache))

    buf = fill_buf.clone_type()

    # assign copied/input values
    assign_value(buf.noTree_loopIndex, event_number)
    assign_value(buf.run, helper.run)
    assign_value(buf.fileno, helper.fileno)
    assign_value(buf.timestamp, timestamp)
    assign_value(buf.timestamp_seconds, timestamp_seconds)
    assign_value(buf.timestamp_nanoseconds, timestamp_nanoseconds)
    assign_value(buf.triggerType, triggerType)
    assign_value(buf.detector, detector)
    assign_value(buf.site, site)
    assign_value(buf.nHit, nHit)
    assign_value(buf.charge, charge)
    assign_value(buf.fMax, fMax)
    assign_value(buf.fQuad, fQuad)
    assign_value(buf.fPSD_t1, fPSD_t1)
    assign_value(buf.fPSD_t2, fPSD_t2)
    assign_value(buf.f2inch_maxQ, f2inch_maxQ)
    assign_value(buf.energy, energy)
    assign_value(buf.x, x)
    assign_value(buf.y, y)
    assign_value(buf.z, z)

    # Initialize dt_next_ADevent to -1
    assign_value(buf.dt_next_ADevent, -1)

    # Compute simple tags and values (those that only require data
    # from the current event)
    event_fID = fID(fMax, fQuad)
    assign_value(buf.fID, event_fID)
    event_fPSD = fPSD(fPSD_t1, fPSD_t2)
    assign_value(buf.fPSD, event_fPSD)
    event_isFlasher = isFlasher(event_fID, event_fPSD, f2inch_maxQ,
            detector)
    assign_value(buf.tag_flasher, event_isFlasher)
    event_isWPevent = isWSMuon(detector, nHit, triggerType)
    assign_value(buf.tag_WPevent, event_isWPevent)
    event_isADEvent = isADEvent(detector, energy)
    assign_value(buf.tag_ADevent, event_isADEvent)

    # Assign the last AD event value since it only depends on this and previous
    # events.
    assign_value(buf.dt_last_ADevent, timestamp -
            helper.last_ADevent_time[detector])


    # Update the dt_previous_* and dt_next_* values

    # This comes after values are assigned to the buffer because we
    # don't want dt_previous_* to be 0.
    if event_isADEvent and event_isFlasher in (0, 2):
        logging.debug("isADEvent")
        helper.last_ADevent_time[detector] = timestamp
        # Decide which events in the event cache to examine, based on whether
        # the event cache has been modified recently
        if helper.go_through_whole_cache:
            # Then the cache has changed substantially since the last
            # iteration, so we will just eat the cost of going through the
            # whole thing.
            pass
        else:
            # Save time by skipping the first chunk of events in the cache that
            # have already been assigned a next_WSMuon. This is achieved by
            # "rotating" the event_cache (a deque type) which means shifting
            # the first n entries to the end, so we start with the new
            # events that need to be labeled. When we reach an already-labeled
            # event, we will be able to exit the loop and un-rotate the cache.
            safe_start_index = helper.first_index_without_next_ADevent
            helper.event_cache.rotate(-1 * safe_start_index)
        # Assign dt_next_WSMuon to the events in the event_cache
        for i, cached_event in enumerate(helper.event_cache):
            # Look for events that have already been assigned a next_WSMuon
            if cached_event.dt_next_ADevent[0] != -1:
                if helper.go_through_whole_cache:
                    # We are going through every event in the cache, even if
                    # it's already been assigned
                    continue
                else:
                    # We have reached the end of the "new" events in the cache
                    # and now we can exit the loop.
                    break
            if cached_event.detector[0] == detector:
                logging.debug('taggedNextADevent%d', cached_event.timestamp[0])
                assign_value(cached_event.dt_next_ADevent,
                        timestamp - cached_event.timestamp[0])
        # Save the last index reached as the new starting location for next
        # time
        if helper.go_through_whole_cache:
            try:
                helper.first_index_without_next_ADevent = i + 1
                helper.go_through_whole_cache = False
            except UnboundLocalError:
                helper.first_index_without_next_ADevent = 0
                helper.go_through_whole_cache = True
        else:
            helper.event_cache.rotate(safe_start_index)
            try:
                helper.first_index_without_next_ADevent = i + safe_start_index
            except UnboundLocalError:
                helper.first_index_without_next_ADevent = 0
                helper.go_through_whole_cache = True

    # Determine which of the oldest events are ready to go into the
    # new TTree. It is possible (due to different ADs) that the
    # oldest event might not be ready but another event might be. In
    # that case, to preserve event order, we will wait for the
    # oldest event to be ready anyways.
    num_to_delete = 0
    all_done_with_cache = True
    cache_size = len(helper.event_cache)
    while all_done_with_cache and num_to_delete < cache_size:
        cached_event = helper.event_cache[num_to_delete]
        if done_with_cache(cached_event):
            logging.debug('event is done with cache')
            num_to_delete += 1
        else:
            logging.debug('event is not done with cache')
            all_done_with_cache = False
    logging.debug('deleting %d events', num_to_delete)
    if num_to_delete > 0:
        helper.go_through_whole_cache = True
    # Remove the oldest events from the cache and fill them into the
    # new TTree, computing the final IBD tags while we're at it.
    for _ in range(num_to_delete):
        cached_event = helper.event_cache.popleft()
        e = cached_event
        coincidence = adevent.hasCoincidence(
                e.tag_ADevent[0],
                e.dt_last_ADevent[0],
                e.dt_next_ADevent[0],
                e.tag_flasher[0])
        high_energy = adevent.hasHighEnergy(e.tag_ADevent[0], e.energy[0])
        logging.debug('cached event: %s', cached_event)
        logging.debug('hasCoincidence: %s', coincidence)
        logging.debug('hasHighEnergy: %s', high_energy)
        logging.debug('isWPevent: %s', cached_event.tag_WPevent[0])
        assign_value(cached_event.tag_coincidence, coincidence)
        callback(cached_event)
        if coincidence or high_energy or cached_event.tag_WPevent[0]:
            logging.debug('saving event')
            cached_event.copyTo(fill_buf)
            outdata.Fill()

    helper.event_cache.append(buf)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--events', type=int, default=-1)
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.events, args.debug)
