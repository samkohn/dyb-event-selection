import argparse
import logging
from collections import deque
import os
from functools import lru_cache

from common import *
import process
import translate
import flashers
import first_pass
from translate import initialize_indata_onefile as initialize
from root_util import *

class RawFileAdapter():
    ATTR_LOOKUP = {
            'Q1': ('Q1', float),
            'Q2': ('Q2', float),
    }
    def __init__(self, ttree, run, fileno, site):
        from ROOT import TFile
        self.ttree = ttree
        self.run = run
        self.fileno = fileno
        self.site = site
        self.ttree.SetBranchStatus('execNumber', 1)
        for branch_name, _ in self.ATTR_LOOKUP.values():
            self.ttree.SetBranchStatus(branch_name, 1)

    def GetEntry(self, index):
        self.ttree.GetEntry(index)
        self.__getattr__.cache_clear()

    def GetEntries(self):
        return self.ttree.GetEntries()

    @lru_cache(maxsize=32)
    def __getattr__(self, name):
        attr_name = self.ATTR_LOOKUP.get(name, None)
        if attr_name is None:
            raise AttributeError('No attribute "{}"'.format(name))
        return fetch_value(self.ttree, attr_name[0], attr_name[1])

    def SetBranchStatus(self, *args):
        self.ttree.SetBranchStatus(*args)


class EventsCombiner:
    """Load up a set of events files and interleave them in event order.

    >>> # input arg is a dict of ad # to TTree object
    >>> events = EventsCombiner({1: t1, 2: t2})
    >>> # Iterate
    >>> for ad, loopIndex in events:
    ...    nuwa_ttree.GetEntry(loopIndex)
    ...    # Do things with nuwa_ttree
    """
    def __init__(self, ttrees):
        self.ttrees = ttrees  # {ad: ttree}
        self.ads = list(ttrees.keys())
        for ttree in ttrees.values():
            ttree.GetEntry(0)
        self.last_loopindexes = {ad: -1 for ad in self.ads}
        self.entries = {ad: 0 for ad in self.ads}
        self.last_lowest_ad = None
        self.finished = {ad: False for ad in self.ads}

    def __iter__(self):
        return self

    def __next__(self):
        if all(finished for finished in self.finished.values()):
            raise StopIteration
        if len(self.ttrees) == 1:
            ad = self.ads[0]
            status = self.ttrees[ad].GetEntry(self.entries[ad] + 1)
            if status <= 0:
                raise StopIteration
            return (ad, self.ttrees[ad].loopIndex)
        new_lowest_loopindex = 1e50
        new_lowest_ad = None
        for ad, loopIndex in self.last_loopindexes.items():
            if self.finished[ad]:
                continue
            if loopIndex == -1:
                # must be at the beginning
                self.ttrees[ad].GetEntry(0)
                if self.ttrees[ad].loopIndex < new_lowest_loopindex:
                    new_lowest_loopindex = self.ttrees[ad].loopIndex
                    new_lowest_ad = ad
            elif ad == self.last_lowest_ad:
                # increment the entry in this AD
                new_entry = self.entries[ad] + 1
                status = self.ttrees[ad].GetEntry(new_entry)
                if status <= 0:
                    self.finished[ad] = True
                    if all(finished for finished in self.finished.values()):
                        raise StopIteration
                    continue  # already overran the length of the TTree
                self.entries[ad] = new_entry
                if self.ttrees[ad].loopIndex < new_lowest_loopindex:
                    new_lowest_loopindex = self.ttrees[ad].loopIndex
                    new_lowest_ad = ad
            else:
                # this AD wasn't the lowest last time but might be this time
                if self.ttrees[ad].loopIndex < new_lowest_loopindex:
                    new_lowest_loopindex = self.ttrees[ad].loopIndex
                    new_lowest_ad = ad
        self.last_loopindexes[new_lowest_ad] = new_lowest_loopindex
        self.last_lowest_ad = new_lowest_ad
        return (new_lowest_ad, new_lowest_loopindex)


def main_loop(events, indata, flasher_ttrees, events_ttrees, ads, debug):
    merged_events = EventsCombiner(events_ttrees)
    for ad, loopIndex in merged_events:
        if loopIndex > events:
            break
        indata.GetEntry(loopIndex)
        ttree, fill_buf = flasher_ttrees[ad]
        load_flasher_buf(fill_buf, indata, loopIndex)
        ttree.Fill()
    return


def load_flasher_buf(buf, indata, loopIndex):
    assign_value(buf.loopIndex, loopIndex)
    assign_value(buf.Q1, indata.Q1)
    assign_value(buf.Q2, indata.Q2)
    return


def create_flashers_TTree(host_file):
    from ROOT import TTree
    host_file.cd()
    title = 'flasher values by Sam Kohn (git: {})'.format(
            translate.git_describe())
    out = TTree('flashers', title)
    buf = TreeBuffer()
    buf.loopIndex = unsigned_int_value()
    buf.Q1 = float_value()
    buf.Q2 = float_value()

    def branch(name, typecode):
        out.Branch(name, getattr(buf, name), '{}/{}'.format(name,
            typecode))
        return

    branch('loopIndex', 'i')
    branch('Q1', 'F')
    branch('Q2', 'F')
    return out, buf

def create_outfiles(out_location, run, fileno, ads):
    from ROOT import TFile
    flashers_name = 'q1q2_ad{}_{}_{:>04}.root'.format('{}', run, fileno)
    flashersFiles = {ad: TFile(os.path.join(out_location, flashers_name.format(ad)),
        'RECREATE') for ad in ads}
    return flashersFiles

def main(num_events, infile, out_location, runfilesite, events_filenames, debug):
    from ROOT import TFile
    run, fileno, site = runfilesite
    infile = TFile(infile, 'READ')
    calibStats, adSimple = initialize(infile, 'AdSimpleNL')
    indata = RawFileAdapter(calibStats, run, fileno, site)
    ads = dets_for(indata.site, run)
    outfiles = create_outfiles(out_location, run, fileno, ads)
    flasher_ttrees = {ad: create_flashers_TTree(f) for ad, f in outfiles.items()}

    # Load lists of AD events so we can skip WS, RPC, etc. events
    events_files = [TFile(filename, 'READ') for filename in events_filenames]
    events_ttrees = {ad: f.Get('events') for ad, f in zip(ads, events_files)}

    if num_events == -1:
        num_events = indata.GetEntries()
    main_loop(num_events, indata, flasher_ttrees, events_ttrees, ads, debug)
    for f in outfiles.values():
        f.Write()
        f.Close()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file name')
    parser.add_argument('-o', '--output', help='Output location')
    parser.add_argument(
        '-r', '--runfilesite', nargs=3, type=int, help='<run> <fileno> <site>'
    )
    parser.add_argument(
        '--events-files', nargs='+', help='first_pass outputs for this run'
    )
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('-n', '--num-events', type=int, default=-1)
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    main(args.num_events, args.input, args.output, args.runfilesite, args.events_files, args.debug)
