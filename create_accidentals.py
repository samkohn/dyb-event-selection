from __future__ import print_function

import argparse
from collections import deque
import math
import os
import random

import numpy as np
np.seterr('raise')

import common
import delayeds
from root_util import assign_value

def distance(a, b):
    '''Distance where coordinates are specified as 'x', 'y', 'z' keys.
    '''
    return math.sqrt(
            (a['x'] - b['x'])**2
            + (a['y'] - b['y'])**2
            + (a['z'] - b['z'])**2)

def get_event(computed):
    event = dict(
            loopIndex=computed.loopIndex,
            timestamp=computed.timestamp,
            triggerNumber=computed.triggerNumber,
            triggerType=computed.triggerType,
            nHit=computed.nHit,
            charge=computed.charge,
            fQuad=computed.fQuad,
            fMax=computed.fMax,
            fPSD_t1=computed.fPSD_t1,
            fPSD_t2=computed.fPSD_t2,
            f2inch_maxQ=computed.f2inch_maxQ,
            energy=computed.energy,
            x=computed.x,
            y=computed.y,
            z=computed.z,
    )
    return event

def sequential_pairing(computed):
    entries = computed.GetEntries()
    # Ignore the last event if there are an odd number of events
    if entries % 2 == 1:
        entries -= 1
    halfway = entries//2
    index = 0
    first_events = [None] * halfway
    while index < entries:
        # Increment until we get a good single
        computed.GetEntry(index)
        event = get_event(computed)
        if index < halfway:
            first_events[halfway - index - 1] = event
        else:
            first = first_events[2 * halfway - index - 1]
            yield first, event
            first_events[2 * halfway - index - 1] = None
        index += 1

def random_pairing_N(computed):
    entries = computed.GetEntries()
    halfway = entries//2
    first_events = [None] * halfway
    second_events = [None] * halfway
    if entries % 2 == 0:
        ordering = list(range(entries))
    else:
        ordering = list(range(entries - 1))
    random.shuffle(ordering)
    for entry, order in enumerate(ordering):
        computed.GetEntry(entry)
        event = get_event(computed)
        if order % 2 == 0:
            first_events[order // 2] = event
            if second_events[order // 2] is not None:
                yield first_events[order // 2], second_events[order // 2]
                first_events[order // 2] = None
                second_events[order // 2] = None
        else:
            second_events[order // 2] = event
            if first_events[order // 2] is not None:
                yield first_events[order // 2], second_events[order // 2]
                first_events[order // 2] = None
                second_events[order // 2] = None

def random_pairing_many(computed, num_samples, cut=lambda entry:False, singles_rate=None,
        seed=1):
    if singles_rate is None:
        singles_rate = 0.01
    entries = computed.GetEntries()
    events = []
    computed.SetBranchStatus('*', 0)
    computed.SetBranchStatus('x', 1)
    computed.SetBranchStatus('y', 1)
    computed.SetBranchStatus('z', 1)
    computed.SetBranchStatus('energy', 1)
    for entry in computed:
        if cut(entry):
            continue
        events.append((entry.x, entry.y, entry.z))
    events = np.array(events)
    rng = np.random.default_rng(seed)
    num_repeats = (num_samples // len(events)) + 1
    event_drs = np.array([])
    for cycle in range(num_repeats):
        first_events = events.copy()
        second_events = events.copy()
        rng.shuffle(first_events, axis=0)
        rng.shuffle(second_events, axis=0)
        event_drs = np.concatenate((event_drs, np.linalg.norm(first_events - second_events, axis=1)))
    event_drs = event_drs[:num_samples]
    dt_value_choices = np.linspace(1e-6, delayeds._NH_THU_MAX_TIME/1e9, 1500)
    # Exponential distribution
    dt_value_probabilities = singles_rate * np.exp(-singles_rate * dt_value_choices)
    dt_value_probabilities /= sum(dt_value_probabilities)
    event_dts_s = rng.choice(dt_value_choices, p=dt_value_probabilities, size=num_samples)
    event_dts = (1e9 * event_dts_s).astype(int)
    event_DTs = delayeds.nH_THU_DT(event_drs, event_dts)
    return event_DTs

def only_DT_eff(infilename, ttree_name, pairing, update_db, **kwargs):
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    computed = infile.Get(ttree_name)
    computed.GetEntry(0)
    run = computed.run
    detector = computed.detector
    site = computed.site
    num_pairs = kwargs['num_pairs']
    seed = kwargs['seed']
    if computed.GetEntries() < 61200:
        efficiency = None
        error = None
    else:
        if pairing == 'random_many':
            event_DTs = random_pairing_many(computed, num_pairs, seed=seed)
            DT_CUT = delayeds._NH_THU_DIST_TIME_MAX
            num_passes_cut = np.count_nonzero(event_DTs < DT_CUT)
            efficiency = num_passes_cut / num_pairs
            error = math.sqrt(num_pairs * efficiency * (1 - efficiency)) / num_pairs
        elif pairing == 'random_many_resid_flasher':
            def cut(entry):
                x, y, z = entry.x, entry.y, entry.z
                return z > 2200 and np.sqrt(x*x + y*y) > 500
            event_DTs = random_pairing_many(computed, num_pairs, cut, seed=seed)
            DT_CUT = delayeds._NH_THU_DIST_TIME_MAX
            num_passes_cut = np.count_nonzero(event_DTs < DT_CUT)
            efficiency = num_passes_cut / num_pairs
            error = math.sqrt(num_pairs * efficiency * (1 - efficiency)) / num_pairs
        elif pairing == 'energy_lt_2MeV':
            def cut(entry):
                return entry.energy < 2
            event_DTs = random_pairing_many(computed, num_pairs, cut, seed=seed)
            DT_CUT = delayeds._NH_THU_DIST_TIME_MAX
            num_passes_cut = np.count_nonzero(event_DTs < DT_CUT)
            efficiency = num_passes_cut / num_pairs
            error = math.sqrt(num_pairs * efficiency * (1 - efficiency)) / num_pairs
        elif pairing == 'energy_gt_2MeV':
            def cut(entry):
                return entry.energy >= 2
            event_DTs = random_pairing_many(computed, num_pairs, cut, seed=seed)
            DT_CUT = delayeds._NH_THU_DIST_TIME_MAX
            num_passes_cut = np.count_nonzero(event_DTs < DT_CUT)
            efficiency = num_passes_cut / num_pairs
            error = math.sqrt(num_pairs * efficiency * (1 - efficiency)) / num_pairs
        else:
            raise NotImplemented(pairing)
        try:
            percent_error = 100 * error / efficiency
        except ZeroDivisionError:
            percent_error = 0
    pairing = pairing + '_expo_time'
    if update_db is None:
        print(f'Pairing type: {pairing}')
        print(f'Efficiency: {efficiency:.6f} +/- {error:.6f} ({percent_error:.1f}%)')
        print(f'Total pairs: {num_pairs}')
        print(f'Passed DT cut: {num_passes_cut}')
    else:
        with common.get_db(update_db) as conn:
            cursor = conn.cursor()
            cursor.execute('''INSERT OR REPLACE INTO distance_time_eff_study
                VALUES (?, ?, ?, ?, ?, ?)''',
                (run, detector, pairing, efficiency, error, num_pairs))
    return


def is_complete(infilename, outfilename):
    """Check to see if the output file has at least 48% as many entries
    as the input.
    """
    if not os.path.isfile(outfilename):
        return False
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    in_singles = infile.Get('singles')
    in_entries = in_singles.GetEntries()
    infile.Close()
    outfile = ROOT.TFile(outfilename, 'READ')
    out_acc = outfile.Get('all_pairs')
    if not out_acc:  # PyROOT-speak for null pointer test
        return False
    out_entries = out_acc.GetEntries()
    PAIRING_FACTOR = 0.48
    threshold = in_entries * PAIRING_FACTOR
    if out_entries < threshold:
        return False
    return True


def main(infilename, outfile, ttree_name, pairing_algorithm, pairing_note, update_db):
    import process
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    outfile = ROOT.TFile(outfile, 'RECREATE')
    acc, acc_buf = process.create_computed_TTree('accidentals', outfile,
            'nh_THU', 'Accidentals sample (git: %s)')
    all_acc, all_acc_buf = process.create_computed_TTree('all_pairs',
            outfile, 'nh_THU', 'All paired singles (git: %s)')
    outfile.cd()
    acc_spectrum_hist = ROOT.TH2F('acc_spectrum', 'acc_spectrum',
            210, 1.5, 12, 210, 1.5, 12)
    eps_DT_hist = ROOT.TH2F('DT_cut_fails', 'DT_cut_fails',
            210, 1.5, 12, 210, 1.5, 12)
    computed = infile.Get(ttree_name)
    computed.SetBranchStatus('*', 0)
    computed.SetBranchStatus('run', 1)
    computed.SetBranchStatus('detector', 1)
    computed.SetBranchStatus('site', 1)
    computed.SetBranchStatus('energy', 1)
    computed.SetBranchStatus('x', 1)
    computed.SetBranchStatus('y', 1)
    computed.SetBranchStatus('z', 1)
    computed.GetEntry(0)
    run = computed.run
    detector = computed.detector
    site = computed.site
    computed.SetBranchStatus('run', 0)
    computed.SetBranchStatus('detector', 0)
    computed.SetBranchStatus('site', 0)
    if pairing_algorithm == 'sequential':
        #first_events, second_events = sequential_pairing(computed)
        generator = sequential_pairing(computed)
    elif pairing_algorithm == 'random_N':
        #first_events, second_events = random_pairing_N(computed)
        generator = random_pairing_N(computed)
    else:
        raise NotImplemented(pairing_algorithm)
    #for (first_event, second_event) in zip(first_events, second_events):
    event_DTs = []
    for (first_event, second_event) in generator:
        # Compare event distances and fill the tree
        event_dr = distance(first_event, second_event)
        event_dt = random.randint(1000, delayeds._NH_THU_MAX_TIME)
        all_acc_buf.multiplicity[0] = 2
        assign_value(all_acc_buf.run, run)
        assign_value(all_acc_buf.site, site)
        assign_value(all_acc_buf.loopIndex, first_event['loopIndex'], 0)
        assign_value(all_acc_buf.timestamp, first_event['timestamp'], 0)
        assign_value(all_acc_buf.timestamp_seconds,
                first_event['timestamp'] // 1000000000, 0)
        assign_value(all_acc_buf.timestamp_nanoseconds,
                first_event['timestamp'] % 1000000000, 0)
        assign_value(all_acc_buf.detector, detector, 0)
        assign_value(all_acc_buf.triggerNumber, first_event['triggerNumber'], 0)
        assign_value(all_acc_buf.triggerType, first_event['triggerType'], 0)
        assign_value(all_acc_buf.nHit, first_event['nHit'], 0)
        assign_value(all_acc_buf.charge, first_event['charge'], 0)
        assign_value(all_acc_buf.fQuad, first_event['fQuad'], 0)
        assign_value(all_acc_buf.fMax, first_event['fMax'], 0)
        assign_value(all_acc_buf.fPSD_t1, first_event['fPSD_t1'], 0)
        assign_value(all_acc_buf.fPSD_t2, first_event['fPSD_t2'], 0)
        assign_value(all_acc_buf.f2inch_maxQ, first_event['f2inch_maxQ'], 0)
        assign_value(all_acc_buf.energy, first_event['energy'], 0)
        assign_value(all_acc_buf.x, first_event['x'], 0)
        assign_value(all_acc_buf.y, first_event['y'], 0)
        assign_value(all_acc_buf.z, first_event['z'], 0)
        assign_value(all_acc_buf.loopIndex, second_event['loopIndex'], 1)
        assign_value(all_acc_buf.timestamp, second_event['timestamp'], 1)
        assign_value(all_acc_buf.timestamp_seconds,
                second_event['timestamp'] // 1000000000, 1)
        assign_value(all_acc_buf.timestamp_nanoseconds,
                second_event['timestamp'] % 1000000000, 1)
        assign_value(all_acc_buf.detector, detector, 0)
        assign_value(all_acc_buf.triggerNumber, second_event['triggerNumber'], 1)
        assign_value(all_acc_buf.triggerType, second_event['triggerType'], 1)
        assign_value(all_acc_buf.nHit, second_event['nHit'], 1)
        assign_value(all_acc_buf.charge, second_event['charge'], 1)
        assign_value(all_acc_buf.fQuad, second_event['fQuad'], 1)
        assign_value(all_acc_buf.fMax, second_event['fMax'], 1)
        assign_value(all_acc_buf.fPSD_t1, second_event['fPSD_t1'], 1)
        assign_value(all_acc_buf.fPSD_t2, second_event['fPSD_t2'], 1)
        assign_value(all_acc_buf.f2inch_maxQ, second_event['f2inch_maxQ'], 1)
        assign_value(all_acc_buf.energy, second_event['energy'], 1)
        assign_value(all_acc_buf.x, second_event['x'], 1)
        assign_value(all_acc_buf.y, second_event['y'], 1)
        assign_value(all_acc_buf.z, second_event['z'], 1)
        all_acc_buf.dr_to_prompt[0] = 0
        all_acc_buf.dr_to_prompt[1] = event_dr
        all_acc_buf.dt_to_prompt[0] = 0
        all_acc_buf.dt_to_prompt[1] = event_dt
        event_DT = delayeds.nH_THU_DT(event_dr, event_dt)
        event_DTs.append(event_DT)
        if event_DT < delayeds._NH_THU_DIST_TIME_MAX:
            acc_spectrum_hist.Fill(first_event['energy'],
                    second_event['energy'])
            acc_spectrum_hist.Fill(second_event['energy'],
                    first_event['energy'])
            all_acc_buf.copyTo(acc_buf)
            acc.Fill()
        else:
            eps_DT_hist.Fill(first_event['energy'],
                    second_event['energy'])
            eps_DT_hist.Fill(second_event['energy'],
                    first_event['energy'])
        all_acc.Fill()
    if update_db is not None:
        event_DTs = np.array(event_DTs)
        num_pairs = len(event_DTs)
        DT_CUT = delayeds._NH_THU_DIST_TIME_MAX
        num_passes_cut = np.count_nonzero(event_DTs < DT_CUT)
        efficiency = num_passes_cut / num_pairs
        error = math.sqrt(num_pairs * efficiency * (1 - efficiency)) / num_pairs
        with common.get_db(update_db) as conn:
            cursor = conn.cursor()
            cursor.execute('''INSERT OR REPLACE INTO distance_time_eff_study
                VALUES (?, ?, ?, ?, ?, ?)''',
                (run, detector, f'{pairing_algorithm}; {pairing_note}', efficiency, error, num_pairs))
    outfile.Write()
    outfile.Close()
    infile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('--ttree-name', default='singles')
    parser.add_argument('--pairing', default='sequential')
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--only-DT-eff', action='store_true')
    parser.add_argument('--update-db')
    parser.add_argument('--num-pairs', type=int)
    parser.add_argument('--pairing-note', default='')
    args = parser.parse_args()

    random.seed(args.seed)
    if args.only_DT_eff:
        only_DT_eff(args.infile, args.ttree_name, args.pairing,
                args.update_db, num_pairs=args.num_pairs, seed=args.seed)
    else:
        main(args.infile, args.outfile, args.ttree_name, args.pairing, args.pairing_note,
                args.update_db)
