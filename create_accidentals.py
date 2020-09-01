from __future__ import print_function

import argparse
import math
import random

import delayeds

def distance(a, b):
    '''Distance where coordinates are specified as 'x', 'y', 'z' keys.
    '''
    return math.sqrt(
            (a['x'] - b['x'])**2
            + (a['y'] - b['y'])**2
            + (a['z'] - b['z'])**2)

def get_event(computed):
    event = {
            'energy': computed.energy,
            'detector': computed.detector,
            'x': computed.x,
            'y': computed.y,
            'z': computed.z,
            }
    return event

def sequential_pairing(computed):
    entries = computed.GetEntries()
    halfway = entries//2
    index = 0
    first_events = []
    second_events = []
    while index < entries:
        # Increment until we get a good single
        computed.GetEntry(index)
        event = get_event(computed)
        if index < halfway:
            first_events.append(event)
        else:
            second_events.append(event)
        index += 1
    return first_events, second_events

def random_pairing_N(computed):
    entries = computed.GetEntries()
    halfway = entries//2
    first_events = [None] * halfway
    second_events = [None] * halfway
    if entries % 2 == 0:
        ordering = list(range(entries))
    else:
        ordering = list(range(entries - 1))
    for entry, order in enumerate(ordering):
        computed.GetEntry(entry)
        event = get_event(computed)
        if order % 2 == 0:
            first_events[order // 2] = event
        else:
            second_events[order // 2] = event
    return first_events, second_events

def random_pairing_N_repeats(computed):
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
        else:
            second_events[order // 2] = event
    return first_events, second_events

def main(infilename, outfile, ttree_name, pairing_algorithm, repeat):
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
    computed.SetBranchStatus('energy', 1)
    computed.SetBranchStatus('x', 1)
    computed.SetBranchStatus('y', 1)
    computed.SetBranchStatus('z', 1)
    if pairing_algorithm == 'sequential':
        first_events, second_events = sequential_pairing(computed)
    elif pairing_algorithm == 'random_N':
        if repeat > 1:
            efficiency_list = []
            for _ in range(repeat):
                first_events, second_events = random_pairing_N(computed)
                num_passes_DT = 0
                for (first_event, second_event) in zip(first_events, second_events):
                    # Compare event distances and fill the tree
                    event_dr = distance(first_event, second_event)
                    event_dt = random.randint(1000, delayeds._NH_THU_MAX_TIME)
                    if delayeds.nH_THU_DT(event_dr, event_dt) < delayeds._NH_THU_DIST_TIME_MAX:
                        num_passes_DT += 1
                print(num_passes_DT)
                efficiency = num_passes_DT / len(first_events)
                print(efficiency)
                efficiency_list.append(efficiency)
            print(efficiency_list)
            return
        else:
            first_events, second_events = random_pairing_N(computed)
    else:
        raise NotImplemented(pairing_algorithm)
    for (first_event, second_event) in zip(first_events, second_events):
        # Compare event distances and fill the tree
        event_dr = distance(first_event, second_event)
        event_dt = random.randint(1000, delayeds._NH_THU_MAX_TIME)
        all_acc_buf.multiplicity[0] = 2
        all_acc_buf.energy[0] = first_event['energy']
        all_acc_buf.energy[1] = second_event['energy']
        all_acc_buf.detector[0] = first_event['detector']
        all_acc_buf.detector[1] = second_event['detector']
        all_acc_buf.dr_to_prompt[0] = 0
        all_acc_buf.dr_to_prompt[1] = event_dr
        all_acc_buf.dt_to_prompt[0] = 0
        all_acc_buf.dt_to_prompt[1] = event_dt
        if delayeds.nH_THU_DT(event_dr, event_dt) < delayeds._NH_THU_DIST_TIME_MAX:
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
    outfile.Write()
    outfile.Close()
    infile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('--ttree-name', default='singles')
    parser.add_argument('--pairing', default='sequential')
    parser.add_argument('--repeat', type=int, default=-1)
    args = parser.parse_args()

    main(args.infile, args.outfile, args.ttree_name, args.pairing, args.repeat)
