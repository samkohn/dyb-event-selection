from __future__ import print_function

import argparse
import math

def distance(a, b):
    '''Distance where coordinates are specified as 'x', 'y', 'z' keys.
    '''
    return math.sqrt(
            (a['x'] - b['x'])**2
            + (a['y'] - b['y'])**2
            + (a['z'] - b['z'])**2)

def is_single(ttree_event):
    'Applies the single event criteria to the loaded entry in this TTree.'
    return (
            ttree_event.energy < 12 and ttree_event.energy > 1.5
            and not (ttree_event.tag_WSMuonVeto
                or ttree_event.tag_ADMuonVeto
                or ttree_event.tag_ShowerMuonVeto)
            and ttree_event.coincidence_number == 1
            and ttree_event.num_ADevents_800us_all == 1
            and ttree_event.dts_ADevents_800us_all[0] > 400e3)

def main(infilename, outfile, AD, ttree_name):
    import process
    import ROOT
    infile = ROOT.TFile(infilename, 'READ')
    infile2 = ROOT.TFile(infilename, 'READ')
    outfile = ROOT.TFile(outfile, 'RECREATE')
    acc, acc_buf = process.create_computed_TTree('accidentals', outfile,
            'nh_THU', 'Accidentals sample (git: %s)')
    outfile.cd()
    acc_spectrum_hist = ROOT.TH2F('acc_spectrum', 'acc_spectrum',
            210, 1.5, 12, 210, 1.5, 12)
    eps_distance_hist = ROOT.TH2F('distance_cut_fails', 'distance_cut_fails',
            210, 1.5, 12, 210, 1.5, 12)
    computed = infile.Get(ttree_name)
    computed.SetBranchStatus('*', 0)
    computed.SetBranchStatus('detector', 1)
    computed.SetBranchStatus('energy', 1)
    computed.SetBranchStatus('tag_WSMuonVeto', 1)
    computed.SetBranchStatus('tag_ADMuonVeto', 1)
    computed.SetBranchStatus('tag_ShowerMuonVeto', 1)
    computed.SetBranchStatus('coincidence_number', 1)
    computed.SetBranchStatus('num_ADevents_800us_all', 1)
    computed.SetBranchStatus('dts_ADevents_800us_all', 1)
    computed.SetBranchStatus('x', 1)
    computed.SetBranchStatus('y', 1)
    computed.SetBranchStatus('z', 1)
    computed2 = infile2.Get(ttree_name)
    computed2.SetBranchStatus('*', 0)
    computed2.SetBranchStatus('detector', 1)
    computed2.SetBranchStatus('energy', 1)
    computed2.SetBranchStatus('tag_WSMuonVeto', 1)
    computed2.SetBranchStatus('tag_ADMuonVeto', 1)
    computed2.SetBranchStatus('tag_ShowerMuonVeto', 1)
    computed2.SetBranchStatus('coincidence_number', 1)
    computed2.SetBranchStatus('num_ADevents_800us_all', 1)
    computed2.SetBranchStatus('dts_ADevents_800us_all', 1)
    computed2.SetBranchStatus('x', 1)
    computed2.SetBranchStatus('y', 1)
    computed2.SetBranchStatus('z', 1)
    entries = computed.GetEntries()
    halfway = entries//2
    first_half_index = 0
    second_half_index = halfway
    phase = 1
    first_event = {}
    second_event = {}
    while first_half_index < halfway and second_half_index < entries:
        if phase == 1:
            # Increment first_half_index until we get a good single
            computed.GetEntry(first_half_index)
            if (computed.detector == AD
                    and is_single(computed)):
                first_event['energy'] = computed.energy
                first_event['detector'] = computed.detector
                first_event['x'] = computed.x
                first_event['y'] = computed.y
                first_event['z'] = computed.z
                phase = 2
            else:
                first_half_index += 1
        if phase == 2:
            # Increment second_half_index until we get a good single
            computed2.GetEntry(second_half_index)
            if (computed2.detector == AD
                    and is_single(computed2)):
                second_event['energy'] = computed2.energy
                second_event['detector'] = computed2.detector
                second_event['x'] = computed2.x
                second_event['y'] = computed2.y
                second_event['z'] = computed2.z
                phase = 3
            else:
                second_half_index += 1
        if phase == 3:
            # Compare event distances and fill the tree
            event_dr = distance(first_event, second_event)
            if event_dr < 500:
                acc_spectrum_hist.Fill(first_event['energy'],
                        second_event['energy'])
                acc_spectrum_hist.Fill(second_event['energy'],
                        first_event['energy'])
                acc_buf.energy[0] = first_event['energy']
                acc_buf.energy_previous_PromptLike[0] = second_event['energy']
                acc_buf.detector[0] = first_event['detector']
                acc.Fill()
                print('Found an accidental')
            else:
                eps_distance_hist.Fill(first_event['energy'],
                        second_event['energy'])
                eps_distance_hist.Fill(second_event['energy'],
                        first_event['energy'])
            first_half_index += 1
            second_half_index += 1
            phase = 1
    outfile.Write()
    outfile.Close()
    infile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('ad', type=int)
    parser.add_argument('--ttree-name', default='computed')
    args = parser.parse_args()

    main(args.infile, args.outfile, args.ad, args.ttree_name)
