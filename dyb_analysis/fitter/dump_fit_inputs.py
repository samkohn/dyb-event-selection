"""Create a portable representation of fitter inputs for sharing."""

import argparse

import ROOT

import prediction as pred

def spec_to_ROOT_hist(values, bins, name):
    num_bins = len(values)
    low_bound = bins[0]
    up_bound = bins[-1]
    new_hist = ROOT.TH1D(name, name, num_bins, low_bound, up_bound)
    new_hist.SetBins(num_bins, bins)
    for i, value in enumerate(values):
        new_hist.SetBinContent(i + 1, value)
    return new_hist

def constants_to_dump_str(constants):
    """Convert the FitConstants object into a string suitable for dumping."""
    dump_notes = """# Fitter inputs
# Skip lines starting with #
# Each comment/label applies to the *following* line
# Each data line has comma-separated values in AD order:
# EH1-AD1,EH1-AD2,EH2-AD1,EH2-AD2,EH3-AD1,EH3-AD2,EH3-AD3,EH3-AD4
"""
    num_ads = 8
    names = []
    rows = []
    names.append('Total livetime in days')
    rows.append([time_s/60/60/24 for time_s in constants.daq_livetimes.values()])
    names.append('Muon veto efficiency')
    rows.append(list(constants.muon_eff.values()))
    names.append('Multiplicity efficiency')
    rows.append(list(constants.multiplicity_eff.values()))
    names.append('Relative uncertainty in detection efficiency')
    rows.append([constants.efficiency_err] * num_ads)
    names.append('Relative uncertainty in reactor power')
    rows.append([constants.reactor_err] * num_ads)
    names.append('Relative uncertainty in relative energy scale')
    rows.append([constants.rel_escale_err] * num_ads)
    names.append('Total target protons relative to EH1-AD1')
    rows.append([mass/constants.masses[1, 1] for mass in constants.masses.values()])
    final_string = dump_notes
    for name, row in zip(names, rows):
        final_string += f'# {name}\n'
        final_string += ','.join(map(str,row))
        final_string += '\n'
    return final_string


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fit_config')
    parser.add_argument('outfile_root')
    parser.add_argument('outfile_txt')
    args = parser.parse_args()
    outfile_root = ROOT.TFile(args.outfile_root, 'RECREATE')
    constants = pred.load_constants(args.fit_config)
    observed_hists = {}
    accidentals_hists = {}
    bins = constants.reco_bins
    for hall, det in pred.all_ads:
        num_coincidences = constants.observed_candidates[hall, det]
        name = f'observed_counts_EH{hall}_AD{det}'
        observed_hists[hall, det] = spec_to_ROOT_hist(num_coincidences, bins, name)
        num_accidentals = constants.nominal_bgs['accidental'][hall, det]
        name = f'accidentals_EH{hall}_AD{det}'
        accidentals_hists[hall, det] = spec_to_ROOT_hist(num_accidentals, bins, name)
    outfile_root.Write()
    outfile_root.Close()
    txt_dump = constants_to_dump_str(constants)
    with open(args.outfile_txt, 'w') as outfile_txt:
        outfile_txt.write(txt_dump)


