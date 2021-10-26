'''
Take output JSON files and create the rate-only analysis.

'''
from __future__ import print_function
import json
import math

from rate_only import sin22theta13, r_measured

NS_PER_DAY = 1e9*60*60*24

# load up JSON files
filelist = []
with open('outfile_list.txt', 'r') as f:
    for line in f:
        filelist.append(line.strip())

DETS = (1, 2, 3, 4, 5, 6, 7, 8)
ibds = {det: 0 for det in DETS}
tdaq = {det: 0 for det in DETS}
tmuveto = {det: 0 for det in DETS}
multeff_weighted = {det: 0 for det in DETS}
bkg_weighted = {det: 0 for det in DETS}


for outfile in filelist:
    with open(outfile, 'r') as f:
        result = json.load(f)
    site = result['site']
    if site == 1:
        dets = {1: 1, 2: 2}
    elif site == 2:
        dets = {1: 3, 2: 8}
    elif site == 3:
        dets = {1: 4, 2: 5, 3: 6, 4: 7}
    else:
        raise ValueError('Invalid site %s' % site)

    for ehdet, addet in dets.items():
        ibds[addet] += result['number_IBDs'][str(ehdet)]
        tdaq[addet] += result['daq_livetime']
        tmuveto[addet] += result['usable_livetime'][str(ehdet)]
        multeff_weighted[addet] += (
                result['multiplicity_efficiency'][str(ehdet)]
                * result['daq_livetime'])
        bkg_weighted[addet] += (
                result['accidental_rate_perday'][str(ehdet)]
                * result['usable_livetime'][str(ehdet)]
                * result['multiplicity_efficiency'][str(ehdet)])

ibdrate = {det: 0 for det in DETS}
ibdrate_error = {det: 0 for det in DETS}
mueff = {det: 0 for det in DETS}
multeff = {det: 0 for det in DETS}
bkg = {det: 0 for det in DETS}
for det in DETS:
    multeff[det] = multeff_weighted[det]/tdaq[det]
    mueff[det] = float(tmuveto[det])/tdaq[det]
    bkg[det] = bkg_weighted[det]/(tdaq[det] * mueff[det] * multeff[det])
    ibdrate[det] = ibds[det]/(tdaq[det]/NS_PER_DAY * mueff[det] * multeff[det]) - bkg[det]
    ibdrate_error[det] = math.sqrt(ibds[det])/(tdaq[det]/NS_PER_DAY * mueff[det] * multeff[det])
    print('AD%d: livetime %.1f | mueff %.1f | multeff %.1f' % (det,
        tdaq[det]/NS_PER_DAY, 100*mueff[det],
        100*multeff[det]))

ibd_rates_final = {
        'EH1': {1: ibdrate[1], 2: ibdrate[2]},
        'EH2': {1: ibdrate[3], 2: ibdrate[8]},
        'EH3': {1: ibdrate[4], 2: ibdrate[5], 3: ibdrate[6], 4: ibdrate[7]}
        }
livetimes_final = {
        'EH1': {1: tdaq[1], 2: tdaq[2]},
        'EH2': {1: tdaq[3], 2: tdaq[8]},
        'EH3': {1: tdaq[4], 2: tdaq[5], 3: tdaq[6], 4: tdaq[7]}
        }
ibd_rate_error_final = {
        'EH1': {1: ibdrate_error[1], 2: ibdrate_error[2]},
        'EH2': {1: ibdrate_error[3], 2: ibdrate_error[8]},
        'EH3': {1: ibdrate_error[4], 2: ibdrate_error[5], 3: ibdrate_error[6], 4: ibdrate_error[7]}
        }

for det in DETS:
    print('AD%d: ' % det, end='')
    print('%.2f +/- %.2f' % (ibdrate[det], ibdrate_error[det]))
print('R = %.3f +/- %.3f' % r_measured(ibd_rates_final, livetimes_final,
    ibd_rate_error_final))
theta13, theta13_error = sin22theta13(ibd_rates_final, livetimes_final, ibd_rate_error_final)
print('sin^2(2theta_13) = %.3f +/- %.3f' % (theta13, theta13_error))
