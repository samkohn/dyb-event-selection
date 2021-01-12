"""Loop through all "hadded" data files and save the coincidences based on
which reconstruction's distance-time cut is satisfied."""

import argparse
import itertools as it
import multiprocessing

import common
import delayeds

ADSIMPLE_DT_CUT = delayeds._NH_THU_DIST_TIME_CUT_STR
ADTIME_DT_CUT = (
    f'(dr_to_prompt_AdTime[1] + ({delayeds._NH_THU_DIST_TIME_CONST} '
    f'* dt_to_prompt[1]) < {delayeds._NH_THU_DIST_TIME_MAX})'
)

def one_file(run_key, data_file_path, energy_lookup):
    import ROOT
    run, site, ad = run_key
    if run % 10 == 0:
        print(run_key)
    data_file_name = data_file_path.format(run=run, site=site, ad=ad)
    data_file = ROOT.TFile(data_file_name, 'READ')
    ad_events = data_file.Get('ad_events')
    delayed_min_adsimple, delayed_max_adsimple = energy_lookup['nominal', site, ad]
    delayed_min_adtime, delayed_max_adtime = energy_lookup['adtime', site, ad]
    #HIST_NBINS = 105
    #hists_for_ad = (
        #ROOT.TH1F(f'both_eh{site}ad{ad}', f'both_eh{site}ad{ad}', HIST_NBINS, 1.5, 12),
        #ROOT.TH1F(f'adsimple_only_eh{site}ad{ad}', f'adsimple_only_eh{site}ad{ad}',
            #HIST_NBINS, 1.5, 12)
        #ROOT.TH1F(f'adtime_only_eh{site}ad{ad}', f'adtime_only_eh{site}ad{ad}',
            #HIST_NBINS, 1.5, 12)
    #)
    PASSES_ADSIMPLE = (
        f'(energy[1] > {delayed_min_adsimple} && '
        f'energy[1] < {delayed_max_adtime} && '
        f'{ADSIMPLE_DT_CUT})'
    )
    PASSES_ADTIME = (
        f'(energy[1] > {delayed_min_adtime} && '
        f'energy[1] < {delayed_max_adtime} && '
        f'{ADTIME_DT_CUT})'
    )
    num_passes_both = ad_events.Draw(
        f'energy[0] >> both_eh{site}ad{ad}',
        f'{PASSES_ADSIMPLE} && {PASSES_ADTIME}',
        'goff'
    )
    #print('num_passes_both:', num_passes_both)
    num_passes_only_adsimple = ad_events.Draw(
        f'energy[0] >> adsimple_only_eh{site}ad{ad}',
        f'{PASSES_ADSIMPLE} && !{PASSES_ADTIME}',
        'goff'
    )
    #print('num_passes_only_adsimple:', num_passes_only_adsimple)
    num_passes_only_adtime = ad_events.Draw(
        f'energy[0] >> adtime_only_eh{site}ad{ad}',
        f'!{PASSES_ADSIMPLE} && {PASSES_ADTIME}',
        'goff'
    )
    #print('num_passes_only_adtime:', num_passes_only_adtime)
    num_passes_neither = ad_events.Draw('energy[0]',
        f'!{PASSES_ADSIMPLE} && !{PASSES_ADTIME}',
        'goff'
    )
    #print('num_passes_neither:', num_passes_neither)
    # Construct return arrays
    #both_histo = []
    #adsimple_histo = []
    #adtime_histo = []
    #both, adsimple, adtime = hists_for_ad
    #for i in range(HIST_NBINS):
        #both_histo.append(both.GetBinContent(i + 1))
        #adsimple_histo.append(adsimple.GetBinContent(i + 1))
        #adtime_histo.append(adtime.GetBinContent(i + 1))
    return (
        run,
        ad,
        num_passes_both,
        num_passes_only_adsimple,
        num_passes_only_adtime,
        num_passes_neither,
        #both_histo,
        #adsimple_histo,
        #adtime_histo,
    )

def main(database, data_file_path, update_db):
    import ROOT
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        # Fetch all triplets of RunNo, Hall, DetNo to use to find files
        cursor.execute('''SELECT RunNo, Hall, DetNo
            FROM runs NATURAL JOIN muon_rates
            ORDER BY RunNo, Hall, DetNo''')
        run_keys = cursor.fetchall()
        # Look up delayed energy cuts
        cursor.execute('''SELECT Hall, DetNo, Peak - 3 * Resolution,
            Peak + 3 * Resolution
            FROM delayed_energy_fits
            WHERE Source = "nominal"''')
        nominal_energy_bounds = cursor.fetchall()
        cursor.execute('''SELECT Hall, DetNo, Peak - 3 * Resolution,
            Peak + 3 * Resolution
            FROM delayed_energy_fits
            WHERE Source = "adtime"''')
        adtime_energy_bounds = cursor.fetchall()

    energy_lookup = {}
    for site, ad, low_bound, up_bound in nominal_energy_bounds:
        energy_lookup['nominal', site, ad] = (low_bound, up_bound)
    for site, ad, low_bound, up_bound in adtime_energy_bounds:
        energy_lookup['adtime', site, ad] = (low_bound, up_bound)

    with multiprocessing.Pool() as pool:
        results = pool.starmap(one_file, zip(run_keys, it.repeat(data_file_path),
            it.repeat(energy_lookup)))
    if update_db:
        # add Label field
        rows = [(*result, '1/11/21') for result in results]
        with common.get_db(database) as conn:
            cursor = conn.cursor()
            cursor.executemany('''INSERT OR REPLACE INTO adtime_DT_cut_breakdown
                VALUES (?, ?, ?, ?, ?, ?, ?)''', rows)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('--data-file-path',
        help='File path template for run-by-run data files, e.g. '
        '/home/skohn/EH{site}/hadded_ad{ad}/out_ad{ad}_{run}.root'
    )
    parser.add_argument('--update-db', action='store_true')
    args = parser.parse_args()
    main(args.database, args.data_file_path, args.update_db)
