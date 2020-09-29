"""Near-far prediction.

Based on the procedure outlined in DocDB-8774 Section 3.
"""

import argparse
from dataclasses import dataclass
from datetime import datetime
import json
import logging
import pdb
import sqlite3
from typing import Any

import numpy as np

@dataclass
class InputOscParams:
    theta12: float
    m2_21: float
    m2_ee: float

default_osc_params = InputOscParams(33.6469*np.pi/180, 7.53e-5, 2.48e-3)
no_osc_params = InputOscParams(0, 0, 0)
near_ads = [(1, 1), (1, 2), (2, 1), (2, 2)]
far_ads = [(3, 1), (3, 2), (3, 3), (3, 4)]
all_ads = near_ads + far_ads

def ad_dict(initial_value, halls='all'):
    if halls == 'all':
        ads = all_ads
    elif halls == 'near':
        ads = near_ads
    elif halls == 'far':
        ads = far_ads
    else:
        raise ValueError(f'Invalid halls: {halls} (must be "all", "near", or "far")')
    return dict(zip(ads, [initial_value]*len(ads)))

@dataclass
class FitConstants:
    detector_response: dict
    true_bins_response: np.array
    reco_bins_response: np.array
    observed_candidates: dict
    nominal_bgs: dict
    reco_bins: np.array
    total_emitted_by_AD: dict
    true_bins_spectrum: np.array
    input_osc_params: InputOscParams
    muon_eff: dict
    multiplicity_eff: dict

@dataclass
class FitParams:
    theta13: float
    # m2_ee: float
    pull_bg: dict
    pull_near_stat: dict
    @classmethod
    def from_list(cls, in_list):
        theta13 = in_list[0]
        # m2_ee = in_list[1]
        # Background pull parameters
        pull_bg = {}
        for pull, halldet in zip(in_list[1:9], all_ads):
            pull_bg[halldet] = pull
        # Near statistics pull parameters
        pull_near_stat = {}
        for pull, halldet in zip(in_list[9:13], near_ads):
            pull_near_stat[halldet] = pull
        return cls(theta13, pull_bg, pull_near_stat)
    def to_list(self):
        return (
            [self.theta13]
            + [self.pull_bg[halldet] for halldet in all_ads]
            + [self.pull_near_stat[halldet] for halldet in near_ads]
        )

@dataclass
class Config:
    database: str
    period: Any
    backgrounds: Any
    backgrounds_source: str
    mult_eff: Any
    muon_eff: Any
    num_coincs: Any
    num_coincs_source: str
    reco_bins: Any
    det_response_source: str

@dataclass
class ADPeriod:
    name: str
    start_week: int
    end_week: int
    start_run: int
    end_run: int
    start_time_s: int
    end_time_s: int
    live_ads: list

# Hard-code AD Period constants
period_6ad = ADPeriod(
        name="6ad",
        start_week=0,
        end_week=31,  # TODO double-check
        start_run=21221,
        end_run=26694,
        start_time_s=1324684800,
        end_time_s=1344038400,
        live_ads=near_ads[:-1] + far_ads[:-1],
)
period_8ad = ADPeriod(
        name="8ad",
        start_week=42,
        end_week=260,
        start_run=34523,
        end_run=67012,
        start_time_s=1350086400,
        end_time_s=1482537600,
        live_ads=all_ads,
)
period_7ad = ADPeriod(
        name="7ad",
        start_week=265,
        end_week=296,
        start_run=67625,
        end_run=72455,
        start_time_s=1484956800,
        end_time_s=1504310400,
        live_ads=all_ads[1:],
)
ad_periods = {
        "6ad": period_6ad,
        "8ad": period_8ad,
        "7ad": period_7ad,
}


def load_constants(config_file):
    with open(config_file, 'r') as f:
        config_dict = json.load(f)
        config = Config(**config_dict)
    database = config.database
    ad_period = ad_periods[config.period]
    source_det_resp = config.det_response_source
    matrix, true_bins_response, reco_bins_response = true_to_reco_energy_matrix(
            database, source_det_resp
    )
    total_emitted_by_AD, true_bins_spectrum = total_emitted(
            database, slice(ad_period.start_week, ad_period.end_week+1)
    )

    # Parse num coincidences: Nominal, 0, hard-coded, or alternate database
    coincs_source = config.num_coincs_source
    if config.num_coincs is True:
        num_coincidences, reco_bins = num_coincidences_per_AD(database, coincs_source)
    elif config.num_coincs is False:
        num_coincidences = ad_dict(0)
        reco_bins = config.reco_bins
    elif isinstance(config.num_coincs, list):
        coinc_arrays = [np.array(coincs) for coincs in config.num_coincs]
        num_coincidences = dict(zip(all_ads, coinc_arrays))
        reco_bins = config.reco_bins
    elif isinstance(config.num_coincs, str):
        num_coincidences, reco_bins = num_coincidences_per_AD(
                config.num_coincs, coincs_source
        )
    else:
        raise ValueError(
                f"Invalid num coincidences specification: {config.num_coincs}"
        )

    # Parse backgrounds: Nominal, 0, hard-coded, or alternate database
    backgrounds_source = config.backgrounds_source
    if config.backgrounds is True:
        nominal_bgs, reco_bins = backgrounds_per_AD(database, backgrounds_source)
    elif config.backgrounds is False:
        nominal_bgs = ad_dict(0)
    elif isinstance(config.backgrounds, list):
        bg_arrays = [np.array(bg) for bg in config.backgrounds]
        nominal_bgs = dict(zip(all_ads, bg_arrays))
    elif isinstance(config.backgrounds, str):
        nominal_bgs, reco_bins = backgrounds_per_AD(config.backgrounds,
                backgrounds_source)
    else:
        raise ValueError(f"Invalid backgrounds specification: {config.backgrounds}")

    # Parse muon_eff: Nominal, 1, hard-coded, or alternate database
    if config.muon_eff is True:
        muon_eff = muon_veto_efficiency(database)
    elif config.muon_eff is False:
        muon_eff = ad_dict(1)
    elif isinstance(config.muon_eff, list):
        muon_eff = dict(zip(all_ads, config.muon_eff))
    elif isinstance(config.muon_eff, str):
        muon_eff = muon_veto_efficiency(config.muon_eff)
    else:
        raise ValueError(f"Invalid muon_eff specification: {config.muon_eff}")

    # Parse mult_eff: Nominal, 1, hard-coded, or alternate database
    if config.mult_eff is True:
        multiplicity_eff = multiplicity_efficiency(database)
    elif config.mult_eff is False:
        multiplicity_eff = ad_dict(1)
    elif isinstance(config.mult_eff, list):
        multiplicity_eff = dict(zip(all_ads, config.mult_eff))
    elif isinstance(config.mult_eff, str):
        multiplicity_eff = multiplicity_efficiency(config.mult_eff)
    else:
        raise ValueError(f"Invalid mult_eff specification: {config.mult_eff}")

    return FitConstants(
            matrix,
            true_bins_response,
            reco_bins_response,
            num_coincidences,
            nominal_bgs,
            reco_bins,
            total_emitted_by_AD,
            true_bins_spectrum,
            default_osc_params,
            muon_eff,
            multiplicity_eff
    )


def survival_probability(L, E, theta13, m2_ee, input_osc_params=default_osc_params):
    theta12 = input_osc_params.theta12
    m2_21 = input_osc_params.m2_21
    cos4theta13 = np.power(np.cos(theta13), 4)
    sin2_2theta12 = np.power(np.sin(2*theta12), 2)
    sin2_2theta13 = np.power(np.sin(2*theta13), 2)
    delta_21 = delta_ij(L, E, m2_21)
    delta_ee = delta_ij(L, E, m2_ee)
    sin2_delta_21 = np.power(np.sin(delta_21), 2)
    sin2_delta_ee = np.power(np.sin(delta_ee), 2)

    return (
        1 - cos4theta13 * sin2_2theta12 * sin2_delta_21 - sin2_2theta13 * sin2_delta_ee
    )


def delta_ij(L, E, m2_ij):
    return 1.267 * m2_ij * L / E

def muon_veto_efficiency(database):
    """Compute the DAQ-livetime-weighted muon livetime efficiency by AD.

    Return a dict of (hall, det) -> efficiency (as a float).
    """
    run_by_run = {}
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        for hall, det in all_ads:
            cursor.execute('''SELECT Livetime_ns, Efficiency
            FROM muon_rates NATURAL JOIN runs
            WHERE Hall = ? AND DetNo = ?''', (hall, det))
            run_by_run[hall, det] = np.array(cursor.fetchall())
    to_return = {}
    for (hall, det), by_run in run_by_run.items():
        # Weighted average by total livetime.
        # IMPORTANT: `Livetime_ns` is already multiplied by efficiency.
        # True DAQ livetime = Livetime_ns/Efficiency
        total_livetime_ns = by_run[:, 0]/by_run[:, 1]
        efficiency = np.average(by_run[:, 1], weights=total_livetime_ns)
        to_return[hall, det] = efficiency
    return to_return

def multiplicity_efficiency(database):
    """Compute the DAQ-livetime-weighted multiplicity efficiency by AD.

    Return a dict of (hall, det) -> efficiency (as a float).
    """
    run_by_run = {}
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        for hall, det in all_ads:
            cursor.execute('''SELECT Livetime_ns, Efficiency,
                MultiplicityVetoEfficiency
            FROM (muon_rates NATURAL JOIN runs)
                INNER JOIN singles_rates USING (RunNo, DetNo)
            WHERE Hall = ? AND DetNo = ?''', (hall, det))
            run_by_run[hall, det] = np.array(cursor.fetchall())
    to_return = {}
    for (hall, det), by_run in run_by_run.items():
        # Weighted average by total livetime.
        # IMPORTANT: `Livetime_ns` is already multiplied by efficiency.
        # True DAQ livetime = Livetime_ns/Efficiency
        total_livetime_ns = by_run[:, 0]/by_run[:, 1]
        efficiency = np.average(by_run[:, 2], weights=total_livetime_ns)
        to_return[hall, det] = efficiency
    return to_return


def reactor_spectrum(database, core):
    """Returns (spectrum, weekly_time_bins, energy_bins).

    spectrum[i, j] has energy bin i and time/week bin j,
    and units of nuebar/MeV/s.
    """
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT U235, U238, Pu239, Pu241 FROM thermal_energy_per_fission
            WHERE Source = "Phys. Rev. C 88, 014605 (2013)"''')
        energy_per_fission = np.array(cursor.fetchone())
        cursor.execute('''SELECT StartTime_s, EndTime_s, FractionalPower, FracU235, FracU238, FracPu239,
            FracPu241 FROM reactor_power_fissionfrac
            WHERE Core = ? AND Source = "DybBerkFit/Matt/WeeklyAvg_P17B_by_Beda.txt"
            ORDER BY StartTime_s''', (core,))
        power_fissionfrac_history = np.array(cursor.fetchall())
        cursor.execute('''SELECT Energy, U235, U238, Pu239, Pu241 FROM reactor_nuebar_spectrum
            WHERE Source = "DybBerkFit/Matt/C. Lewis/Huber/PRD70, 053011 (2004)"
            ORDER BY Energy''')
        nuebar_spectrum_with_bins = np.array(cursor.fetchall())
        energy_bins = nuebar_spectrum_with_bins[:, 0]
        nuebar_spectrum = nuebar_spectrum_with_bins[:, 1:5]
    weekly_time_bins = power_fissionfrac_history[:, 0:2]
    weekly_power_frac = power_fissionfrac_history[:, 2]
    weekly_fissionfracs = power_fissionfrac_history[:, 3:7]
    MEV_PER_GIGAJOULE = 6.24150907e21
    NOMINAL_POWER_GW = 2.9
    NOMINAL_POWER_MEV_PER_S = NOMINAL_POWER_GW * MEV_PER_GIGAJOULE
    weekly_num_fissions = NOMINAL_POWER_MEV_PER_S * weekly_power_frac / np.sum(weekly_fissionfracs *
            energy_per_fission, axis=1)
    weekly_isotope_spectrum = np.array([
            fissionfracs * nuebar_spectrum for fissionfracs in weekly_fissionfracs
    ])  # indexes: (week, energy, isotope)
    spectrum_per_fission = weekly_isotope_spectrum.sum(axis=2)  # sum over isotopes
    spectrum = weekly_num_fissions * spectrum_per_fission.T
    return spectrum, weekly_time_bins, energy_bins

def livetime_by_week(database, weekly_time_bins):
    """Retrieve a 1D array of # seconds livetime for each week."""
    ordered = {}
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        for hall, det in all_ads:
            cursor.execute('''SELECT Hall, DetNo, Start_time, Livetime_ns/Efficiency
                FROM muon_rates NATURAL JOIN runs
                WHERE Hall = ? AND DetNo = ?
                ORDER BY Start_time''', (hall, det))
            ordered[(hall, det)] = np.array(cursor.fetchall())
    by_week = {}
    for (hall, det), by_run in ordered.items():
        start_times_s = by_run[:, 2]/1e9
        livetimes_s = by_run[:, 3]/1e9
        end_times_s = start_times_s + livetimes_s
        weekly_livetimes_s = np.zeros((len(weekly_time_bins),))
        time_index = 0
        week_start_s, week_end_s = weekly_time_bins[time_index]
        for start_time_s, end_time_s, livetime_s in zip(start_times_s, end_times_s,
                livetimes_s):
            if start_time_s < week_start_s:
                raise ValueError(f'Got out of order with run starting {start_time_s},'
                        f'week starts {week_start_s}')
            if start_time_s > week_end_s:
                while start_time_s > week_end_s and time_index < len(weekly_time_bins):
                    time_index += 1
                    week_start_s, week_end_s = weekly_time_bins[time_index]
                if time_index == len(weekly_time_bins):
                    logging.warn('Reached end of weekly time bins without finishing all runs')
            if end_time_s <= week_end_s:
                weekly_livetimes_s[time_index] += livetime_s
            elif start_time_s < week_end_s and end_time_s > week_end_s:
                in_next_week = end_time_s - week_end_s
                within_current_week = livetime_s - in_next_week
                weekly_livetimes_s[time_index] += within_current_week
                time_index += 1
                week_start_s, week_end_s = weekly_time_bins[time_index]
                weekly_livetimes_s[time_index] += in_next_week

        by_week[(hall, det)] = weekly_livetimes_s
    return by_week

def total_emitted(database, week_range):
    """Return a dict of ((hall, det), core) -> spectrum[i] with energy bin i,
    and an array of energy bins (dict, array as a tuple).

    Units are nuebar/MeV emitted during the livetime of a particular AD.

    ret[0][(hall, det)] corresponds roughly to phi_j in Eq. 1 in DocDB-8774.
    The only difference is that this return value is adjusted
    to compensate for the livetime for each AD
    (and how it overlaps with reactor power, fission fraction, etc.).
    I believe that phi_j should also be adjusted in that way.
    """
    total_spectrum_by_AD = {}
    for core in range(1, 7):
        spectrum, time_bins, energies = reactor_spectrum(database, core)
        by_week = livetime_by_week(database, time_bins)
        for (hall, det), AD_livetime_by_week in by_week.items():
            spectrum_by_week_AD = spectrum * AD_livetime_by_week
            total_spectrum_by_AD[(hall, det), core] = np.sum(spectrum_by_week_AD[:, week_range], axis=1)
    return total_spectrum_by_AD, energies

def cross_section(database):
    """Return (xsec, energy) in units of cm^2, MeV."""
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT Energy, CrossSection FROM xsec
                WHERE Source = "DybBerkFit/toySpectra//Wei on 7/11/2012"
                ORDER BY Energy''')
        values = np.array(cursor.fetchall())
    xsec = values[:, 1]
    energy = values[:, 0]
    return xsec, energy

def xsec_weighted_spec(database):
    """Return a tuple (dict of weighted spectrum, energy bins).

    The dict has keys ((hall, det), core) with core indexed from 1.

    weighted spectrum is a 1D array with units IBDs * cm^2/MeV.
    It needs to be converted to IBDs/MeV by being multiplied by P_sur(L/E)/L^2.
    """
    to_return = {}
    total_spectrum, energy_bins_spec = total_emitted(database, core)
    for ((hall, det), core), spec in total_spectrum.items():
            xsec, energy_bins_xsec = cross_section(database)
            if not np.array_equal(energy_bins_spec, energy_bins_xsec):
                raise ValueError("Energy bins don't line up :(")
            to_return[(hall, det), core] = total_spectrum * xsec
    return to_return, energy_bins_xsec

def flux_fraction(constants, fit_params, week_range=slice(None, None, None)):
    """Return a tuple (dict of flux fractions, energy bins).

    The dict has keys ((hall, det), core) with core indexed from 1.
    """
    to_return = {}
    for (hall, det) in all_ads:
        numerators = np.zeros((len(constants.true_bins_spectrum), 6))
        for core in range(1, 7):
            spec = constants.total_emitted_by_AD[(hall, det), core]
            # Split the label "D1", "L2", etc.
            # into D or L (core_group) and the core index within the group.
            core_group = distance_conversion[core][0]
            core_index = int(distance_conversion[core][1])-1
            distance_m = distances[core_group][core_index][f'EH{hall}'][det-1]
            p_osc = survival_probability(
                    distance_m,
                    constants.true_bins_spectrum,
                    fit_params.theta13,
                    constants.input_osc_params.m2_ee,
                    input_osc_params=constants.input_osc_params
            )
            numerators[:, core-1] = spec * p_osc / distance_m**2
        denominators = np.sum(numerators, axis=1)
        for core in range(1, 7):
            if any(denominators == 0):
                to_return[(hall, det), core] = np.zeros_like(denominators)
            else:
                to_return[(hall, det), core] = numerators[:, core-1] / denominators
    return to_return, constants.true_bins_spectrum

def extrapolation_factor(constants, fit_params):
    """Return a tuple (dict of extrapolation factors, energy bins).

    The dict has keys ((far_hall, far_det), core, (near_hall, near_det))
    with core indexed from 1.
    """
    to_return = {}
    for (near_hall, near_det) in near_ads:
        denominators = np.zeros((len(constants.true_bins_spectrum), 6))
        for core in range(1, 7):
            spec = constants.total_emitted_by_AD[(near_hall, near_det), core]
            # Split the label "D1", "L2", etc.
            # into D or L (core_group) and the core index within the group.
            core_group = distance_conversion[core][0]
            core_index = int(distance_conversion[core][1])-1
            distance_m = distances[core_group][core_index][f'EH{near_hall}'][near_det-1]
            p_osc = survival_probability(
                    distance_m,
                    constants.true_bins_spectrum,
                    fit_params.theta13,
                    constants.input_osc_params.m2_ee,
                    input_osc_params=constants.input_osc_params
            )
            denominators[:, core-1] = spec * p_osc / distance_m**2
            numerators = np.zeros_like(denominators)
            for (far_hall, far_det) in far_ads:
                spec = constants.total_emitted_by_AD[(far_hall, far_det), core]
                # Split the label "D1", "L2", etc.
                # into D or L (core_group) and the core index within the group.
                core_group = distance_conversion[core][0]
                core_index = int(distance_conversion[core][1])-1
                distance_m = distances[core_group][core_index][f'EH{far_hall}'][far_det-1]
                p_osc = survival_probability(
                        distance_m,
                        constants.true_bins_spectrum,
                        fit_params.theta13,
                        constants.input_osc_params.m2_ee,
                        input_osc_params=constants.input_osc_params
                )
                numerators[:, core-1] = spec * p_osc / distance_m**2
                to_return[(far_hall, far_det), core, (near_hall, near_det)] = (
                        numerators[:, core-1]/denominators[:, core-1]
                )
    return to_return, constants.true_bins_spectrum


def num_coincidences_per_AD(database, source):
    """Retrieve the number of observed IBD candidates in each AD, *including backgrounds*.

    Returns a tuple of (dict, energy_bins),
    where the dict maps (hall, det) to a 1D array of N_coincidences.
    This method assumes that all ADs have the same binning.
    """
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT Hall, DetNo, NumCoincidences_binned, BinEdges
        FROM num_coincidences
        WHERE Source = ?
        ORDER BY Hall, DetNo''', (source,))
        full_data = cursor.fetchall()
    results = {}
    # Parse binned records
    for hall, det, coincidence_str, bins_str in full_data:
        results[hall, det] = np.array(json.loads(coincidence_str))
    bins = np.array(json.loads(bins_str), dtype=float)
    # convert keV to MeV
    bins /= 1000
    return results, bins

def fixed_efficiency_weighted_counts(constants, fit_params):
    """Adjust the observed counts by dividing by muon & multiplicity efficiencies.

    Returns a dict mapping (hall, det) to a 1D array of N_coincidences.
    """
    num_coincidences = constants.observed_candidates
    mult_effs = constants.multiplicity_eff
    muon_effs = constants.muon_eff
    pull_near_stat = fit_params.pull_near_stat
    to_return = {}
    for halldet, coincidences in num_coincidences.items():
        mult_eff = mult_effs[halldet]
        muon_eff = muon_effs[halldet]
        # Account for near hall statistics pull parameter
        #if halldet in pull_near_stat:
            #pulled_coincidences = (1 + pull_near_stat[halldet]) * coincidences
        #else:
        pulled_coincidences = coincidences
        to_return[halldet] = pulled_coincidences / (mult_eff * muon_eff)
    return to_return

def backgrounds_per_AD(database, source):
    """Retrieve the number of predicted background events in each AD.

    Returns a tuple of (dict, energy_bins),
    where the dict maps (hall, det) to a 1D array of N_bg_events
    for that AD, over all background types.
    This method assumes that all ADs and background sources
    have the same binning.

    Currently-included backgrounds:

    - Accidentals
    """
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT Hall, DetNo, NumAccs_binned, BinEdges
        FROM num_accs
        WHERE Source = ?
        ORDER BY Hall, DetNo''', (source,))
        acc_data = cursor.fetchall()
    results = {}
    # Parse binned records
    for hall, det, num_accs_str, bins_str in acc_data:
        results[hall, det] = np.array(json.loads(num_accs_str))
    bins = np.array(json.loads(bins_str), dtype=float)
    # convert keV to MeV
    bins /= 1000
    return results, bins

def num_bg_subtracted_IBDs_per_AD(constants, fit_params):
    """Compute the bg-subtracted IBD spectra for each AD.

    Returns the same style values as num_coincidences_per_AD.
    """
    num_candidates = fixed_efficiency_weighted_counts(constants, fit_params)
    predicted_bg = constants.nominal_bgs
    results = {}
    for (hall, det), bg_spec in predicted_bg.items():
        pull = fit_params.pull_bg[hall, det]
        results[hall, det] = num_candidates[hall, det] - bg_spec * (1 + pull)
    return results, constants.reco_bins

def true_to_reco_energy_matrix(database, source):
    """Retrieve the conversion matrix from true to/from reconstructed energy.

    Returns a tuple of (conversion, true_bins, reco_bins).
    The conversion is a 2D array indexed by [true_bin, reco_bin].
    The true and reco bins give the low bin edge value
    for the corresponding axis.
    Note that the true bins are generally taken to be much finer
    than the reco bins.

    There is no particular normalization to the returned array.
    """
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT Matrix, RecoBins, TrueBins
            FROM detector_response
            WHERE Source = ?''', (source,))
        matrix_str, reco_bins_str, true_bins_str = cursor.fetchone()
    # Matrix is stored transposed so we need to un-transpose it
    matrix = np.array(json.loads(matrix_str)).T
    # Divide by 1000 because we store the energies in keV
    reco_bins = np.array(json.loads(reco_bins_str))/1000
    true_bins = np.array(json.loads(true_bins_str))/1000
    return matrix, true_bins, reco_bins

def reco_to_true_energy(constants, fit_params):
    """Apply the detector response to get the "true" observed spectrum.

    This is handled independently for each reconstructed energy bin,
    so the returned spectrum arrays are 2-dimensional.

    Returns a tuple of (true_spectrum_dict, true_bins, reco_bins).
    The units of the spectrum are IBDs per MeV.
    The dict maps (hall, det) to a 2D array with indexes [true_index, reco_index].
    """
    # Normalize the matrix so that each column sums to 1.
    # (Thus each reco event will be spread among true energy bins
    # in a unitary manner).
    weights_per_reco_bin = constants.detector_response.sum(axis=0, keepdims=True)
    normalized_response = constants.detector_response/weights_per_reco_bin
    if not np.array_equal(constants.reco_bins_response, constants.reco_bins):
        print(constants.reco_bins_response)
        print(constants.reco_bins)
        raise NotImplementedError("Different reco binnings")
    true_energies = {}
    bg_subtracted, bins = num_bg_subtracted_IBDs_per_AD(constants, fit_params)
    for (hall, det), num_IBDs in bg_subtracted.items():
        true_energies[hall, det] = normalized_response * np.expand_dims(num_IBDs, axis=0)
    return true_energies, constants.true_bins_response, constants.reco_bins_response

def num_IBDs_from_core(constants, fit_params):
    """Compute the number of IBDs in each AD that come from a given core.

    Return a dict of ((hall, det), core) -> N_ij, and the binnings,
    as a tuple (N_ij_dict, true_bins, reco_bins).
    (N_ij from Eq. 2 of DocDB-8774.)
    N_ij is a 2D array with index [true_index, reco_index]
    with bins of energy given by the second return value.
    """
    flux_fractions, true_bins_spec = flux_fraction(constants, fit_params)
    true_energies, true_bins_resp, reco_bins = reco_to_true_energy(constants, fit_params)
    # reactor flux bins don't include the upper edge of 12 MeV
    true_bins_spec = np.concatenate((true_bins_spec, [12]))
    n_ij = {}
    for (halldet, core), flux_frac in flux_fractions.items():
        rebinned_fluxfrac = average_bins(flux_frac, true_bins_spec, true_bins_resp)
        # Must expand dims so that the shape of the arrays is
        # (N_true, 1) * (N_true, N_reco)
        n_ij[halldet, core] = np.expand_dims(rebinned_fluxfrac, axis=1) * true_energies[halldet]
    return n_ij, true_bins_resp, reco_bins

def predict_IBD_true_energy(constants, fit_params):
    """Compute the predicted number of IBDs for a given pair of ADs, by true and reco
    energy.

    Return a tuple (f_kji, true_bins, reco_bins) where f_kji is a dict of
    ((far_hall, far_det), core, (near_hall, near_det)) -> N,
    with core indexed from 1 and N[true_index, reco_index].
    """
    num_from_core, true_bins, reco_bins = num_IBDs_from_core(constants, fit_params)
    extrap_factor, true_bins_spec = extrapolation_factor(constants, fit_params)
    # reactor flux bins don't include the upper edge of 12 MeV
    true_bins_spec = np.concatenate((true_bins_spec, [12]))
    f_kji = {}
    for (far_halldet, core, near_halldet), extrap_fact in extrap_factor.items():
        rebinned_extrap_fact = average_bins(extrap_fact, true_bins_spec, true_bins)
        n_ij = num_from_core[near_halldet, core]
        f_kji[far_halldet, core, near_halldet] = (
                np.expand_dims(rebinned_extrap_fact, axis=1) * n_ij
        )
    return f_kji, true_bins, reco_bins

def predict_ad_to_ad_IBDs(constants, fit_params):
    """Compute the predicted number of IBDs for a given pair of ADs, summed over all
    cores.

    Return a tuple (f_ki, reco_bins) where f_ki is a dict of
    ((far_hall, far_det), (near_hall, near_det)) -> N,
    and N indexed by reco_index.
    """
    f_kji, true_bins, reco_bins = predict_IBD_true_energy(constants, fit_params)
    f_ki = {}
    for (far_halldet, core, near_halldet), n in f_kji.items():
        if (far_halldet, near_halldet) in f_ki:
            f_ki[far_halldet, near_halldet] += n.sum(axis=0)
        else:
            f_ki[far_halldet, near_halldet] = n.sum(axis=0)
    return f_ki, reco_bins

def predict_ad_to_ad_obs(constants, fit_params):
    """Compute the number of observed coincidences for a given pair of ADs.

    The background events are added back to the far halls
    and the muon-veto and multiplicity-veto efficiencies are re-applied.

    Return a tuple (f_ki, reco_bins) where f_ki is a dict of
    ((far_hall, far_det), (near_hall, near_det)) -> N,
    and N indexed by reco_index.

    The return value of this function should be directly comparable
    to the number of observed coincidences at the far hall ADs,
    without adjusting the observed counts for backgrounds or any efficiencies.
    """
    f_ki, reco_bins = predict_ad_to_ad_IBDs(constants, fit_params)
    predicted_bg = constants.nominal_bgs
    muon_effs = constants.muon_eff
    mult_effs = constants.multiplicity_eff
    results = f_ki.copy()
    for (far_hall, far_det), near_halldet in f_ki:
        bg_for_far = predicted_bg[far_hall, far_det]
        pull = fit_params.pull_bg[far_hall, far_det]
        results[(far_hall, far_det), near_halldet] += bg_for_far * (1 + pull)
    for (far_halldet, near_halldet), result in results.items():
        muon_eff = muon_effs[far_halldet]
        mult_eff = mult_effs[far_halldet]
        results[far_halldet, near_halldet] = result * muon_eff * mult_eff
    return results, constants.reco_bins


def predict_halls(constants, fit_params):
    """Compute the predicted number of IBDs in EH3 based on EH1 or EH2.

    Return a tuple (f_pred, reco_bins) where f_pred is a dict with keys
    1 and 2 corresponding to EH1 and EH2, and values of a 1D array
    indexed by reco_index.

    There are 8 AD pairs from a near hall to the far hall,
    and therefore 8 predictions.
    The predictions are all summed to represent combining the far-hall ADs
    and then halved to represent averaging over the 2 near AD predictions.
    """
    f_ki, reco_bins = predict_ad_to_ad_IBDs(constants, fit_params)
    prediction = {
            1: np.zeros_like(reco_bins[:-1]),
            2: np.zeros_like(reco_bins[:-1])
    }
    for (far_halldet, (near_hall, near_det)), ad_prediction in f_ki.items():
        prediction[near_hall] += ad_prediction/2
    return prediction, reco_bins


def average_bins(values_fine, bins_fine, bins_coarse):
    """Average the values from the fine bins into a coarser binning.

    Weighted average with weight = bin width.

    Assumptions:

    - ``len(values_fine) + 1 == len(bins_fine)``
    - ``bins_fine[0] == bins_coarse[0]``
    - ``bins_fine[-1] == bins_coarse[-1]``
    - ``len(bins_coarse) < len(bins_fine)``
    - bins_coarse is a subset of bins_fine (i.e. no fine bins cross a coarse bin edge)
    - bins_coarse and bins_fine are strictly increasing
    - and that all inputs are 1D arrays.
    """
    # Test assumptions
    #test_bins = len(values_fine) + 1 == len(bins_fine)
    #test_start_bin = np.isclose(bins_fine[0], bins_coarse[0])
    #test_end_bin = np.isclose(bins_fine[-1], bins_coarse[-1])
    #test_fine_really_fine = len(bins_coarse) < len(bins_fine)
    #test_no_crossings = True  # Trust the user!
    bin_weights = np.diff(bins_fine)
    #test_increasing = np.all(bin_weights > 0)
    #test_1D = True  # Trust the user!
    #if not np.all([
        #test_bins, test_start_bin, test_end_bin,
        #test_fine_really_fine, test_no_crossings,
        #test_increasing, test_1D
    #]):
        #print(bins_fine[0], bins_coarse[0])
        #print(bins_fine[-1], bins_coarse[-1])
        #raise ValueError("Invalid inputs. Check assumptions!")
    values_coarse = np.empty((len(bins_coarse)-1,))
    fine_index = 0
    fine_up_edge = bins_fine[0]  # initial value only for first comparison in while loop
    for coarse_index, coarse_up_edge in enumerate(bins_coarse[1:]):
        next_coarse_value = 0
        next_coarse_weights_sum = 0
        #while not np.isclose(fine_up_edge, coarse_up_edge):
        while fine_up_edge != coarse_up_edge:
            #if fine_up_edge > coarse_up_edge:
                #print(fine_up_edge, coarse_up_edge)
                #return
            fine_up_edge = bins_fine[fine_index + 1]
            fine_value = values_fine[fine_index]
            fine_weight = bin_weights[fine_index]
            next_coarse_value += fine_value * fine_weight
            next_coarse_weights_sum += fine_weight
            fine_index += 1
        # Then we can finalize the bin
        next_coarse_value /= next_coarse_weights_sum
        values_coarse[coarse_index] = next_coarse_value
    return values_coarse






distance_conversion = {
    1: 'D1',
    2: 'D2',
    3: 'L1',
    4: 'L2',
    5: 'L3',
    6: 'L4',
}


distances = {
        'D': [{
            'EH1': [362.38, 357.94],
            'EH2': [1332.48, 1337.43],
            'EH3': [1919.63, 1917.52, 1925.26, 1923.15]
            }, {
            'EH1': [371.76, 368.41],
            'EH2': [1358.15, 1362.88],
            'EH3': [1894.34, 1891.98, 1899.86, 1897.51]
            }],
        'L': [{
            'EH1': [903.47, 903.35],
            'EH2': [467.57, 472.97],
            'EH3': [1533.18, 1534.92, 1538.93, 1540.67]
            }, {
            'EH1': [817.16, 816.90],
            'EH2': [489.58, 495.35],
            'EH3': [1533.63, 1535.03, 1539.47, 1540.87]
            }, {
            'EH1': [1353.62, 1354.23],
            'EH2': [557.58, 558.71],
            'EH3': [1551.38, 1554.77, 1556.34, 1559.72]
            }, {
            'EH1': [1265.32, 1265.89],
            'EH2': [499.21, 501.07],
            'EH3': [1524.94, 1528.05, 1530.08, 1533.18]
            }]
        }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config')
    parser.add_argument('-d', '--debug', action='store_true')
    args = parser.parse_args()
    constants = load_constants(args.config)
    fit_params = FitParams(0.15, 2.5e-3, ad_dict(0), ad_dict(0))
    print(predict_halls(constants, fit_params))
