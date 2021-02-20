"""Near-far prediction.

Based on the procedure outlined in DocDB-8774 Section 3.
"""

import argparse
from collections import defaultdict
from dataclasses import dataclass, field as dc_field
from datetime import datetime
import itertools
import json
import logging
import pdb
import sqlite3
from typing import Any, ClassVar

import numpy as np
import scipy.interpolate

import common
import compute_dataset_summary

@dataclass
class InputOscParams:
    theta12: float
    m2_21: float
    hierarchy: bool = True  # True for normal, False for inverted
    m2_ee_conversion: float = 5.17e-5  # Conversion from m2_ee to m2_32

default_osc_params = InputOscParams(33.6469*np.pi/180, 7.53e-5, True, 5.17e-5)
no_osc_params = InputOscParams(0, 0, True, 0)
near_ads = [(1, 1), (1, 2), (2, 1), (2, 2)]
far_ads = [(3, 1), (3, 2), (3, 3), (3, 4)]
all_ads = near_ads + far_ads
sites = [1, 2, 3]
cores = range(1, 7)

def ad_dict(initial_value, halls='all', factory=False):
    if halls == 'all':
        ads = all_ads
    elif halls == 'near':
        ads = near_ads
    elif halls == 'far':
        ads = far_ads
    else:
        raise ValueError(f'Invalid halls: {halls} (must be "all", "near", or "far")')
    if factory:
        val_array = [initial_value() for _ in ads]
    else:
        val_array = [initial_value for _ in ads]
    return dict(zip(ads, val_array))

def core_dict(initial_value):
    return dict(zip(range(1, 7), [initial_value] * 6))

def site_dict(initial_value):
    val_array = [initial_value] * len(sites)
    return dict(zip(sites, val_array))

class RebinCache:
    """A cache to save time in rebinning flux fractor and extrap factor.

    These matrix multiplications / rebins are the slowest operations in
    the prediction calculation, by far.
    """
    def __init__(self):
        self.flux_fractions = {}
        self.extrap_factors = {}
        return

    def _save_flux_frac(self, orig_hash, rebinned):
        self.flux_fractions[orig_hash] = rebinned
        return

    def _save_extrap_fact(self, orig_hash, rebinned):
        self.extrap_factors[orig_hash] = rebinned
        return

    def flux_frac(self, orig, rebin_matrix):
        """Compute or retrieve the rebinned flux fraction.

        Relies on the cache and if the original flux fraction
        has already been rebinned, will return the saved version.
        If the original is new, then this method will
        compute the rebinning, save it in the internal cache,
        and return it.
        """
        orig_hash = hash(tuple(orig))
        if orig_hash not in self.flux_fractions:
            self._save_flux_frac(
                orig_hash,
                np.matmul(rebin_matrix, orig)
            )
        return self.flux_fractions[orig_hash]

    def extrap_fact(self, orig, rebin_matrix):
        """Compute or retrieve the rebinned extrapolation factor.

        Relies on the cache and if the original extrapolation factor
        has already been rebinned, will return the saved version.
        If the original is new, then this method will
        compute the rebinning, save it in the internal cache,
        and return it.
        """
        orig_hash = hash(tuple(orig))
        if orig_hash not in self.extrap_factors:
            self._save_extrap_fact(
                orig_hash,
                np.matmul(rebin_matrix, orig)
            )
        return self.extrap_factors[orig_hash]


@dataclass
class FitConstants:
    detector_response: dict
    true_bins_response: np.array
    reco_bins_response: np.array
    observed_candidates: dict
    nominal_bgs: dict  # keys = different types of bg, values = dict by AD/Hall
    bg_errors: dict  # keys = different types of bg, values = dict by AD/Hall
    reco_bins: np.array
    total_emitted_by_AD: dict
    total_emitted_by_week_AD: dict
    livetime_by_week_AD: dict
    daq_livetimes: dict
    true_bins_spectrum: np.array
    input_osc_params: InputOscParams
    theta12_err: float
    m2_21_err: float
    m2_ee_err: float
    muon_eff: dict
    multiplicity_eff: dict
    dist_time_eff: dict
    efficiency_err: float
    reactor_err: float
    rel_escale_err: float
    masses: dict
    standard_mass: float
    cross_section: np.array
    rebin_matrix: np.array
    lbnl_comparison: bool
    rel_escale_parameters: dict
    prompt_eff_baseline_corrs: dict
    _rebin_cache: RebinCache

@dataclass
class FitParams:
    """All the parameters that can be adjusted during the fit process.

    Attributes
    ----------
    theta13
        Not sin-squared, but the actual value, in radians.
    m2_ee
    pull_bg
    pull_near_stat
    pull_reactor
    pull_efficiency
        This pull parameter does *not* include the uncertainties due to
        prompt energy cut, which is fully correlated with
        the ``pull_rel_escale`` (relative energy scale) parameter.
    pull_rel_escale
        Pull parameter for relative energy scale uncertainty.
        This parameter impacts both the shape of the prompt spectrum
        and the efficiency (via prompt energy cut).
    """
    num_pulls: ClassVar = 6 + 8 + 8 + 8 + 3 + 3 + 1 + 8 + 8 + 1 + 1 + 1 + 4 * 34  # TODO binning
    theta13: float
    m2_ee: float
    pull_reactor: dict = dc_field(default_factory=lambda: core_dict(0))
    pull_efficiency: dict = dc_field(default_factory=lambda: ad_dict(0))
    pull_rel_escale: dict = dc_field(default_factory=lambda: ad_dict(0))
    pull_accidental: dict = dc_field(default_factory=lambda: ad_dict(0))
    pull_li9: dict = dc_field(default_factory=lambda: site_dict(0))
    pull_fast_neutron: dict = dc_field(default_factory=lambda: site_dict(0))
    pull_amc: float = 0
    pull_alpha_n: dict = dc_field(default_factory=lambda: ad_dict(0))
    pull_rad_n: dict = dc_field(default_factory=lambda: ad_dict(0))
    pull_theta12: float = 0
    pull_m2_21: float = 0
    pull_m2_ee: float = 0
    pull_near_stat: dict = dc_field(
        # initialize an ad_dict with independent numpy arrays
        default_factory=lambda: ad_dict(lambda: np.zeros(34), halls='near', factory=True)
        # TODO binning
    )
    @property
    def sin2_2theta13(self):
        return np.power(np.sin(2 * self.theta13), 2)

    @sin2_2theta13.setter
    def sin2_2theta13(self, new):
        self.theta13 = 0.5 * np.arcsin(np.sqrt(new))

    def clone(self):
        return FitParams.from_list(self.to_list())

    @staticmethod
    def index_map():
        to_return = {}
        to_return['theta13'] = 0
        to_return['m2_ee'] = 1
        to_return['reactor'] = slice(2, 8)
        to_return['efficiency'] = slice(8, 16)
        to_return['rel-escale'] = slice(16, 24)
        to_return['accidental'] = slice(24, 32)
        to_return['li9'] = slice(32, 35)
        to_return['fast-neutron'] = slice(35, 38)
        to_return['amc'] = 38
        to_return['alpha-n'] = slice(39, 47)
        to_return['rad-n'] = slice(47, 55)
        to_return['pull_theta12'] = 55
        to_return['pull_m2_21'] = 56
        to_return['pull_m2_ee'] = 57
        n_bins = 34  # TODO hard-coded
        n_near_ads = 4
        num_near_stat = n_bins * n_near_ads
        to_return['near-stat'] = slice(58, 58 + num_near_stat)
        return to_return

    @classmethod
    def from_list(cls, in_list):
        indexes = cls.index_map()
        theta13 = in_list[indexes['theta13']]
        m2_ee = in_list[indexes['m2_ee']]
        # Reactor pull parameters
        pull_reactor = {}
        for i, pull in enumerate(in_list[indexes['reactor']]):
            pull_reactor[i + 1] = pull
        # Detection efficiency pull parameters
        pull_efficiency = {}
        for pull, halldet in zip(in_list[indexes['efficiency']], all_ads):
            pull_efficiency[halldet] = pull
        # Relative energy scale pull parameters
        pull_rel_escale = {}
        for pull, halldet in zip(in_list[indexes['rel-escale']], all_ads):
            pull_rel_escale[halldet] = pull
        # Background pull parameters
        pull_acc = {}
        for pull, halldet in zip(in_list[indexes['accidental']], all_ads):
            pull_acc[halldet] = pull
        pull_li9 = {}
        for pull, site in zip(in_list[indexes['li9']], sites):
            pull_li9[site] = pull
        pull_fastn = {}
        for pull, site in zip(in_list[indexes['fast-neutron']], sites):
            pull_fastn[site] = pull
        pull_amc = in_list[indexes['amc']]
        pull_alpha_n = {}
        for pull, halldet in zip(in_list[indexes['alpha-n']], all_ads):
            pull_alpha_n[halldet] = pull
        pull_rad_n = {}
        for pull, halldet in zip(in_list[indexes['rad-n']], all_ads):
            pull_rad_n[halldet] = pull
        pull_theta12 = in_list[indexes['pull_theta12']]
        pull_m2_21 = in_list[indexes['pull_m2_21']]
        pull_m2_ee = in_list[indexes['pull_m2_ee']]
        # Near statistics pull parameters
        near_stat_slice = indexes['near-stat']
        first_entry = near_stat_slice.start
        n_near_pulls = near_stat_slice.stop - first_entry
        n_bins = n_near_pulls // 4
        pull_near_stat = {}
        last_entry = first_entry + n_bins
        for halldet in near_ads:
            pull_near_stat[halldet] = np.array(in_list[first_entry:last_entry])
            first_entry += n_bins
            last_entry += n_bins
        return cls(
            theta13,
            m2_ee,
            pull_reactor,
            pull_efficiency,
            pull_rel_escale,
            pull_acc,
            pull_li9,
            pull_fastn,
            pull_amc,
            pull_alpha_n,
            pull_rad_n,
            pull_theta12,
            pull_m2_21,
            pull_m2_ee,
            pull_near_stat,
        )

    def to_list(self):
        return (
            [self.theta13, self.m2_ee]
            + [self.pull_reactor[core] for core in range(1, 7)]
            + [self.pull_efficiency[halldet] for halldet in all_ads]
            + [self.pull_rel_escale[halldet] for halldet in all_ads]
            + [self.pull_accidental[halldet] for halldet in all_ads]
            + [self.pull_li9[site] for site in sites]
            + [self.pull_fast_neutron[site] for site in sites]
            + [self.pull_amc]
            + [self.pull_alpha_n[halldet] for halldet in all_ads]
            + [self.pull_rad_n[halldet] for halldet in all_ads]
            + [self.pull_theta12, self.pull_m2_21, self.pull_m2_ee]
            + [x for halldet in near_ads for x in self.pull_near_stat[halldet]]
        )

    def pull_bg(self, name, all=False):
        """Lookup the corresponding pull param dict using a string.

        If all is True, then name is ignored and a dict is returned
        that contains all background types.
        """
        bg_dict = {
            'accidental': self.pull_accidental,
            'li9': self.pull_li9,
            'fast-neutron': self.pull_fast_neutron,
            'amc': self.pull_amc,
            'alpha-n': self.pull_alpha_n,
            'rad-n': self.pull_rad_n,
        }
        if all:
            return bg_dict
        else:
            return bg_dict[name]

    def pulled_osc_params(self, nominal_params):
        """Return an InputOscParams object with the pulled parameter values."""
        pulled_theta12 = nominal_params.theta12 * (1 + self.pull_theta12)
        pulled_m2_21 = nominal_params.m2_21 * (1 + self.pull_m2_21)
        hierarchy = nominal_params.hierarchy
        m2_ee_conversion = nominal_params.m2_ee_conversion
        return InputOscParams(pulled_theta12, pulled_m2_21, hierarchy, m2_ee_conversion)



@dataclass
class Config:
    database: str
    period: Any
    backgrounds: Any
    background_counts_source: str
    background_spectra_source: str
    background_types: list
    background_errors: list
    mult_eff: Any
    mult_eff_label: str
    muon_eff: Any
    muon_eff_label: str
    dist_time_eff: Any
    dist_time_eff_label: str
    masses: Any
    num_coincs: Any
    num_coincs_source: str
    num_coincs_format_version: int
    reco_bins: Any
    det_response_source: str
    rel_escale_source: str
    sin2_2theta12_error: Any
    m2_21_error: Any
    m2_ee_relative_error: Any
    binning_id: int
    efficiency_error: float
    reactor_error: float
    rel_escale_error: float
    prompt_eff_corr_source: str
    lbnl_comparison: bool = False


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
    # total_emitted_by_AD, true_bins_spectrum, total_emitted_by_week_AD = total_emitted(
            # database, slice(ad_period.start_week, ad_period.end_week+1)
    # )
    # Sum over all periods TODO
    total_emitted_by_AD = defaultdict(int)
    for period in ad_periods.values():
        total_emitted_for_period, true_bins_spectrum = total_emitted_shortcut(
            database, period.name
        )
        for halldetcore, spec in total_emitted_for_period.items():
            total_emitted_by_AD[halldetcore] += spec
    _, time_bins, _ = reactor_spectrum(database, 1)
    # LBNL comparison removes the AD-to-AD livetime dependence from the reactor flux
    if config.lbnl_comparison:
        livetime_by_week_by_AD = livetime_by_week(database, time_bins)
        livetime_by_AD_for_periods = {
                (halldet, period.name):
                sum(livetime_by_week_by_AD[halldet][
                    period.start_week:period.end_week+1])
                for period in [period_6ad, period_8ad, period_7ad]
                for halldet in all_ads
        }
        for key in total_emitted_by_AD:
            halldet, core = key
            total_emitted_by_AD[key] /= livetime_by_AD_for_periods[halldet,
                    ad_period.name]

    livetimes_by_week = livetime_by_week(database, time_bins)
    daq_livetimes = dict(zip(all_ads, compute_dataset_summary.daq_livetime_s(
        database if not isinstance(config.muon_eff, str) else config.muon_eff,
        config.muon_eff_label,
    )))
    rebin_matrix = generate_bin_averaging_matrix(
            np.concatenate((true_bins_spectrum, [12])),
            true_bins_response
    )

    # Parse num coincidences: Nominal, 0, hard-coded, or alternate database
    coincs_source = config.num_coincs_source
    coincs_format = config.num_coincs_format_version
    if config.num_coincs is True:
        num_coincidences, reco_bins = num_coincidences_per_AD(database, coincs_source,
            version=coincs_format, binning_id=config.binning_id)
    elif config.num_coincs is False:
        num_coincidences = ad_dict(0)
        reco_bins = config.reco_bins
    elif isinstance(config.num_coincs, list):
        coinc_arrays = [np.array(coincs) for coincs in config.num_coincs]
        num_coincidences = dict(zip(all_ads, coinc_arrays))
        reco_bins = config.reco_bins
    elif isinstance(config.num_coincs, str):
        num_coincidences, reco_bins = num_coincidences_per_AD(
            config.num_coincs, coincs_source, version=coincs_format,
            binning_id=config.binning_id
        )
    else:
        raise ValueError(
                f"Invalid num coincidences specification: {config.num_coincs}"
        )

    # Parse backgrounds: Nominal, 0, hard-coded, or alternate database
    bg_counts_source = config.background_counts_source
    bg_spec_source = config.background_spectra_source
    bg_types = config.background_types
    if config.backgrounds is True:
        nominal_bgs, bg_errors = backgrounds_per_AD(
            database, bg_counts_source, bg_spec_source, types=bg_types
        )
    elif config.backgrounds is False:
        nominal_bgs, bg_errors = backgrounds_per_AD(
            None, None, None, return_zero=True, types=bg_types
        )
    elif isinstance(config.backgrounds, list):
        # List of lists to get each AD's spectrum
        # List of floats to get each AD's absolute error (then divide by sum(count) for
        # relative)
        bg_arrays = [np.array(bg) for bg in config.backgrounds]
        bg_errors_input = [
            err / sum(count) for err, count in zip(config.background_errors, bg_arrays)
        ]
        nominal_bgs, bg_errors = backgrounds_per_AD(
            None, None, None, return_zero=True, types=bg_types[0]
        )
        # Assign supplied backgrounds to be the first type given in background_types
        nominal_bgs[bg_types[0]] = dict(zip(all_ads, bg_arrays))
        bg_errors[bg_types[0]] = dict(zip(all_ads, bg_errors_input))
    elif isinstance(config.backgrounds, str):
        nominal_bgs, bg_errors = backgrounds_per_AD(
            config.backgrounds,
            bg_counts_source,
            bg_spec_source,
            types=bg_types,
        )
    else:
        raise ValueError(f"Invalid backgrounds specification: {config.backgrounds}")

    def zip_ads(list_of_values):
        return dict(zip(all_ads, list_of_values))

    # Parse muon_eff: Nominal, 1, hard-coded, or alternate database
    if config.muon_eff is True:
        muon_eff = zip_ads(
            compute_dataset_summary.muon_efficiency(
                database, config.muon_eff_label
            )
        )
    elif config.muon_eff is False:
        muon_eff = ad_dict(1)
    elif isinstance(config.muon_eff, list):
        muon_eff = zip_ads(config.muon_eff)
    elif isinstance(config.muon_eff, str):
        muon_eff = zip_ads(
            compute_dataset_summary.muon_efficiency(
                config.muon_eff, config.muon_eff_label
            )
        )
    else:
        raise ValueError(f"Invalid muon_eff specification: {config.muon_eff}")

    # Parse mult_eff: Nominal, 1, hard-coded, or alternate database
    if config.mult_eff is True:
        multiplicity_eff = zip_ads(
            compute_dataset_summary.multiplicity_efficiency(
                database, config.mult_eff_label
            )
        )
    elif config.mult_eff is False:
        multiplicity_eff = ad_dict(1)
    elif isinstance(config.mult_eff, list):
        multiplicity_eff = zip_ads(config.mult_eff)
    elif isinstance(config.mult_eff, str):
        multiplicity_eff = zip_ads(
            compute_dataset_summary.multiplicity_efficiency(
                config.mult_eff, config.mult_eff_label
            )
        )
    else:
        raise ValueError(f"Invalid mult_eff specification: {config.mult_eff}")

    # Parse dist_time_eff: Nominal, 1, hard-coded, or alternate database
    if config.dist_time_eff is True:
        dist_time_eff = distance_time_efficiency(database, config.mult_eff_label)
    elif config.dist_time_eff is False:
        dist_time_eff = ad_dict(1)
    elif isinstance(config.dist_time_eff, list):
        dist_time_eff = config.dist_time_eff
    elif isinstance(config.dist_time_eff, str):
        dist_time_eff = distance_time_efficiency(
            config.dist_time_eff, config.dist_time_eff_label
        )
    else:
        raise ValueError(f"Invalid dist_time_eff specification: {config.dist_time_eff}")

    masses = dict(zip(all_ads, config.masses))

    cross_sec = cross_section(database)

    rel_escale_params = rel_escale_parameters(database, config.rel_escale_source)

    input_osc_params = default_osc_params

    # Convert from error on sin2_2theta21 (easy to input)
    # to relative error on theta21 (easy to compute with)
    sin2_2theta12 = np.sin(2 * input_osc_params.theta12)**2
    sin2_2theta12_up = sin2_2theta12 + config.sin2_2theta12_error
    sin2_2theta12_low = sin2_2theta12 - config.sin2_2theta12_error
    theta12_up = 0.5 * np.arcsin(np.sqrt(sin2_2theta12_up))
    theta12_low = 0.5 * np.arcsin(np.sqrt(sin2_2theta12_low))
    theta12_err = max(
        theta12_up - input_osc_params.theta12, input_osc_params.theta12 - theta12_low
    )
    theta12_rel_err = theta12_err/input_osc_params.theta12

    # Compute m2_21 relative error
    m2_21_rel_err = config.m2_21_error/input_osc_params.m2_21

    prompt_eff_baseline_corrs = load_prompt_eff_baseline_corr(
        database, config.prompt_eff_corr_source, config.binning_id
    )

    return FitConstants(
            matrix,
            true_bins_response,
            reco_bins_response,
            num_coincidences,
            nominal_bgs,
            bg_errors,
            reco_bins,
            total_emitted_by_AD,
            #total_emitted_by_week_AD,
            None,
            livetimes_by_week,
            daq_livetimes,
            true_bins_spectrum,
            input_osc_params,
            theta12_rel_err,
            m2_21_rel_err,
            config.m2_ee_relative_error,
            muon_eff,
            multiplicity_eff,
            dist_time_eff,
            config.efficiency_error,
            config.reactor_error,
            config.rel_escale_error,
            masses,
            masses[1, 1],
            cross_sec,
            rebin_matrix,
            config.lbnl_comparison,
            rel_escale_params,
            prompt_eff_baseline_corrs,
            RebinCache(),
    )


def survival_probability(L, E, theta13, m2_ee, input_osc_params=default_osc_params):
    theta12 = input_osc_params.theta12
    m2_21 = input_osc_params.m2_21
    # Hierarchy factor = 2 * int(hierarchy) - 1 (+1 if True, else -1)
    hierarchy = 2 * int(input_osc_params.hierarchy) - 1
    m2_32 = m2_ee - hierarchy * input_osc_params.m2_ee_conversion
    m2_31 = m2_32 + hierarchy * m2_21

    cos4theta13 = np.power(np.cos(theta13), 4)
    cos2theta12 = np.power(np.cos(theta12), 2)
    sin2theta12 = np.power(np.sin(theta12), 2)
    sin2_2theta12 = np.power(np.sin(2*theta12), 2)
    sin2_2theta13 = np.power(np.sin(2*theta13), 2)
    delta_21 = delta_ij(L, E, m2_21)
    delta_31 = delta_ij(L, E, m2_31)
    delta_32 = delta_ij(L, E, m2_32)

    return (
        1
        - cos4theta13 * sin2_2theta12 * delta_21
        - cos2theta12 * sin2_2theta13 * delta_31
        - sin2theta12 * sin2_2theta13 * delta_32
    )


def delta_ij(L, E, m2_ij):
    phase_arg = 1.267 * m2_ij * L / E
    sin_val = np.sin(phase_arg)
    return sin_val * sin_val

def muon_veto_efficiency(database):
    """Compute the DAQ-livetime-weighted muon livetime efficiency by AD.

    Return a dict of (hall, det) -> efficiency (as a float).
    """
    run_by_run = {}
    with common.get_db(database) as conn:
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
    with common.get_db(database) as conn:
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

def distance_time_efficiency(database, source):
    """Retrieve the distance-time (DT) efficiency."""
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''
            SELECT
                Efficiency
            FROM
                distance_time_cut_efficiency
            WHERE
                Source = ?
            ORDER BY
                Hall,
                DetNo
            ''',
            (source,)
        )
        result = np.array(cursor.fetchall()).reshape(-1)
    return dict(zip(all_ads, result))

def rel_escale_parameters(database, source):
    """Load the shape distortion parameters for the relative energy scale.

    Return a dict keyed by hall
    whose values are arrays with rows corresponding to reco bins and columns
    corresponding to 0: distortion for plus, 1: distortion for minus.
    The scale of the distortions is normalized so that the factor
    (1 + rel_escale_offset * rel_escale_parameters) gives the appropriate
    correction if rel_escale_offset is e.g. 0.001 for 0.1%.
    """
    result = {}
    halls = [1, 2, 3]
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        for hall in halls:
            cursor.execute('''
                SELECT
                    CoefficientPlus,
                    CoefficientMinus,
                    NominalUncertainty
                FROM
                    rel_energy_scale_shape
                WHERE
                        BinningId = 1
                    AND
                        Hall = ?
                    AND
                        Source = ?
                ORDER BY
                    BinIndex
                ''',
                (hall, source)
            )
            params = np.array(cursor.fetchall())
            result[hall] = params[:, :2]/params[:, 2:]  # divide by last column
    return result


def _prompt_eff_baseline_corr_zeros(sin2, dm2):
    """Return only zeros."""
    return 0

def load_prompt_eff_baseline_corr(database, source, binning_id):
    """Load the corrections to prompt efficiency due to baseline / oscillation.

    Return a dict keyed by (halldet, core) whose values are
    scipy.interpolate.BivariateSpline objects which, given (sin2_2theta13, dm2_ee),
    produce the approximate correction for the given baseline.

    >>> corrections_for_baseline = result[(1, 2), 3]  # EH1-AD2 core 3
    >>> specific_correction = corrections_for_baseline(0.734, 0.00273)  # sin2, dm2
    >>> corrected_number = num_events/(1 + specific_correction)

    To convert from observed counts to comparable counts,
    divide by (1 + correction).
    """
    result = {}
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        for hall, det in all_ads:
            for core in cores:
                if source is None:
                    result[(hall, det), core] = _prompt_eff_baseline_corr_zeros
                    continue
                cursor.execute('''
                    SELECT
                        SinSq2T13
                    FROM
                        prompt_eff_osc_corrections
                    WHERE
                            BinningId = ?
                        AND
                            Source = ?
                        AND
                            Hall = ?
                        AND
                            DetNo = ?
                        AND
                            Core = ?
                    GROUP BY
                        SinSq2T13
                    ORDER BY
                        SinSq2T13
                    ''',
                    (binning_id, source, hall, det, core)
                )
                sin2_values = np.array(cursor.fetchall()).reshape(-1)
                cursor.execute('''
                    SELECT
                        DM2ee
                    FROM
                        prompt_eff_osc_corrections
                    WHERE
                            BinningId = ?
                        AND
                            Source = ?
                        AND
                            Hall = ?
                        AND
                            DetNo = ?
                        AND
                            Core = ?
                    GROUP BY
                        DM2ee
                    ORDER BY
                        DM2ee
                    ''',
                    (binning_id, source, hall, det, core)
                )
                m2_ee_values = np.array(cursor.fetchall()).reshape(-1)
                cursor.execute('''
                    SELECT
                        SinSq2T13,
                        DM2ee,
                        Correction
                    FROM
                        prompt_eff_osc_corrections
                    WHERE
                            BinningId = ?
                        AND
                            Source = ?
                        AND
                            Hall = ?
                        AND
                            DetNo = ?
                        AND
                            Core = ?
                    ORDER BY
                        SinSq2T13,
                        DM2ee
                    ''',
                    (binning_id, source, hall, det, core)
                )
                rows = np.array(cursor.fetchall())
                # Validate that there is a (dense) grid in parameter space
                pairs_in_rows = set(tuple(pair) for pair in rows[:, :2])
                for grid_pair in itertools.product(sin2_values, m2_ee_values):
                    if grid_pair not in pairs_in_rows:
                        raise ValueError(
                            f"Couldn't find {grid_pair} in data for "
                            f"EH{hall}-AD{det} Core {core} prompt eff corrections"
                        )
                # For each (halldet, core) pair, construct an interpolator
                gridded_corrections = rows[:, 2].reshape((len(sin2_values),
                    len(m2_ee_values)))
                interpolator = scipy.interpolate.RectBivariateSpline(
                    sin2_values,
                    m2_ee_values,
                    gridded_corrections,
                    kx=1, ky=1,  # interpolation degree = 1 ==> linear
                )
                result[(hall, det), core] = interpolator
    return result


def reactor_spectrum(database, core):
    """Returns (spectrum, weekly_time_bins, energy_bins).

    spectrum[i, j] has energy bin i and time/week bin j,
    and units of nuebar/MeV/s.
    """
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        # cursor.execute('''SELECT U235, U238, Pu239, Pu241 FROM thermal_energy_per_fission
            # WHERE Source = "Phys. Rev. C 88, 014605 (2013)"''')
        cursor.execute('''SELECT U235, U238, Pu239, Pu241 FROM thermal_energy_per_fission
            WHERE Source = "DybBerkFit / DocDB-7580"''')
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
        cursor.execute('''SELECT U235, U238, Pu239, Pu241 FROM nonequilibrium
            WHERE Source = "DybBerkFit make_combined_spectra_P17B_unblinded.C"
            ORDER BY Energy''')
        nonequilibrium = np.array(cursor.fetchall())
        cursor.execute('''SELECT CorrectionFactor FROM snf_correction
            WHERE Source = "DybBerkFit make_combined_spectra_P17B_unblinded.C"
            AND Core = ?
            ORDER BY Energy''', (core,))
        snf_correction = np.array(cursor.fetchall())
    nuebar_spectrum = nuebar_spectrum * nonequilibrium * snf_correction
    weekly_time_bins = power_fissionfrac_history[:, 0:2]
    weekly_power_frac = power_fissionfrac_history[:, 2]
    weekly_fissionfracs = power_fissionfrac_history[:, 3:7]
    MEV_PER_GIGAJOULE = 6.24150907e21
    NOMINAL_POWER_GW = 2.95
    NOMINAL_POWER_MEV_PER_S = NOMINAL_POWER_GW * MEV_PER_GIGAJOULE
    weekly_num_fissions = NOMINAL_POWER_MEV_PER_S * weekly_power_frac / np.sum(weekly_fissionfracs *
            energy_per_fission, axis=1)
    weekly_isotope_spectrum = np.array([
            fissionfracs * nuebar_spectrum for fissionfracs in weekly_fissionfracs
    ])  # indexes: (week, energy, isotope)
    spectrum_per_fission = weekly_isotope_spectrum.sum(axis=2)  # sum over isotopes
    spectrum = weekly_num_fissions * spectrum_per_fission.T
    return spectrum, weekly_time_bins, energy_bins

def livetime_for_period(weekly_livetimes, period):
    """Sum the appropriate weeks' livetimes for the given data period.

    weekly_livetimes should be a 1D array of livetime for each week
    for the desired AD.

    period should be an ADPeriod object.
    """
    return sum(weekly_livetimes[period.start_week : period.end_week+1])

def livetime_by_week(database, weekly_time_bins):
    """Retrieve a 1D array of # seconds livetime for each week."""
    ordered = {}
    with common.get_db(database) as conn:
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
        # TODO mc_corrections
    mc_corrections = {(1, 1): 4.445342402310752e-06,
 (1, 2): 4.5343298147423175e-06,
 (2, 1): 1.0546587931834448e-05,
 (2, 2): 1.0523020263115215e-05,
 (3, 1): 5.837502639759078e-05,
 (3, 2): 5.8327610278713265e-05,
 (3, 3): 6.88640992131794e-05,
 (3, 4): 5.8191618278611e-05}
    for halldet in by_week:
        by_week[halldet] = by_week[halldet] * (1 + mc_corrections[halldet])

    return by_week

def total_emitted_shortcut(database, data_period):
    """Same as total_emitted but use saved database values
    and don't break down by week."""

    total_spectrum_by_AD = {}
    _, time_bins, _ = reactor_spectrum(database, 1)
    livetime_by_week_by_AD = livetime_by_week(database, time_bins)
    livetime_by_AD_for_periods = {
            (halldet, period.name):
            sum(livetime_by_week_by_AD[halldet][period.start_week:period.end_week+1])
            for period in [period_6ad, period_8ad, period_7ad]
            for halldet in all_ads
    }
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        for core in range(1, 7):
            cursor.execute('''SELECT Energy, NuPerMeVPerSec
            FROM reactor_emitted_spectrum
            WHERE Core = ? AND DataPeriod = ?
                AND Source = "DybBerkFit/ReactorPowerCalculator/isotope_spectra_by_Beda"
                AND Energy <= 10
            ORDER BY Energy''',
            (core, data_period))
            result = np.array(cursor.fetchall())
            spectrum = result[:, 1]
            spectrum_mid_bins = spectrum[:-1] + 0.5 * np.diff(spectrum)
            for halldet in all_ads:
                total_spectrum_by_AD[halldet, core] = (
                    spectrum_mid_bins
                    * livetime_by_AD_for_periods[halldet, data_period]
                )

            energies = result[:-1, 0]  # Last entry is upper bin boundary
    return total_spectrum_by_AD, energies


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
    spectrum_by_week_AD = {}
    for core in range(1, 7):
        spectrum, time_bins, energies = reactor_spectrum(database, core)
        by_week = livetime_by_week(database, time_bins)
        for (hall, det), AD_livetime_by_week in by_week.items():
            spectrum_by_week_this_AD = spectrum * AD_livetime_by_week
            spectrum_by_week_AD[(hall, det), core] = spectrum_by_week_this_AD
            total_spectrum_by_AD[(hall, det), core] = np.sum(
                    spectrum_by_week_this_AD[:, week_range], axis=1
            )
    return total_spectrum_by_AD, energies, spectrum_by_week_AD

def cross_section(database):
    """Return (xsec, energy) in units of cm^2, MeV."""
    with common.get_db(database) as conn:
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
    total_spectrum, energy_bins_spec, _ = total_emitted(database, core)
    for ((hall, det), core), spec in total_spectrum.items():
            xsec, energy_bins_xsec = cross_section(database)
            if not np.array_equal(energy_bins_spec, energy_bins_xsec):
                raise ValueError("Energy bins don't line up :(")
            to_return[(hall, det), core] = total_spectrum * xsec
    return to_return, energy_bins_xsec


def flux_fraction(constants, fit_params, week_range=slice(None, None, None),
        include_osc=True):
    """Return a dict of flux fractions for near ADs.

    The dict has keys ((hall, det), core) with core indexed from 1.
    """
    to_return = {}
    get_to_bin_centers_hack = np.diff(constants.true_bins_spectrum)[0]/2
    pulled_osc_params = fit_params.pulled_osc_params(constants.input_osc_params)
    for (hall, det) in near_ads:
        numerators = np.zeros((len(constants.true_bins_spectrum), 6))
        for core in range(1, 7):
            if week_range == slice(None):
                spec = constants.total_emitted_by_AD[(hall, det), core]
            else:
                spec = np.sum(constants.total_emitted_by_week_AD[(hall, det), core][
                        :, week_range], axis=1)
            # Apply the reactor pull parameter
            pull = fit_params.pull_reactor[core]
            spec = (1 + pull) * spec
            distance_m = distances[core][f'EH{hall}'][det-1]
            if include_osc:
                p_osc = survival_probability(
                        distance_m,
                        constants.true_bins_spectrum + get_to_bin_centers_hack,
                        fit_params.theta13,
                        fit_params.m2_ee * (1 + fit_params.pull_m2_ee),
                        input_osc_params=pulled_osc_params
                )
            else:
                p_osc = 1
            numerators[:, core-1] = spec * p_osc / distance_m**2
        denominators = np.sum(numerators, axis=1)
        for core in range(1, 7):
            if any(denominators == 0):
                to_return[(hall, det), core] = np.zeros_like(denominators)
            else:
                to_return[(hall, det), core] = numerators[:, core-1] / denominators
    return to_return

def extrapolation_factor(constants, fit_params):
    """Return a dict of extrapolation factors.

    The dict has keys ((far_hall, far_det), core, (near_hall, near_det))
    with core indexed from 1.
    """
    to_return = {}
    far_osc_probs = {}
    get_to_bin_centers_hack = np.diff(constants.true_bins_spectrum)[0]/2
    pulled_osc_params = fit_params.pulled_osc_params(constants.input_osc_params)
    for (near_hall, near_det) in near_ads:
        denominators = np.zeros((len(constants.true_bins_spectrum), 6))
        for core in range(1, 7):
            spec = constants.total_emitted_by_AD[(near_hall, near_det), core]
            # Apply the reactor pull parameter
            pull = fit_params.pull_reactor[core]
            spec = (1 + pull) * spec
            distance_m = distances[core][f'EH{near_hall}'][near_det-1]
            p_osc = survival_probability(
                    distance_m,
                    constants.true_bins_spectrum + get_to_bin_centers_hack,
                    fit_params.theta13,
                    fit_params.m2_ee * (1 + fit_params.pull_m2_ee),
                    input_osc_params=pulled_osc_params
            )
            denominators[:, core-1] = spec * p_osc / distance_m**2
            numerators = np.zeros_like(denominators)
            for (far_hall, far_det) in far_ads:
                spec = constants.total_emitted_by_AD[(far_hall, far_det), core]
                # Apply the reactor pull parameter
                pull = fit_params.pull_reactor[core]
                spec = (1 + pull) * spec
                distance_m = distances[core][f'EH{far_hall}'][far_det-1]
                if ((far_hall, far_det), core) in far_osc_probs:
                    p_osc = far_osc_probs[(far_hall, far_det), core]
                else:
                    p_osc = survival_probability(
                            distance_m,
                            constants.true_bins_spectrum + get_to_bin_centers_hack,
                            fit_params.theta13,
                            fit_params.m2_ee * (1 + fit_params.pull_m2_ee),
                            input_osc_params=pulled_osc_params
                    )
                    far_osc_probs[(far_hall, far_det), core] = p_osc
                numerators[:, core-1] = spec * p_osc / distance_m**2
                to_return[(far_hall, far_det), core, (near_hall, near_det)] = (
                        numerators[:, core-1]/denominators[:, core-1]
                )
    return to_return


def num_coincidences_per_AD(database, source, version, binning_id):
    """Retrieve the number of observed IBD candidates in each AD, *including backgrounds*.

    Returns a tuple of (dict, energy_bins),
    where the dict maps (hall, det) to a 1D array of N_coincidences.
    This method assumes that all ADs have the same binning.

    Version 0 has the entire histogram for an AD in 1 row
    stored as a JSON-ified list.

    Version 1 has each bin as a row.
    """
    if version == 0:
        with common.get_db(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''SELECT Hall, DetNo, NumCoincidences_binned
            FROM num_coincidences
            WHERE Source = ? AND BinningId = ?
            ORDER BY Hall, DetNo''', (source, binning_id))
            full_data = cursor.fetchall()
            cursor.execute('''SELECT BinEdgeEnergy_keV FROM reco_binnings
            WHERE Id = ?
            ORDER BY BinEdgeIndex''', (binning_id,))
            bins = cursor.fetchall()  #  [(b0,), (b1,), ...]
        results = {}
        # Parse binned records
        for hall, det, coincidence_str in full_data:
            results[hall, det] = np.array(json.loads(coincidence_str))
        bins = np.array(bins, dtype=float).reshape(-1)  # flatten bins
        # convert keV to MeV
        bins /= 1000
    return results, bins


def efficiency_weighted_counts(constants, fit_params, halls='all'):
    """Adjust the observed counts by dividing by efficiencies.

    This includes the muon and multiplicity efficiencies,
    detection efficiency pull parameter, target masses,
    and changes due to relative energy scale variation.
    The actual absolute detection efficiency is not included since it is common to all ADs.
    The change in prompt efficiency due to oscillation effects is applied in
    num_near_IBDs_from_core since it varies by core.

    Returns a dict mapping (hall, det) to a 1D array of N_coincidences.
    """
    num_coincidences = num_bg_subtracted_IBDs_per_AD(constants, fit_params, halls=halls)
    mult_effs = constants.multiplicity_eff
    muon_effs = constants.muon_eff
    dist_time_effs = constants.dist_time_eff
    mass_effs = constants.masses
    pull_near_stat = fit_params.pull_near_stat
    to_return = {}
    for halldet, coincidences in num_coincidences.items():
        mult_eff = mult_effs[halldet]
        muon_eff = muon_effs[halldet]
        dist_time_eff = dist_time_effs[halldet]
        mass_eff = mass_effs[halldet] / constants.standard_mass
        pull_eff = fit_params.pull_efficiency[halldet]
        pull_rel_escale = fit_params.pull_rel_escale[halldet]
        hall, det = halldet
        if pull_rel_escale > 0:
            rel_escale_params = constants.rel_escale_parameters[hall][:, 0]
        else:
            rel_escale_params = constants.rel_escale_parameters[hall][:, 1]
        # NB: combined_eff has a bin dependence due to rel_escale_params
        combined_eff = (
            mult_eff * muon_eff * dist_time_eff * mass_eff * (1 + pull_eff)
            * (1 + pull_rel_escale * rel_escale_params)
        )
        # Account for near hall statistics pull parameter
        if halldet in pull_near_stat:
            pulled_coincidences = (1 + pull_near_stat[halldet]) * coincidences
        else:
            pulled_coincidences = coincidences
        to_return[halldet] = pulled_coincidences / combined_eff
    return to_return


def backgrounds_per_AD(
    database, counts_source, spectra_sources, return_zero=False, types=None
):
    """Retrieve the number of predicted background events and rate errors in each AD.

    Returns a tuple of (spectra, errors).
    Spectra is a dict of (bg_name -> dict of (hall, det) -> 1D array).
    Errors is a dict of (bg_name -> dict of (hall, det) ->
    *relative* error on the normalization of the corresponding
    spectrum) (i.e. 0.1 = 10%).
    This method assumes that all ADs and background sources
    have the same binning.

    Currently-included backgrounds:

    - Accidentals
    - Li9
    - Fast-neutron
    - AmC
    - Alpha-n

    To get a dict with the correct format but 0 background
    and a nominal error of 1, set return_zero=True.
    The error is set to 1 so that in the chi-square expression
    we can divide by the error without dividing by zero.

    The types parameter, if provided, limits the returned dict to only
    certain types of background. The list can contain any of:

     - "accidental"
     - "li9"
     - "fast-neutron"
     - "amc"
     - "alpha-n"
     - "rad-n" (radiogenic neutrons)

    If types is None, then the default of all types is used.

    If spectra_sources is a str, use it for all types of background.
    If it is a list, then for each type of background, try using
    successive sources until the first one that has more than 0 entries.
    So if one spectrum is updated, use spectra_sources = [newer, older, oldest].
    """
    if return_zero:
        if types is None:
            return (
                {
                    'accidental': ad_dict(0),
                    'li9': ad_dict(0),
                    'fast-neutron': ad_dict(0),
                    'amc': ad_dict(0),
                    'alpha-n': ad_dict(0),
                    'rad-n': ad_dict(0),
                },
                {
                    'accidental': 1,
                    'li9': 1,
                    'fast-neutron': 1,
                    'amc': 1,
                    'alpha-n': 1,
                    'rad-n': 1,
                },
            )
        else:
            return (
                {bg: ad_dict(0) for bg in types},
                {bg: 1 for bg in types},
            )
    if isinstance(spectra_sources, str):
        spectra_sources = [spectra_sources]
    result = {}
    errors = {}
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        if 'accidental' in types:
            acc_data = {}
            acc_errors = {}
            for hall, det in all_ads:
                for spectra_source in spectra_sources:
                    cursor.execute('''
                        SELECT
                            Spectrum
                        FROM
                            accidentals_spectrum
                        WHERE
                            Label = ?
                            AND Hall = ?
                            AND DetNo = ?
                        ORDER BY BinIndex
                        ''',
                        (spectra_source, hall, det),
                    )
                    acc_data[hall, det] = np.array(cursor.fetchall()).reshape(-1)  # flatten
                    if len(acc_data[hall, det]) == 0:
                        continue
                    else:
                        break
                else:  # nobreak
                    raise ValueError('No valid sources given for accidentals: '
                        + str(spectra_sources)
                    )
                cursor.execute('''
                SELECT
                    Count, Error
                FROM
                    bg_counts
                WHERE
                    BgName = "accidental"
                    AND Label = ?
                    AND Hall = ?
                    AND DetNo = ?
                ''', (counts_source, hall, det)
                )
                count, error = cursor.fetchone()
                acc_data[hall, det] *= count
                acc_errors[hall, det] = error/count
            result['accidental'] = acc_data
            errors['accidental'] = acc_errors
        if 'li9' in types:
            for spectra_source in spectra_sources:
                cursor.execute('''
                    SELECT
                        Spectrum
                    FROM
                        li9_spectrum
                    WHERE
                        Label = ?
                    ORDER BY BinIndex
                    ''',
                    (spectra_source,),
                )
                li9_spectrum = np.array(cursor.fetchall()).reshape(-1)
                if len(li9_spectrum) == 0:
                    continue
                else:
                    break
            else:  # nobreak
                raise ValueError('No valid sources given for li9: '
                    + str(spectra_sources)
                )
            li9_data = {}
            li9_errors = {}
            for hall, det in all_ads:
                cursor.execute('''
                SELECT
                    Count, Error
                FROM
                    bg_counts
                WHERE
                    BgName = "li9"
                    AND Label = ?
                    AND Hall = ?
                    AND DetNo = ?
                ''', (counts_source, hall, det)
                )
                count, error = cursor.fetchone()
                li9_data[hall, det] = count * li9_spectrum
                li9_errors[hall, det] = error/count
            result['li9'] = li9_data
            errors['li9'] = li9_errors
        if 'fast-neutron' in types:
            for spectra_source in spectra_sources:
                cursor.execute('''
                    SELECT
                        Spectrum
                    FROM
                        fast_neutron_spectrum
                    WHERE
                        Label = ?
                    ORDER BY BinIndex
                    ''',
                    (spectra_source,),
                )
                fast_neutron_spectrum = np.array(cursor.fetchall()).reshape(-1)
                if len(fast_neutron_spectrum) == 0:
                    continue
                else:
                    break
            else:  # nobreak
                raise ValueError('No valid sources given for fast neutron: '
                    + str(spectra_sources)
                )
            fast_neutron_data = {}
            fast_neutron_errors = {}
            for hall, det in all_ads:
                cursor.execute('''
                SELECT
                    Count, Error
                FROM
                    bg_counts
                WHERE
                    BgName = "fast-neutron"
                    AND Label = ?
                    AND Hall = ?
                    AND DetNo = ?
                ''', (counts_source, hall, det)
                )
                count, error = cursor.fetchone()
                fast_neutron_data[hall, det] = count * fast_neutron_spectrum
                fast_neutron_errors[hall, det] = error/count
            result['fast-neutron'] = fast_neutron_data
            errors['fast-neutron'] = fast_neutron_errors
        if 'amc' in types:
            for spectra_source in spectra_sources:
                cursor.execute('''
                    SELECT
                        Spectrum
                    FROM
                        amc_spectrum
                    WHERE
                        Label = ?
                    ORDER BY BinIndex
                    ''',
                    (spectra_source,),
                )
                amc_spectrum = np.array(cursor.fetchall()).reshape(-1)
                if len(amc_spectrum) == 0:
                    continue
                else:
                    break
            else:  # nobreak
                raise ValueError('No valid sources given for amc: '
                    + str(spectra_sources)
                )
            amc_data = {}
            amc_errors = {}
            for hall, det in all_ads:
                cursor.execute('''
                SELECT
                    Count, Error
                FROM
                    bg_counts
                WHERE
                    BgName = "amc"
                    AND Label = ?
                    AND Hall = ?
                    AND DetNo = ?
                ''', (counts_source, hall, det)
                )
                count, error = cursor.fetchone()
                amc_data[hall, det] = count * amc_spectrum
                amc_errors[hall, det] = error/count
            result['amc'] = amc_data
            errors['amc'] = amc_errors
        if 'alpha-n' in types:
            alpha_n_data = {}
            alpha_n_errors = {}
            for hall, det in all_ads:
                for spectra_source in spectra_sources:
                    cursor.execute('''
                        SELECT
                            Spectrum
                        FROM
                            alpha_n_spectrum
                        WHERE
                            Label = ?
                            AND Hall = ?
                            AND DetNo = ?
                        ORDER BY BinIndex
                        ''',
                        (spectra_source, hall, det),
                    )
                    alpha_n_data[hall, det] = np.array(cursor.fetchall()).reshape(-1)
                    if len(alpha_n_data[hall, det]) == 0:
                        continue
                    else:
                        break
                else:  # nobreak
                    raise ValueError('No valid sources given for alpha-n: '
                        + str(spectra_sources)
                    )
                cursor.execute('''
                SELECT
                    Count, Error
                FROM
                    bg_counts
                WHERE
                    BgName = "alpha-n"
                    AND Label = ?
                    AND Hall = ?
                    AND DetNo = ?
                ''', (counts_source, hall, det)
                )
                count, error = cursor.fetchone()
                alpha_n_data[hall, det] *= count
                alpha_n_errors[hall, det] = error/count
            result['alpha-n'] = alpha_n_data
            errors['alpha-n'] = alpha_n_errors
        if 'rad-n' in types:
            rad_n_data = {}
            rad_n_errors = {}
            for hall, det in all_ads:
                for spectra_source in spectra_sources:
                    cursor.execute('''
                        SELECT
                            Spectrum
                        FROM
                            rad_n_spectrum
                        WHERE
                            Label = ?
                            AND Hall = ?
                            AND DetNo = ?
                        ORDER BY BinIndex
                        ''',
                        (spectra_source, hall, det),
                    )
                    rad_n_data[hall, det] = np.array(cursor.fetchall()).reshape(-1)
                    if len(rad_n_data[hall, det]) == 0:
                        continue
                    else:
                        break
                else:  # nobreak
                    raise ValueError('No valid sources given for rad-n: '
                        + str(spectra_sources)
                    )
                cursor.execute('''
                SELECT
                    Count, Error
                FROM
                    bg_counts
                WHERE
                    BgName = "rad-n"
                    AND Label = ?
                    AND Hall = ?
                    AND DetNo = ?
                ''', (counts_source, hall, det)
                )
                count, error = cursor.fetchone()
                rad_n_data[hall, det] *= count
                rad_n_errors[hall, det] = error/count
            result['rad-n'] = rad_n_data
            errors['rad-n'] = rad_n_errors
    return (result, errors)


def bg_with_pulls(constants, fit_params, halls='all'):
    """Apply the pull parameters to each background.

    This is non-trivial because some pull parameters are per-AD
    but others are per-site or totally correlated across all ADs.

    The halls parameter could also be 'near' or 'far'.
    """
    if halls == 'all':
        ads = all_ads
    elif halls == 'near':
        ads = near_ads
    elif halls == 'far':
        ads = far_ads
    nominal_bg = constants.nominal_bgs
    pulled_bg = {}
    if 'accidental' in nominal_bg:
        # Uncorrelated between ADs
        bg_specs = nominal_bg['accidental']
        bg_pulls = fit_params.pull_accidental
        bg_dict = {}
        for halldet in ads:
            bg_spec = bg_specs[halldet]
            pull = bg_pulls[halldet]
            bg_dict[halldet] = bg_spec * (1 + pull)
        pulled_bg['accidental'] = bg_dict
    if 'li9' in nominal_bg:
        # Correlated within site, uncorrelated between sites
        bg_specs = nominal_bg['li9']
        bg_pulls = fit_params.pull_li9
        bg_dict = {}
        for (hall, det) in ads:
            bg_spec = bg_specs[hall, det]
            pull = bg_pulls[hall]
            bg_dict[hall, det] = bg_spec * (1 + pull)
        pulled_bg['li9'] = bg_dict
    if 'fast-neutron' in nominal_bg:
        # Correlated within site, uncorrelated between sites
        bg_specs = nominal_bg['fast-neutron']
        bg_pulls = fit_params.pull_fast_neutron
        bg_dict = {}
        for (hall, det) in ads:
            bg_spec = bg_specs[hall, det]
            pull = bg_pulls[hall]
            bg_dict[hall, det] = bg_spec * (1 + pull)
        pulled_bg['fast-neutron'] = bg_dict
    if 'amc' in nominal_bg:
        # Correlated among all ADs and sites
        bg_specs = nominal_bg['amc']
        pull = fit_params.pull_amc
        bg_dict = {}
        for halldet in ads:
            bg_spec = bg_specs[halldet]
            bg_dict[halldet] = bg_spec * (1 + pull)
        pulled_bg['amc'] = bg_dict
    if 'alpha-n' in nominal_bg:
        # Uncorrelated between ADs
        bg_specs = nominal_bg['alpha-n']
        bg_pulls = fit_params.pull_alpha_n
        bg_dict = {}
        for halldet in ads:
            bg_spec = bg_specs[halldet]
            pull = bg_pulls[halldet]
            bg_dict[halldet] = bg_spec * (1 + pull)
        pulled_bg['alpha-n'] = bg_dict
    if 'rad-n' in nominal_bg:
        # Uncorrelated between ADs
        bg_specs = nominal_bg['rad-n']
        bg_pulls = fit_params.pull_rad_n
        bg_dict = {}
        for halldet in ads:
            bg_spec = bg_specs[halldet]
            pull = bg_pulls[halldet]
            bg_dict[halldet] = bg_spec * (1 + pull)
        pulled_bg['rad-n'] = bg_dict
    return pulled_bg


def num_bg_subtracted_IBDs_per_AD(constants, fit_params, halls='all'):
    """Compute the bg-subtracted IBD spectra for each AD.

    Returns a dict mapping (hall, det) to a 1D array of IBD counts.

    Currently applies the same pull parameter to each background type.

    Assumes that the background counts are the raw number of events
    contaminating the observed number of counts.
    This means e.g. the number of Li9 events that survived the muon veto
    and were not multiplicity-vetoed.
    These counts will be corrected for these effects *later*.
    """
    num_candidates = constants.observed_candidates
    predicted_bg = constants.nominal_bgs
    pulled_bg = bg_with_pulls(constants, fit_params, halls=halls)
    total_bg = ad_dict(0, halls=halls)
    for bg_type, bg_dict in pulled_bg.items():
        for halldet, bg_spec in bg_dict.items():
            total_bg[halldet] += bg_spec
    results = {}
    for halldet, bg in total_bg.items():
        results[halldet] = num_candidates[halldet] - bg
    return results


def true_to_reco_energy_matrix(database, source):
    """Retrieve the conversion matrix from true to/from reconstructed energy.

    Returns a tuple of (conversion, true_bins, reco_bins).
    The conversion is a 2D array indexed by [true_bin, reco_bin].
    The true and reco bin arrays give the bin edge values
    for the corresponding axis, and are of length nbins+1
    (where nbins is the number of bins for the given axis).
    Note that the true bins are generally taken to be much finer
    than the reco bins.

    There is no particular normalization to the returned array.
    """
    with common.get_db(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT Matrix, RecoBins, TrueBins
            FROM detector_response
            WHERE Source = ?''', (source,))
        matrix_str, reco_bins_str, true_bins_str = cursor.fetchone()
    # Matrix is stored transposed so we need to un-transpose it
    drm_matrix = np.array(json.loads(matrix_str)).T
    # Divide by 1000 because we store the energies in keV
    reco_bins = np.array(json.loads(reco_bins_str))/1000
    true_bins = np.array(json.loads(true_bins_str))/1000
    return drm_matrix, true_bins, reco_bins

def true_to_reco_near_matrix_with_osc(constants, fit_params):
    """Compute the normalized detector response matrix including oscillations.

    Returns a dict mapping (near_hall, near_det) -> a 2D array drm, with
    drm[true_bin, reco_bin] being the probability that an IBD in the
    reco_bin was caused by an antineutrino with an energy in true_bin.

    For a given bin of reconstructed energy,
    oscillations are applied by splitting into the
    approximate reactor contributions using the flux fraction
    (assuming no oscillation)
    and then multiplying by the survival probability for that baseline.
    The decomposed oscillated spectra are then added back together
    and normalized so that the total integral of the reconstructed energy bin
    is unity.

    The original implementation of this procedure is in the LBNL fitter
    Prediction.cc file, MakePrediction() first few dozen lines.
    """
    true_bin_fluxfrac_centers = (
        constants.true_bins_spectrum[:-1]
        + 0.5 * np.diff(constants.true_bins_spectrum)
    )
    #TODO hardcoded last bin...standardize that bins should include upper edge!
    true_bin_fluxfrac_centers = np.concatenate((
        true_bin_fluxfrac_centers,
        [9.975]
    ))
    flux_fractions_no_osc = flux_fraction(constants, fit_params, include_osc=False)
    decomposed_response = {}
    un_normalized_response = {halldet: 0 for halldet in near_ads}
    pulled_osc_params = fit_params.pulled_osc_params(constants.input_osc_params)
    for ((hall, det), core), flux_frac_no_osc in flux_fractions_no_osc.items():
        distance_m = distances[core][f'EH{hall}'][det-1]
        p_sur = survival_probability(
            distance_m,
            true_bin_fluxfrac_centers,
            fit_params.theta13,
            fit_params.m2_ee * (1 + fit_params.pull_m2_ee),
            input_osc_params=pulled_osc_params
        )
        weighted_fluxfrac = flux_frac_no_osc * p_sur
        rebinned_fluxfrac = constants._rebin_cache.flux_frac(weighted_fluxfrac,
                constants.rebin_matrix)
        decomposed_response[(hall, det), core] = (
            constants.detector_response * np.expand_dims(rebinned_fluxfrac, axis=1)
        )
        un_normalized_response[hall, det] += decomposed_response[(hall, det), core]
    # Normalize the matrix so that each column sums to 1.
    # (Thus each reco event will be spread among true energy bins
    # in a unitary manner).
    response_pdfs = {}
    for halldet, un_normalized in un_normalized_response.items():
        weights_per_reco_bin = un_normalized.sum(axis=0, keepdims=True)
        response_pdfs[halldet] = un_normalized/weights_per_reco_bin
    if not np.array_equal(constants.reco_bins_response, constants.reco_bins):
        print(constants.reco_bins_response)
        print(constants.reco_bins)
        raise NotImplementedError("Different reco binnings")
    return response_pdfs

def true_energy_of_near_IBDs(constants, fit_params):
    """Apply the detector response to get the "true" observed spectrum for near ADs.

    This is handled independently for each reconstructed energy bin,
    so the returned spectrum arrays are 2-dimensional.

    Returns a dict of the true spectrum.
    The units of the spectrum are IBDs per MeV.
    The dict maps (hall, det) to a 2D array with indexes [true_index, reco_index].
    The binning is given by constants.true_bins_response and
    constants.reco_bins_response.

    """
    true_energies = {}
    bg_subtracted = efficiency_weighted_counts(constants, fit_params, halls='near')
    response_pdfs = true_to_reco_near_matrix_with_osc(constants, fit_params)
    for (hall, det), response_pdf in response_pdfs.items():
        true_energies[hall, det] = (
            response_pdf * np.expand_dims(bg_subtracted[hall, det], axis=0)
        )
    return true_energies

def num_near_IBDs_from_core(constants, fit_params):
    """Compute the number of IBDs in each near AD that come from a given core.

    Return a dict of ((hall, det), core) -> N_ij.
    (N_ij from Eq. 2 of DocDB-8774.)
    N_ij is a 2D array with index [true_index, reco_index]

    This method also applies the prompt efficiency baseline/oscillation
    correction for near ADs.
    """
    true_bins_spec = constants.true_bins_spectrum
    # reactor flux bins don't include the upper edge of 12 MeV
    true_bins_spec = np.concatenate((true_bins_spec, [12]))
    true_bins_resp = constants.true_bins_response
    reco_bins = constants.reco_bins
    flux_fractions = flux_fraction(constants, fit_params)
    true_spectra = true_energy_of_near_IBDs(constants, fit_params)
    n_ij = {}
    for (halldet, core), flux_frac in flux_fractions.items():
        rebinned_fluxfrac = constants._rebin_cache.flux_frac(flux_frac,
                constants.rebin_matrix)
        # Must expand dims so that the shape of the arrays is
        # (N_true, 1) * (N_true, N_reco)
        n_ij[halldet, core] = np.expand_dims(rebinned_fluxfrac, axis=1) * true_spectra[halldet]
        n_ij[halldet, core] /= (1 + constants.prompt_eff_baseline_corrs[halldet, core](
            fit_params.sin2_2theta13, fit_params.m2_ee * (1 + fit_params.pull_m2_ee)
        ))
    return n_ij

def predict_far_IBD_true_energy(constants, fit_params):
    """Compute the predicted number of far IBDs for a given pair of ADs, by true and reco
    energy.

    Return f_kji, a dict of
    ((far_hall, far_det), core, (near_hall, near_det)) -> N,
    with core indexed from 1 and N[true_index, reco_index].

    This method also applies the prompt efficiency baseline/oscillation
    correction for far ADs.
    """
    true_bins = constants.true_bins_response
    reco_bins = constants.reco_bins
    num_from_core = num_near_IBDs_from_core(constants, fit_params)
    extrap_factors = extrapolation_factor(constants, fit_params)
    f_kji = {}
    for (far_halldet, core, near_halldet), extrap_fact in extrap_factors.items():
        rebinned_extrap_fact = constants._rebin_cache.extrap_fact(extrap_fact,
                constants.rebin_matrix)
        n_ij = num_from_core[near_halldet, core]
        f_kji[far_halldet, core, near_halldet] = (
                np.expand_dims(rebinned_extrap_fact, axis=1) * n_ij
        )
        f_kji[far_halldet, core, near_halldet] *= (
            1 + constants.prompt_eff_baseline_corrs[far_halldet, core](
                fit_params.sin2_2theta13, fit_params.m2_ee * (1 + fit_params.pull_m2_ee)
            )
        )
    return f_kji

def predict_ad_to_ad_IBDs(constants, fit_params):
    """Compute the predicted number of IBDs for a given pair of ADs, summed over all
    cores.

    Return f_ki, a dict of
    ((far_hall, far_det), (near_hall, near_det)) -> N,
    and N indexed by reco_index.
    """
    f_kji = predict_far_IBD_true_energy(constants, fit_params)
    f_ki = {}
    for (far_halldet, core, near_halldet), n in f_kji.items():
        if (far_halldet, near_halldet) in f_ki:
            f_ki[far_halldet, near_halldet] += n.sum(axis=0)
        else:
            f_ki[far_halldet, near_halldet] = n.sum(axis=0)
    # LBNL comparison requires that we manually correct the AD-to-AD
    # differences in livetime.
    if constants.lbnl_comparison:
        period = period_8ad  # TODO multiple AD periods together
        for (far_halldet, near_halldet), n in f_ki.items():
            far_livetimes = constants.livetime_by_week_AD[far_halldet]
            near_livetimes = constants.livetime_by_week_AD[near_halldet]
            correction = (
                livetime_for_period(far_livetimes, period)
                / livetime_for_period(near_livetimes, period)
            )
            f_ki[far_halldet, near_halldet] *= correction
    return f_ki

def predict_ad_to_ad_obs(constants, fit_params):
    """Compute the number of observed coincidences for a given pair of ADs.

    First the efficiencies are re-applied to convert back to
    "detected" counts, and then the background events are added back in.

    Return f_ki, a dict of
    ((far_hall, far_det), (near_hall, near_det)) -> N,
    and N indexed by reco_index.

    The return value of this function should be directly comparable
    to the number of observed coincidences at the far hall ADs,
    without adjusting the observed counts for backgrounds or any efficiencies.
    """
    f_ki = predict_ad_to_ad_IBDs(constants, fit_params)
    results = f_ki.copy()
    muon_effs = constants.muon_eff
    mult_effs = constants.multiplicity_eff
    dist_time_effs = constants.dist_time_eff
    masses = constants.masses
    for (far_halldet, near_halldet), result in results.items():
        muon_eff = muon_effs[far_halldet]
        mult_eff = mult_effs[far_halldet]
        dist_time_eff = dist_time_effs[far_halldet]
        pull_eff = fit_params.pull_efficiency[far_halldet]
        mass_far = masses[far_halldet]
        mass_eff = mass_far / constants.standard_mass
        pull_rel_escale = fit_params.pull_rel_escale[far_halldet]
        hall, det = far_halldet
        if pull_rel_escale > 0:
            rel_escale_params = constants.rel_escale_parameters[hall][:, 0]
        else:
            rel_escale_params = constants.rel_escale_parameters[hall][:, 1]
        # NB: combined_eff has a bin dependence due to rel_escale_params
        combined_eff = (
            mult_eff * muon_eff * dist_time_eff * mass_eff * (1 + pull_eff)
            * (1 + pull_rel_escale * rel_escale_params)
        )
        results[far_halldet, near_halldet] = (
                result * combined_eff
        )
    pulled_bg = bg_with_pulls(constants, fit_params, halls='far')
    total_bg = ad_dict(0, halls='far')
    for bg_type, bg_dict in pulled_bg.items():
        for halldet, bg_spec in bg_dict.items():
            total_bg[halldet] += bg_spec
    for (far_hall, far_det), near_halldet in f_ki.keys():
        bg_for_far = total_bg[far_hall, far_det]
        results[(far_hall, far_det), near_halldet] += bg_for_far
    return results


def predict_halls(constants, fit_params):
    """Compute the predicted number of IBDs in EH3 based on EH1 or EH2.

    Return f_pred, is a dict with keys
    1 and 2 corresponding to EH1 and EH2, and values of a 1D array
    indexed by reco_index.

    There are 8 AD pairs from a near hall to the far hall,
    and therefore 8 predictions.
    The predictions are all summed to represent combining the far-hall ADs
    and then halved to represent averaging over the 2 near AD predictions.
    """
    f_ki = predict_ad_to_ad_IBDs(constants, fit_params)
    prediction = {
            1: np.zeros_like(reco_bins[:-1]),
            2: np.zeros_like(reco_bins[:-1])
    }
    for (far_halldet, (near_hall, near_det)), ad_prediction in f_ki.items():
        prediction[near_hall] += ad_prediction/2
    return prediction


def generate_bin_averaging_matrix(bins_fine, bins_coarse):
    """Create a matrix that converts from fine binning to coarse binning.

    Weighted average with weight = fine bin width.

    If n = len(bins_coarse) and m = len(bins_fine), returns a matrix
    with dimensions (n-1) x (m-1), since the bin arrays have one extra
    element compared to the number of elements in the histograms they
    describe. Note that the rows of the matrix each sum to unity.

    To convert an array of fine-binned values to the coarser binning,
    multiply by the conversion matrix: ``coarse = np.matmul(M, fine)``.
    The ``matmul`` function will automatically convert 1-D arrays and
    lists to column vectors, as intended.

    Works on the basis of matrix multiplication::

        [ [ 0.5, 0.5, 0.0, 0.0 ],      [ [ 1 ],        [ [ 1.5 ],
          [ 0.0, 0.0, 0.5, 0.5 ] ]  x    [ 2 ],    =     [ 3.5 ] ]
                                         [ 3 ],
                                         [ 4 ] ]

    >>> np.matmul(
    ...     generate_bin_averaging_matrix([0, 1, 2, 3, 4], [0, 2, 4]),
    ...     [1, -1, 3, 4]
    ... )
    array([0, 3.5])


    Assumptions:

    - ``len(values_fine) + 1 == len(bins_fine)``
    - ``bins_fine[0] == bins_coarse[0]``
    - ``bins_fine[-1] == bins_coarse[-1]``
    - ``len(bins_coarse) < len(bins_fine)``
    - bins_coarse is a subset of bins_fine (i.e. no fine bins cross a coarse bin edge)
    - bins_coarse and bins_fine are strictly increasing
    - and that all inputs are 1D arrays.
    """
    bin_weights = np.diff(bins_fine)
    # values_coarse = np.empty((len(bins_coarse)-1,))
    out_matrix = np.zeros((len(bins_coarse) - 1, len(bins_fine) - 1))
    fine_index = 0
    fine_up_edge = bins_fine[0]  # initial value only for first comparison in while loop
    try:
        for coarse_index, coarse_up_edge in enumerate(bins_coarse[1:]):
            row_weights = [0] * fine_index  # pre-fill already-visited fine bins with 0
            # next_coarse_value = 0
            next_coarse_weights_sum = 0
            #while not np.isclose(fine_up_edge, coarse_up_edge):
            while fine_up_edge != coarse_up_edge:
                if fine_up_edge > coarse_up_edge:
                    print(fine_up_edge, coarse_up_edge)
                    raise ValueError('Bin edge mismatch')
                fine_up_edge = bins_fine[fine_index + 1]
                # fine_value = values_fine[fine_index]
                fine_weight = bin_weights[fine_index]
                # next_coarse_value += fine_value * fine_weight
                next_coarse_weights_sum += fine_weight
                row_weights.append(fine_weight)
                fine_index += 1
            # Then we can finalize the bin
            # next_coarse_value /= next_coarse_weights_sum
            row_weights = np.array(row_weights, dtype=float) / next_coarse_weights_sum
            out_matrix[coarse_index, :len(row_weights)] = row_weights
            # values_coarse[coarse_index] = next_coarse_value
    except:
        # Test assumptions
        # test_bins = len(values_fine) + 1 == len(bins_fine)
        test_start_bin = np.isclose(bins_fine[0], bins_coarse[0])
        test_end_bin = np.isclose(bins_fine[-1], bins_coarse[-1])
        test_fine_really_fine = len(bins_coarse) < len(bins_fine)
        test_no_crossings = True  # Trust the user!
        test_increasing = np.all(bin_weights > 0)
        test_1D = True  # Trust the user!
        print(bins_fine[0], bins_coarse[0])
        print(bins_fine[-1], bins_coarse[-1])
        print(f'''start bins match: {test_start_bin}
        end bins match: {test_end_bin}
        fine finer than coarse: {test_fine_really_fine}
        not sure about no bin crossings?
        bins strictly increasing: {test_increasing}
        not sure about accurate 1D arrays''')
        raise
    return out_matrix




def get_baseline(hall, det, core):
    return distances[core][f'EH{hall}'][det-1]

distance_conversion = {
    1: 'D1',
    2: 'D2',
    3: 'L1',
    4: 'L2',
    5: 'L3',
    6: 'L4',
}


distances = {
        1: {
            'EH1': [362.3804, 357.9404],
            'EH2': [1332.4793, 1337.493],
            'EH3': [1919.6319, 1917.5188, 1925.2550, 1923.1489]
            },
        2: {
            'EH1': [371.7626, 368.4142],
            'EH2': [1358.1482, 1362.8763],
            'EH3': [1894.3376, 1891.9765, 1899.8607, 1897.5072]
            },
        3: {
            'EH1': [903.4664, 903.3467],
            'EH2': [467.5740, 472.9709],
            'EH3': [1533.1804, 1534.9194, 1538.9301, 1540.6667]
            },
        4: {
            'EH1': [817.1578, 816.8962],
            'EH2': [489.5774, 495.3458],
            'EH3': [1533.6275, 1535.0322, 1539.4683, 1540.8715]
            },
        5: {
            'EH1': [1353.6177, 1354.2293],
            'EH2': [557.5792, 558.7073],
            'EH3': [1551.3843, 1554.7671, 1556.3439, 1559.7210]
            },
        6: {
            'EH1': [1265.3153, 1265.8861],
            'EH2': [499.2072, 501.0714],
            'EH3': [1524.9400, 1528.0460, 1530.0787, 1533.1791]
            },
        }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config')
    parser.add_argument('-d', '--debug', action='store_true')
    args = parser.parse_args()
    constants = load_constants(args.config)
    fit_params = FitParams(
            0.15,
            2.48e-3,
    )
    print(predict_ad_to_ad_obs(constants, fit_params))
