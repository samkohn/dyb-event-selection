"""Dump the standard summary table as a tex file."""

import argparse
from dataclasses import dataclass, field as dc_field
import json

import numpy as np

import common
import compute_dataset_summary
import compute_acc_summary
import prediction


@dataclass
class SummaryTable:
    daq_livetime: list
    muon_eff: list
    multiplicity_eff: list
    muon_rate: list
    singles_rate: list
    num_observed_coincidences: list
    num_observed_coincidences_err: list
    num_accidentals: list
    num_accidentals_err: list
    num_correlated: list
    num_correlated_err: list
    rate_accidentals: list
    rate_accidentals_err: list
    rate_li9: list
    rate_li9_err: list
    rate_fast_neutron: list
    rate_fast_neutron_err: list
    rate_amc: list
    rate_amc_err: list
    rate_ibd: list
    rate_ibd_err: list
    relative_num_target_protons: list
    num_target_protons_err: float
    prompt_energy_efficiency_err: float
    delayed_energy_efficiency_err: float
    distance_time_cut_efficiency_err: float
    combined_efficiency_err: float
    combined_ad_uncorr_err: float
    rel_escale_uncertainty: float

    @classmethod
    def from_config(cls, config_file_name):
        constants = prediction.load_constants(config_file_name)
        with open(config_file_name, 'r') as f:
            config_dict = json.load(f)
            config = prediction.Config(**config_dict)
        values = {}
        values['daq_livetime'] = to_list(
            {
                halldet: livetime/60/60/24
                for halldet, livetime in constants.daq_livetimes.items()
            }
        )
        values['muon_eff'] = to_list(constants.muon_eff)
        muon_database = config.muon_eff
        if not isinstance(muon_database, str):
            muon_database = config.database
        muon_label = config.muon_eff_label
        if not isinstance(muon_label, str):
            raise ValueError("Cannot compute muon rate: missing muon_eff_label")
        values['muon_rate'] = compute_dataset_summary.muon_rate_Hz(
            muon_database, muon_label
        ).tolist()
        values['multiplicity_eff'] = to_list(constants.multiplicity_eff)
        mult_database = config.mult_eff
        if not isinstance(mult_database, str):
            mult_database = config.database
        mult_label = config.mult_eff_label
        if not isinstance(mult_label, str):
            raise ValueError("Cannot compute singles rate: missing mult_eff_label")
        values['singles_rate'] = compute_dataset_summary.singles_rate_Hz(
            mult_database, mult_label
        ).tolist()
        values['num_observed_coincidences'] = to_list(
            {halldet: sum(obs) for halldet, obs in constants.observed_candidates.items()}
        )
        values['num_observed_coincidences_err'] = to_list(
            {
                halldet: np.sqrt(sum(obs))
                for halldet, obs in constants.observed_candidates.items()
            }
        )
        values['num_accidentals'] = to_list(
            {
                halldet: sum(acc)
                for halldet, acc in constants.nominal_bgs['accidental'].items()
            }
        )
        values['num_accidentals_err'] = to_list(
            {
                halldet: err * acc
                for (halldet, err), acc in zip(
                    constants.bg_errors['accidental'].items(),
                    values['num_accidentals'],
                )
            }
        )
        values['num_correlated'] = [
            obs - acc
            for obs, acc in zip(
                values['num_observed_coincidences'], values['num_accidentals']
            )
        ]
        values['num_correlated_err'] = [
            np.sqrt(obs_err * obs_err + acc_err * acc_err)
            for obs_err, acc_err in zip(
                values['num_observed_coincidences_err'], values['num_accidentals_err']
            )
        ]
        # Effective livetime = DAQ livetime * mu_eff * mult_eff
        effective_livetime_days = (
            np.array(to_list(constants.daq_livetimes))/60/60/24
            * to_list(constants.muon_eff)
            * to_list(constants.multiplicity_eff)
        )
        lookup = {
            'accidental': 'accidentals',
            'li9': 'li9',
            'fast-neutron': 'fast_neutron',
            'amc': 'amc',
        }
        for bg_name in ('accidental', 'li9', 'fast-neutron', 'amc'):
            counts = to_list(
                {
                    halldet: sum(bg)
                    for halldet, bg in constants.nominal_bgs[bg_name].items()
                }
            )
            err = to_list(
                {
                    halldet: err * bg
                    for (halldet, err), bg in zip(
                        constants.bg_errors[bg_name].items(), counts
                    )
                }
            )
            rate = counts/effective_livetime_days
            rate_err = err/effective_livetime_days
            values[f'rate_{lookup[bg_name]}'] = rate.tolist()
            values[f'rate_{lookup[bg_name]}_err'] = rate_err.tolist()
        num_ibd = (
            np.array(values['num_observed_coincidences'])
            - values['num_accidentals']
            - values['rate_li9'] * effective_livetime_days
            - values['rate_fast_neutron'] * effective_livetime_days
            - values['rate_amc'] * effective_livetime_days
        )
        values['rate_ibd'] = (num_ibd/effective_livetime_days).tolist()
        num_ibd_err = np.sqrt(
            np.array(values['num_observed_coincidences_err'])**2
            + np.array(values['num_accidentals_err'])**2
            + np.array(values['rate_li9_err'] * effective_livetime_days)**2
            + np.array(values['rate_fast_neutron_err'] * effective_livetime_days)**2
            + np.array(values['rate_amc_err'] * effective_livetime_days)**2
        )
        values['rate_ibd_err'] = (num_ibd_err/effective_livetime_days).tolist()
        values['relative_num_target_protons'] = (
            np.array(to_list(constants.masses))/constants.standard_mass - 1
        ).tolist()
        values['num_target_protons_err'] = 0.0037
        values['prompt_energy_efficiency_err'] = 0.001
        values['delayed_energy_efficiency_err'] = 0.002
        values['distance_time_cut_efficiency_err'] = 0.005
        values['combined_efficiency_err'] = np.sqrt(
            values['prompt_energy_efficiency_err']**2
            + values['delayed_energy_efficiency_err']**2
            + values['distance_time_cut_efficiency_err']**2
        )
        values['combined_ad_uncorr_err'] = np.sqrt(
            values['combined_efficiency_err']**2
            + values['num_target_protons_err']**2
        )
        values['rel_escale_uncertainty'] = 0.005
        return cls(**values)

    def to_simple_python(self):
        """Return a dict of nested lists representing various tables."""
        to_return = {}
        rate_table = []
        rate_table.append(self.daq_livetime)
        rate_table.append(self.muon_eff)
        rate_table.append(self.multiplicity_eff)
        rate_table.append(self.muon_rate)
        rate_table.append(self.singles_rate)
        rate_table.append(self.num_observed_coincidences)
        rate_table.append(self.num_observed_coincidences_err)
        rate_table.append(self.num_accidentals)
        rate_table.append(self.num_accidentals_err)
        rate_table.append(self.num_correlated)
        rate_table.append(self.num_correlated_err)
        rate_table.append(self.rate_accidentals)
        rate_table.append(self.rate_accidentals_err)
        rate_table.append(self.rate_li9)
        rate_table.append(self.rate_li9_err)
        rate_table.append(self.rate_fast_neutron)
        rate_table.append(self.rate_fast_neutron_err)
        rate_table.append(self.rate_amc)
        rate_table.append(self.rate_amc_err)
        rate_table.append(self.rate_ibd)
        rate_table.append(self.rate_ibd_err)
        rate_table.append(self.relative_num_target_protons)
        rate_table.append([self.num_target_protons_err] * 8)
        to_return['rates'] = rate_table
        efficiency_table = []
        efficiency_table.append(self.prompt_energy_efficiency_err)
        efficiency_table.append(self.delayed_energy_efficiency_err)
        efficiency_table.append(self.distance_time_cut_efficiency_err)
        efficiency_table.append(self.combined_efficiency_err)
        efficiency_table.append(self.rel_escale_uncertainty)
        efficiency_table.append(self.combined_ad_uncorr_err)
        to_return['efficiencies'] = efficiency_table
        return to_return

    def to_json(self):
        return json.dumps(self.__dict__, indent=4)

    def to_latex(self):
        def create_ad_columns(name, str_vals, prefix=' '*8):
            line = '&'.join(str_vals)
            return prefix + name + '& ' + line + r' \\' + '\n\\hline\n'
        def add_errors(val_strs, err_strs):
            return [
                '$' + val_str + r' \pm ' + err_str + '$'
                for val_str, err_str in zip(val_strs, err_strs)
            ]
        def map_str_n_decimals(vals, n):
            return map(lambda val: f'{val:.{n}f}', vals)
        to_return_parts = [
            r'''\begin{table}[ht]
    \caption{Caption goes here}
    \label{tab:label_goes_here}
    \centering
    \tiny
    \setlength{\tabcolsep}{2.5pt}
    \sisetup{
        per-mode=reciprocal
    }
    \begin{tabular}[t]{rrrrrrrrr}  % 1 label and 8 ADs
        \hline
        & \multicolumn{2}{c}{EH1}
        & \multicolumn{2}{c}{EH2}
        & \multicolumn{4}{c}{EH3} \\''' + '\n'
        ]
        ad_names = [f'EH{hall}-AD{det}' for hall, det in prediction.all_ads]
        to_return_parts.append(create_ad_columns('', ad_names))
        to_return_parts.append(r'\arrayrulecolor{lightgray}' '\n')
        to_return_parts.append(create_ad_columns(
            r'$T_{\text{DAQ}}$ [\si{\day}]',
            map_str_n_decimals(self.daq_livetime, 3),
        ))
        to_return_parts.append(create_ad_columns(
            r'$\Delta N_{\text{p}}$ [\%]',
            map_str_n_decimals(self.relative_num_target_protons * 100, 2),
        ))
        to_return_parts.append(create_ad_columns(
            r'$\varepsilon_{\mu}$',
            map_str_n_decimals(self.muon_eff, 4),
        ))
        to_return_parts.append(create_ad_columns(
            r'$\varepsilon_{m}$',
            map_str_n_decimals(self.multiplicity_eff, 4),
        ))
        to_return_parts.append(create_ad_columns(
            r'$R_{\mu}$ [\si{\Hz}]',
            map_str_n_decimals(self.muon_rate, 2),
        ))
        to_return_parts.append(create_ad_columns(
            r'$R_s$ [\si{\Hz}]',
            map_str_n_decimals(self.singles_rate, 3),
        ))
        to_return_parts.append(create_ad_columns(
            r'$N_{\text{cand}}$',
            map(str, self.num_observed_coincidences),
        ))
        to_return_parts.append(create_ad_columns(
            r'$N_{\text{acc}}$',
            add_errors(
                map_str_n_decimals(self.num_accidentals, 0),
                map_str_n_decimals(self.num_accidentals_err, 0),
            ),
        ))
        to_return_parts.append(create_ad_columns(
            r'$N_{\text{corr}}$',
            add_errors(
                map_str_n_decimals(self.num_correlated, 0),
                map_str_n_decimals(self.num_correlated_err, 0),
            ),
        ))
        to_return_parts.append(create_ad_columns(
            r'$R_{\text{acc}}$ [\si{\per\day}]',
            add_errors(
                map_str_n_decimals(self.rate_accidentals, 2),
                map_str_n_decimals(self.rate_accidentals_err, 2),
            ),
        ))
        to_return_parts.append(create_ad_columns(
            r'$R_{\text{Li9}}$ [\si{\per\day}]',
            add_errors(
                map_str_n_decimals(self.rate_li9, 2),
                map_str_n_decimals(self.rate_li9_err, 2),
            ),
        ))
        to_return_parts.append(create_ad_columns(
            r'$R_{\text{FastN}}$ [\si{\per\day}]',
            add_errors(
                map_str_n_decimals(self.rate_fast_neutron, 2),
                map_str_n_decimals(self.rate_fast_neutron_err, 2),
            ),
        ))
        to_return_parts.append(r'\arrayrulecolor{black}' '\n')
        to_return_parts.append(create_ad_columns(
            r'$R_{\text{AmC}}$ [\si{\per\day}]',
            add_errors(
                map_str_n_decimals(self.rate_amc, 2),
                map_str_n_decimals(self.rate_amc_err, 2),
            ),
        ))
        to_return_parts.append(create_ad_columns(
            r'$R_{\text{IBD}}$ [\si{\per\day}]',
            add_errors(
                map_str_n_decimals(self.rate_ibd, 2),
                map_str_n_decimals(self.rate_ibd_err, 2),
            ),
        ))
        to_return_parts.append('    \\end{tabular}\n\\end{table}\n')

        # Efficiencies table
        def efficiency_row(name, uncertainty):
            label = '    ' + name + ' & '
            unc = f'{uncertainty*100:.2f}'
            return label + unc + r'\\' + '\n'
        to_return_parts.append(
            r'''\begin{table}[ht]
    \caption{Caption goes here}
    \label{tab:label_goes_here2}
    \centering
    \begin{tabular}[t]{rr}
        \hline
        & Uncertainty (\%) \\
        \hline''' + '\n'
        )
        to_return_parts.append(efficiency_row(
            'Relative energy scale uncertainty',
            self.rel_escale_uncertainty,
        ))
        to_return_parts.append(efficiency_row(
            'Prompt energy cut eff.',
            self.prompt_energy_efficiency_err,
        ))
        to_return_parts.append(efficiency_row(
            'Delayed energy cut eff.',
            self.delayed_energy_efficiency_err,
        ))
        to_return_parts.append(efficiency_row(
            'Distance-time (DT) cut eff.',
            self.distance_time_cut_efficiency_err,
        ))
        to_return_parts.append(r'\hline' + '\n')
        to_return_parts.append(efficiency_row(
            'Combined cut eff.',
            self.combined_efficiency_err,
        ))
        to_return_parts.append(efficiency_row(
            'Target proton number',
            self.num_target_protons_err,
        ))
        to_return_parts.append(r'\hline' + '\n')
        to_return_parts.append(efficiency_row(
            'Total AD-uncorrelated uncertainty',
            self.combined_ad_uncorr_err,
        ))
        to_return_parts.append('    \\end{tabular}\n\\end{table}')
        return ''.join(to_return_parts)






def main(database, rates_label, label, coinc_label, output):
    """Rows:

    - T_daq [days]
    - eps_mu
    - eps_mult
    - R_mu [Hz]
    - R_s [Hz]
    - number observed coincidences
    - number accidentals
    - number correlated events (acc-sub)
    - R_acc [/day]
    - R_Li9 [/day]
    - R_FN [/day]
    - R_AmC [/day]
    - R_IBD [/day] (all subtracted)

    Efficiency & identicalness rows:
    - N_p
    - prompt cut efficiency
    - Distance-Time (DT)
    - delayed cut
    - combined
    """
    daq_livetime_days = compute_dataset_summary.daq_livetime_days(database, rates_label)
    eps_mu = compute_dataset_summary.muon_efficiency(database, rates_label)
    eps_mult = compute_dataset_summary.multiplicity_efficiency(database, rates_label)
    muon_rate_Hz = compute_dataset_summary.muon_rate_Hz(database, rates_label)
    singles_rate_Hz = compute_dataset_summary.singles_rate_Hz(database, rates_label)
    num_coincidences = compute_dataset_summary.coincidences_counts(database, coinc_label)
    num_coincidences_error = np.sqrt(num_coincidences)
    num_acc = compute_acc_summary.counts(database, label, rates_label)
    num_acc_error = compute_acc_summary.count_errors(database, label, rates_label)
    num_correlateds = num_coincidences - num_acc
    num_correlateds_error = np.hypot(num_coincidences_error, num_acc_error)





def to_list(ad_dict):
    """Return the values of the dict in standard AD order (EH1-AD1 -> EH3-AD4)."""
    to_return = []
    for halldet in prediction.all_ads:
        try:
            to_return.append(ad_dict[halldet].item())
        except AttributeError:
            to_return.append(ad_dict[halldet])
    return to_return




def main(config, output, use_json):
    table = SummaryTable.from_config(config)
    if use_json:
        result = table.to_json()
    else:
        result = table.to_latex()
    if output is None:
        print(result)
    else:
        with open(output, 'w') as fout:
            fout.write(result)
    return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Dump a summary table to LaTeX format'
    )
    parser.add_argument('config')
    parser.add_argument('-o', '--output')
    parser.add_argument('--json', action='store_true', help='Use JSON format')
    args = parser.parse_args()
    main(args.config, args.output, args.json)
