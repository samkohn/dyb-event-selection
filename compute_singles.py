import argparse
import json
import math
import os
import sqlite3

from scipy.optimize import fsolve

from delayeds import _NH_THU_MAX_TIME


def P_start(R_total, R_mu, Tc):
    """The probability for a coincidence window to start.

    R_total = R_singles + R_correlated

    R_mu = # muon veto windows / muon-subtracted livetime

    Tc = length of coincidence search window, in commensurate units to
    the rates (Hz + s, MHz + us, etc.)
    """
    R_sum = R_total + R_mu
    sum_exp = math.exp(-R_sum * Tc)
    mu_exp = math.exp(-R_mu * Tc)
    term_1 = R_mu / R_sum * (1 - sum_exp) + sum_exp
    term_2 = R_total / R_sum * mu_exp * (1 - sum_exp)
    term_3 = (
        -R_total
        / (2 * R_total + R_mu)
        * mu_exp
        * (1 - math.exp(-(2 * R_total + R_mu) * Tc))
    )
    return term_1 + term_2 + term_3


def poisson(n, mu):
    """The Poisson probability for a count of n and mean mu."""
    return mu ** n * math.exp(-mu) / math.factorial(n)


def single_uncorr_rate(R_single, R_corr, R_mu, Tc):
    """The rate of single events (multiplicity=1) due to uncorrelateds."""
    opportunity_rate = R_single
    prob_of_starting = P_start(R_single + R_corr, R_mu, Tc)
    prob_of_finishing_alone = poisson(0, (R_single + R_corr) * Tc)
    return opportunity_rate * prob_of_starting * prob_of_finishing_alone

def multiplicity_efficiency(R_single, R_corr, R_mu, Tc):
    """The multiplicity veto efficiency.

    Probability of an event at a given time not lying within a previous
    coincidence window and having no uncorrelated events within
    its own coincidence window.
    """
    prob_of_starting = P_start(R_single + R_corr, R_mu, Tc)
    prob_of_finishing_alone = poisson(0, (R_single + R_corr) * Tc)
    return prob_of_starting * prob_of_finishing_alone

def single_e_rate(R_single, R_corr, R_mu, Tc, neutron_efficiency):
    """The rate of single events due to IBD/correlated positrons."""
    if R_corr == 0:
        return 0
    opportunity_rate = R_corr
    prob_of_starting = P_start(R_single + R_corr, R_mu, Tc)
    prob_of_finishing_alone = poisson(0, (R_single + R_corr) * Tc) * (
        1 - neutron_efficiency
    )
    return opportunity_rate * prob_of_starting * prob_of_finishing_alone


def single_n_rate(
    R_single, R_corr, R_mu, Tc, neutron_efficiency, tau_Gd, tau_LS, alpha
):
    """The rate of single events due to IBD/correlated neutrons.

    tau_Gd = capture time on Gd
    tau_LS = capture time on H (LS)
    alph = weighting of Gd vs. LS.
    """  # TODO define alpha
    if R_corr == 0:
        return 0
    exp = math.exp
    lambda_Gd = 1 / tau_Gd
    lambda_LS = 1 / tau_LS
    term_1 = (
        neutron_efficiency
        * R_mu
        * (
            alpha * subterm(R_single + R_mu + lambda_Gd, Tc)
            + (1 - alpha) * subterm(R_single + R_mu + lambda_LS, Tc)
        )
    )
    term_2 = (
        neutron_efficiency
        * exp(-(lambda_Gd + lambda_LS + R_mu + R_single) * Tc)
        * (alpha * exp(Tc / tau_LS) + (1 - alpha) * exp(Tc / tau_Gd))
    )
    term_3 = (
        neutron_efficiency
        * exp(-R_mu * Tc)
        * (
            alpha * lambda_Gd * subterm(R_mu + lambda_Gd, Tc)
            + (1 - alpha) * lambda_LS * subterm(R_mu + lambda_LS, Tc)
            - alpha * lambda_Gd * subterm(R_mu + R_single + lambda_Gd, Tc)
            - (1 - alpha) * lambda_LS * subterm(R_mu + R_single + lambda_LS, Tc)
        )
    )
    return R_corr * (term_1 + term_2 + term_3) * poisson(0, R_single * Tc)


def single_rate(R_single, R_corr, R_mu, Tc, neutron_efficiency, tau_Gd, tau_LS, alpha):
    """The total single (multiplicity-1) event rate."""
    return (
        single_uncorr_rate(R_single, R_corr, R_mu, Tc)
        + single_e_rate(R_single, R_corr, R_mu, Tc, neutron_efficiency)
        + single_n_rate(
            R_single, R_corr, R_mu, Tc, neutron_efficiency, tau_Gd, tau_LS, alpha
        )
    )


def subterm(sum_term, Tc):
    """Subterm 1 of the complicated neutron rate equation."""
    return (1 - math.exp(-sum_term * Tc)) / sum_term


def main(infile, database, update_db, iteration, extra_cut):
    with open(os.path.splitext(infile)[0] + '.json', 'r') as f:
        stats = json.load(f)
    livetime_s = stats['usable_livetime']/1e9
    ad = stats['ad']
    import ROOT
    ch = ROOT.TChain('ad_events')
    ch.Add(infile)
    ch.GetEntry(0)
    runNo = ch.run
    site = ch.site
    start_time = ch.timestamp[0]
    multiplicity_1_count = ch.Draw('energy', f'detector == {ad} && '
        f'multiplicity == 1 && ({extra_cut})', 'goff')
    multiplicity_1_count_error = math.sqrt(multiplicity_1_count)
    multiplicity_1_rate_Hz = multiplicity_1_count / livetime_s
    with sqlite3.Connection(database) as conn:
        cursor = conn.cursor()
        cursor.execute('''SELECT Rate_Hz FROM muon_rates WHERE
            RunNo = ? AND DetNo = ?''', (runNo, ad))
        muon_rate, = cursor.fetchone()

    # convert to seconds and subtract off 1us
    window_size = _NH_THU_MAX_TIME/1e9 - 1e-6
    if iteration > 0:
        raise NotImplementedError("Haven't connected to database")
    else:
        R_corr = 0
        neutron_efficiency = None
        tau_Gd = None
        tau_LS = None
        alpha = None

    parameters = (R_corr, muon_rate, window_size, neutron_efficiency, tau_Gd,
            tau_LS, alpha)
    underlying_uncorr_rate = fsolve(lambda x: single_rate(x, *parameters) -
            multiplicity_1_rate_Hz, multiplicity_1_rate_Hz)[0]
    multiplicity_eff = multiplicity_efficiency(underlying_uncorr_rate,
            R_corr, muon_rate, window_size)
    if update_db:
        with sqlite3.Connection(database) as conn:
            cursor = conn.cursor()
            cursor.execute('''INSERT OR REPLACE INTO singles_rates
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)''', (
                runNo,
                ad,
                iteration,
                underlying_uncorr_rate,
                multiplicity_1_count_error/livetime_s,
                multiplicity_1_count,
                multiplicity_1_rate_Hz,
                R_corr,
                multiplicity_eff,
            ))
    else:
        print(f'multiplicity-1 rate: {multiplicity_1_rate_Hz} Hz')
        print(f'relative error: {100/multiplicity_1_count_error:.2f}%')
        print(f'underlying uncorr. rate: {underlying_uncorr_rate} Hz')
        print(f'multiplicity efficiency: {multiplicity_eff}')
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('database')
    parser.add_argument('--update-db', action='store_true',
            help='If present, store the results in the given db file')
    parser.add_argument('--iteration', type=int, default=0,
            help='If present, look up IBD rate in database and integrate into '
                'the calculation')
    parser.add_argument('--extra-cut', default='1')
    args = parser.parse_args()
    main(args.infile, args.database, args.update_db, args.iteration, args.extra_cut)
