"""A script that uses the toymc package."""

import argparse

from toymc.single import Single
from toymc.correlated import Correlated
from toymc.muon import Muon
from toymc import ToyMC


def main(outfile, runtime, seed):
    """Run the ToyMC with the given configuration."""
    toymc = ToyMC(outfile, runtime, seed=seed)
    # Single(name, rate_Hz, EH, AD)
    single = Single("Single_event", 20, 1, 1)
    single.truth_label = 0
    single.energy_spectrum = lambda rng: rng.uniform(1.5, 3)
    # Correlated(name, EH, AD, rate_Hz, coincidence_time_ns)
    ibd_nGd = Correlated("IBD_nGd", 1, 1, 0.007, 28000)
    ibd_nGd.truth_label_prompt = 1
    ibd_nGd.truth_label_delayed = 2
    ibd_nH = Correlated("IBD_nH", 1, 1, 0.006, 150000)
    ibd_nH.truth_label_prompt = 3
    ibd_nH.truth_label_delayed = 4
    ibd_nH.delayed_energy_spectrum = lambda rng: rng.uniform(1.9, 2.3)
    ibd_nH.prompt_delayed_distance_mm = 100
    # Muon(name, EH, rate_Hz)
    # Muon events include correlated WP, AD, and Shower muons with
    # configurable rate ratios
    muon = Muon("Muon", 1, 200)
    muon.truth_label_WP = 5
    muon.truth_label_AD = 6
    muon.truth_label_shower = 7
    toymc.add_event_type(single)
    toymc.add_event_type(ibd_nGd)
    toymc.add_event_type(ibd_nH)
    toymc.add_event_type(muon)
    toymc.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Daya Bay Toy MC by Sam Kohn")
    parser.add_argument("outfile")
    parser.add_argument("-t", "--runtime", type=int, help="DAQ runtime in seconds")
    parser.add_argument("-s", "--seed", default=None, type=int, help="random seed")
    args = parser.parse_args()
    main(args.outfile, args.runtime, args.seed)
