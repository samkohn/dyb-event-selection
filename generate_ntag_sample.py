"""A script that uses the toymc package."""

import argparse

from toymc.single import Single
from toymc.correlated import Correlated
from toymc.muon import Muon
from toymc import ToyMC, Event

class MuonWithLi9(Muon):
    def __init__(self, name, site, rate_Hz):
        super().__init__(name, site, rate_Hz)
        self.truth_label_ntag_AD = None
        self.truth_label_ntag_shower = None
        self.truth_label_li9_AD_ntag_prompt = None
        self.truth_label_li9_AD_no_ntag_prompt = None
        self.truth_label_li9_shower_ntag_prompt = None
        self.truth_label_li9_shower_no_ntag_prompt = None
        self.truth_label_li9_AD_ntag_delayed = None
        self.truth_label_li9_AD_no_ntag_delayed = None
        self.truth_label_li9_shower_ntag_delayed = None
        self.truth_label_li9_shower_no_ntag_delayed = None
        self.prob_ntag_AD = 0.001
        self.prob_ntag_shower = 0.5
        self.prob_li9_no_ntag = 0.001
        self.prob_li9_ntag = 1.0
        self.li9_lifetime = 0.25723
        self.ntag_energy_spectrum = lambda rng: rng.uniform(1.8, 12)
        self.ntag_position_spectrum = lambda rng: (0, 0, 0)
        self.ntag_time_delay_ns = lambda rng: rng.integers(20000, 200000)
        self.ADMuon_energy_spectrum = lambda rng: rng.uniform(1000, 2000)
        self.li9_prompt_spectrum = lambda rng: rng.uniform(6, 8)
        self.li9_delayed_spectrum = lambda rng: 2.2
        self.li9_prompt_position = lambda rng: (0, 0, 0)
        self.li9_delayed_pos_from_prompt = lambda rng, pos: (1, 1, 1)

        self.prob_WP_and_AD = 0.1995
        self.prob_WP_and_shower = 0.0005

    def labels(self):
        superlabels = super().labels()
        superlabels[self.truth_label_ntag_AD] = self.name + '_neutrontag_ADmuon'
        superlabels[self.truth_label_ntag_shower] = self.name + '_neutrontag_showermuon'
        superlabels[self.truth_label_li9_AD_ntag_prompt] = \
                'Li9_prompt_from_ADmuon_ntag'
        superlabels[self.truth_label_li9_AD_no_ntag_prompt] = \
                'Li9_prompt_from_ADmuon_no_ntag'
        superlabels[self.truth_label_li9_shower_ntag_prompt] = \
                'Li9_prompt_from_showermuon_ntag'
        superlabels[self.truth_label_li9_shower_no_ntag_prompt] = \
                'Li9_prompt_from_showermuon_no_ntag'
        superlabels[self.truth_label_li9_AD_ntag_delayed] = \
                'Li9_delayed_from_ADmuon_ntag'
        superlabels[self.truth_label_li9_AD_no_ntag_delayed] = \
                'Li9_delayed_from_ADmuon_no_ntag'
        superlabels[self.truth_label_li9_shower_ntag_delayed] = \
                'Li9_delayed_from_showermuon_ntag'
        superlabels[self.truth_label_li9_shower_no_ntag_delayed] = \
                'Li9_delayed_from_showermuon_no_ntag'
        return superlabels

    def generate_events(self, rng, duration_s, t0_s):
        superevents = super().generate_events(rng, duration_s, t0_s)
        actual_number = self.actual_event_count(rng, duration_s, self.rate_hz)
        number_admuons = int(actual_number * self.prob_WP_and_AD)
        number_showermuons = int(actual_number * self.prob_WP_and_shower)
        # The first actual_number events are WP muons. The next
        # number_admuons are AD muons, and the last number_showermuons
        # are shower muons.
        number_ntag_AD = int(self.prob_ntag_AD * number_admuons)
        number_ntag_shower = int(self.prob_ntag_shower * number_showermuons)
        #Li9
        number_li9_AD_ntag = int(self.prob_li9_ntag * number_ntag_AD)
        number_li9_AD_no_ntag = int(self.prob_li9_no_ntag * (number_admuons -
            number_ntag_AD))
        number_li9_shower_ntag = int(self.prob_li9_ntag * number_ntag_shower)
        number_li9_shower_no_ntag = int(self.prob_li9_no_ntag * (number_showermuons -
            number_ntag_shower))

        print(number_ntag_AD)
        print(number_ntag_shower)
        # The list of events is broken down into sections based on ntag and li9:
        # [ WP Muons | AD ntag li9 | AD ntag | AD no_ntag Li9 | AD | ...
        #   Shower ntag li9 | shower ntag | shower no_ntag Li9 | shower ]
        wp_end = actual_number
        ad_ntag_li9_end = wp_end + number_li9_AD_ntag
        ad_ntag_end = actual_number + number_ntag_AD
        ad_no_ntag_li9_end = ad_ntag_end + number_li9_AD_no_ntag
        ad_end = actual_number + number_admuons
        shower_ntag_li9_end = ad_end + number_li9_shower_ntag
        shower_ntag_end = ad_end + number_ntag_shower
        shower_no_ntag_li9_end = shower_ntag_end + number_li9_shower_no_ntag
        # Lists of the base muon events for each category
        ad_ntag_li9_muons = superevents[wp_end:ad_ntag_li9_end]
        ad_ntag_no_li9_muons = superevents[ad_ntag_li9_end:ad_ntag_end]
        ad_no_ntag_li9_muons = superevents[ad_ntag_end:ad_no_ntag_li9_end]
        shower_ntag_li9_muons = superevents[ad_end:shower_ntag_li9_end]
        shower_ntag_no_li9_muons = superevents[shower_ntag_li9_end:shower_ntag_end]
        shower_no_ntag_li9_muons = superevents[shower_ntag_end:shower_no_ntag_li9_end]
        print(len(ad_ntag_li9_muons))
        print(len(ad_ntag_no_li9_muons))
        print(len(shower_ntag_li9_muons))
        print(len(shower_ntag_no_li9_muons))

        new_events = []
        for event in ad_ntag_li9_muons:
            ntag_event = self.new_ntag_event(rng, event.timestamp, event.detector,
                    self.truth_label_ntag_AD)
            li9_events = self.new_li9_pair(rng, event.timestamp, event.detector,
                    self.truth_label_li9_AD_ntag_prompt,
                    self.truth_label_li9_AD_ntag_delayed)
            new_events.append(ntag_event)
            new_events.extend(li9_events)
        for event in ad_ntag_no_li9_muons:
            ntag_event = self.new_ntag_event(rng, event.timestamp, event.detector,
                    self.truth_label_ntag_AD)
            new_events.append(ntag_event)
        for event in ad_no_ntag_li9_muons:
            li9_events = self.new_li9_pair(rng, event.timestamp, event.detector,
                    self.truth_label_li9_AD_no_ntag_prompt,
                    self.truth_label_li9_AD_no_ntag_delayed)
            new_events.extend(li9_events)
        for event in shower_ntag_li9_muons:
            ntag_event = self.new_ntag_event(rng, event.timestamp, event.detector,
                    self.truth_label_ntag_shower)
            li9_events = self.new_li9_pair(rng, event.timestamp, event.detector,
                    self.truth_label_li9_shower_ntag_prompt,
                    self.truth_label_li9_shower_ntag_delayed)
            new_events.append(ntag_event)
            new_events.extend(li9_events)
        for event in shower_ntag_no_li9_muons:
            ntag_event = self.new_ntag_event(rng, event.timestamp, event.detector,
                    self.truth_label_ntag_shower)
            new_events.append(ntag_event)
        for event in shower_no_ntag_li9_muons:
            li9_events = self.new_li9_pair(rng, event.timestamp, event.detector,
                    self.truth_label_li9_shower_no_ntag_prompt,
                    self.truth_label_li9_shower_no_ntag_delayed)
            new_events.extend(li9_events)
        return superevents + new_events

    def new_li9_pair(self, rng, muon_timestamp, muon_AD, truth_label_prompt,
            truth_label_delayed):
        prompt_position = self.li9_prompt_position(rng)
        delayed_position = self.li9_delayed_pos_from_prompt(rng, prompt_position)
        prompt_timestamp = int(muon_timestamp + rng.exponential(1e9*self.li9_lifetime))
        delayed_timestamp = prompt_timestamp + 200000
        prompt = Event(
            truth_label_prompt,
            1,
            prompt_timestamp,
            muon_AD,
            self.trigger_type,
            self.site,
            self.li9_prompt_spectrum(rng),
            192,
            500,
            prompt_position[0],
            prompt_position[1],
            prompt_position[2],
            0.1,
            0.1,
            0.99,
            0.99,
            0
        )
        delayed = Event(
            truth_label_delayed,
            1,
            delayed_timestamp,
            muon_AD,
            self.trigger_type,
            self.site,
            self.li9_delayed_spectrum(rng),
            192,
            500,
            delayed_position[0],
            delayed_position[1],
            delayed_position[2],
            0.1,
            0.1,
            0.99,
            0.99,
            0
        )
        return (prompt, delayed)

    def new_ntag_event(self, rng, muon_timestamp, muon_AD, truth_label):
        position = self.ntag_position_spectrum(rng)
        event = Event(
            truth_label,
            1,
            muon_timestamp + self.ntag_time_delay_ns(rng),
            muon_AD,
            self.trigger_type,
            self.site,
            self.ntag_energy_spectrum(rng),
            192,
            500,
            position[0],
            position[1],
            position[2],
            0.1,
            0.1,
            0.99,
            0.99,
            0
        )
        return event


def main(outfile, runtime, t0, seed):
    """Run the ToyMC with the given configuration."""
    toymc = ToyMC(outfile, runtime, t0, seed=seed)
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
    muon = MuonWithLi9("Muon", 1, 200)
    muon.avail_ads = (1,)
    muon.truth_label_WP = 5
    muon.truth_label_AD = 6
    muon.truth_label_shower = 7
    muon.truth_label_ntag_AD = 8
    muon.truth_label_ntag_shower = 9
    muon.truth_label_li9_AD_ntag_prompt = 10
    muon.truth_label_li9_AD_no_ntag_prompt = 11
    muon.truth_label_li9_shower_ntag_prompt = 12
    muon.truth_label_li9_shower_no_ntag_prompt = 13
    muon.truth_label_li9_AD_ntag_delayed = 14
    muon.truth_label_li9_AD_no_ntag_delayed = 15
    muon.truth_label_li9_shower_ntag_delayed = 16
    muon.truth_label_li9_shower_no_ntag_delayed = 17
    toymc.add_event_type(single)
    toymc.add_event_type(ibd_nGd)
    toymc.add_event_type(ibd_nH)
    toymc.add_event_type(muon)
    toymc.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Daya Bay Toy MC by Sam Kohn")
    parser.add_argument("outfile")
    parser.add_argument("-t", "--runtime", type=int, help="DAQ runtime in seconds")
    parser.add_argument("--t0", "--start-time", type=int, help="Start time of run")
    parser.add_argument("-s", "--seed", default=None, type=int, help="random seed")
    args = parser.parse_args()
    main(args.outfile, args.runtime, args.t0, args.seed)
