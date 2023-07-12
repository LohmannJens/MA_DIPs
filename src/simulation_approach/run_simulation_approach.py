'''

'''
import os
import sys

import numpy as np
import matplotlib.pyplot as plt


sys.path.insert(0, "..")
from utils import RESULTSPATH


class DIPClass():
    def __init__(self, name, n_init, prob, r_increase, r_decrease, de_novo_event_prob):
        self.name = name
        self.n = n_init
        self.p = prob
        self.r_incr = r_increase
        self.r_decr = r_decrease
        self.de_novo_event_prob = de_novo_event_prob
        self.candidate_counts = np.array([np.random.randint(1, 20) for _ in range(self.n)])

    def get_size(self):
        return len(self.candidate_counts)

    def get_count(self):
        return np.sum(self.candidate_counts)

class DIPPopulationSimulator():
    def __init__(self, classes, n_max, cand_max):
        self.classes = classes
        self.n_max = n_max
        self.cand_max = cand_max
        self.timeseries = dict({c.name : list() for c in self.classes})
    
    def get_population_size(self):
        return np.sum([c.get_size() for c in self.classes])
    
    def get_population_count(self):
        return np.sum([c.get_count() for c in self.classes])

    def run_simulation(self, steps):
        for i in range(steps):
            loss_event_p = 0.1
            loss_event = True if loss_event_p > np.random.random() else False

            for cls in self.classes:
                new_counts = list()
                for cand in cls.candidate_counts:
                    thresh = np.random.random()

                    p = 0.1 if loss_event else cls.p

                    if p > thresh:
                        change = round(cand * (1 - self.get_population_count() / self.n_max) * cls.r_incr )#* 0.2) # here I have an extra artificial factor included right now
                    else:
                        change = round(cand * (self.get_population_count() / self.n_max) * cls.r_decr)

                    if change <= 0:
                        change = min(change, -3)
                    elif change > 0:
                        change = min(change, 100)

                    cand = cand + change
                    if cand > 0:
                        new_counts.append(cand)

                for _ in range(self.cand_max - self.get_population_size()):
                    de_novo_event = np.random.random()
                    if de_novo_event < cls.de_novo_event_prob:
                        new_counts.append(np.random.randint(1, 20))

                cls.candidate_counts = new_counts    
              #  self.timeseries[cls.name].append(cls.get_size())
                self.timeseries[cls.name].append(cls.get_count())

    def show_results(self):
        fig, axs = plt.subplots(1,1)
        for k, v in self.timeseries.items():
            axs.plot(v, label=k)

        plt.xlabel("time point")
        plt.ylabel("number of candidates")
        plt.ylim(bottom=0)
        plt.legend()
        plt.show()


if __name__ == "__main__":
    steps = 200
    use_abs = True
    if use_abs:
        Gain = DIPClass("gain", 261, 0.5012221023325012, 1.7862913959727234, -0.49710422318156694, 0.08333851479648849)
        Loss = DIPClass("loss", 373, 0.47635938952884715, 1.7891287030150624, -0.6390361485007808, 0.08504530568538506)
    else:
        Gain = DIPClass("gain", 261, 0.5225600824272179, 1.3383349229859336, -0.41407475207207134, 0.08449292411478025)
        Loss = DIPClass("loss", 373, 0.4782325094182225, 1.3973153622239822, -0.5737954358770873, 0.07461418881440231)

    Simulator = DIPPopulationSimulator([Gain, Loss], 100000, 1400)
    Simulator.run_simulation(steps)
    Simulator.show_results()

