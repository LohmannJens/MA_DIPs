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

    def get_n_cand(self):
        return len(self.candidate_counts)

    def get_count(self):
        return np.sum(self.candidate_counts)

class DIPPopulationSimulator():
    def __init__(self, classes, count_max, n_cand_max):
        self.classes = classes
        self.count_max = count_max
        self.n_cand_max = n_cand_max
        self.n_cand_timeseries = dict({c.name : list() for c in self.classes})
        self.count_timeseries = dict({c.name : list() for c in self.classes})
    
    def get_population_n_cand(self):
        return np.sum([c.get_n_cand() for c in self.classes])
    
    def get_population_count(self):
        return np.sum([c.get_count() for c in self.classes])

    def run_simulation(self, steps):
        for i in range(steps):
            if i % 100 == 0:
                print(i)
            loss_event_p = 0.1
            loss_event = True if loss_event_p > np.random.random() else False

            for cls in self.classes:
                new_counts = list()
                for cand_counts in cls.candidate_counts:
                    thresh = np.random.random()

                    p = 0.1 if loss_event else cls.p

                    if p > thresh:
                        change = round(cand_counts * (1 - self.get_population_count() / self.count_max) * cls.r_incr * 0.1) # here I have an extra artificial factor included right now
                    else:
                        change = round(cand_counts * (self.get_population_count() / self.count_max) * cls.r_decr)

                    if change <= 0:
                        change = min(change, -3)

                    cand_counts = cand_counts + change
                    if cand_counts > 0:
                        new_counts.append(cand_counts)

                for _ in range(round((self.count_max - self.get_population_count())/50)):
                    de_novo_event = np.random.random()
                    if de_novo_event < cls.de_novo_event_prob:
                        new_counts.append(np.random.randint(1, 20))

                cls.candidate_counts = new_counts    
                self.n_cand_timeseries[cls.name].append(cls.get_n_cand())
                self.count_timeseries[cls.name].append(cls.get_count())

    def show_results(self):
        fig, axs = plt.subplots(1,2)
        for k, v in self.n_cand_timeseries.items():
            axs[0].plot(v, label=k)
            axs[0].set_xlabel("time point")
            axs[0].set_ylabel("number of candidates")
            axs[0].set_ylim(bottom=0)
            axs[0].legend()
        for k, v in self.count_timeseries.items():
            axs[1].plot(v, label=k)
            axs[1].set_xlabel("time point")
            axs[1].set_ylabel("count of candidates")
            axs[1].set_ylim(bottom=0)
            axs[1].legend()

        plt.show()

    def show_candidate_developement(self):
        fig, axs = plt.subplots(1,2)
        for k, v in self.n_cand_timeseries.items():
            axs[0].plot(v, label=k)
            axs[0].set_xlabel("time point")
            axs[0].set_ylabel("number of candidates")
            axs[0].set_ylim(bottom=0)
            axs[0].legend()
        for k, v in self.count_timeseries.items():
            axs[1].plot(v, label=k)
            axs[1].set_xlabel("time point")
            axs[1].set_ylabel("count of candidates")
            axs[1].set_ylim(bottom=0)
            axs[1].legend()

        plt.show()        



if __name__ == "__main__":
    steps = 1000
    use_abs = False
    if use_abs:
        Gain = DIPClass("gain", 261, 0.5012221023325012, 1.7862913959727234, -0.49710422318156694, 0.0396850070459469)
        Loss = DIPClass("loss", 373, 0.47635938952884715, 1.7891287030150624, -0.6390361485007808, 0.044547541073296934)
    else:
        Gain = DIPClass("gain", 261, 0.5225600824272179, 1.3383349229859336, -0.41407475207207134, 0.044258198345837274)
        Loss = DIPClass("loss", 373, 0.4782325094182225, 1.3973153622239822, -0.5737954358770873, 0.039083622712305977)

    Simulator = DIPPopulationSimulator([Gain, Loss], 50000, 1400)
    Simulator.run_simulation(steps)
    Simulator.show_results()

