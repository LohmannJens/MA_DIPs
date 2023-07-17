'''

'''
import os
import sys
import copy

import numpy as np
import matplotlib.pyplot as plt


sys.path.insert(0, "..")
from utils import RESULTSPATH


class DIPClass():
    def __init__(self, name):
        if name == "gain":
            self.name = name
            self.p = 0.5225600824272179
            self.r_incr = 1.3383349229859336
            self.r_decr = -0.41407475207207134
        elif name == "loss":
            self.name = name
            self.p = 0.4782325094182225
            self.r_incr = 1.3973153622239822
            self.r_decr = -0.5737954358770873
    
    def get_class_name(self):
        return self.name
    
    def get_prob(self):
        return self.p
    
    def get_rate_increase(self):
        return self.r_incr
        
    def get_rate_decrease(self):
        return self.r_decr

    def get_de_novo_event_prob(self):
        return self.de_novo_event_prob

class DIPCandidate():
    def __init__(self, cls, count):
        self.cls = cls
        self.count = count

    def get_count(self):
        return self.count
    
    def get_class(self):
        return self.cls

class DIPPopulation():
    def __init__(self, candidates, count_max, n_cand_max):
        self.candidates = candidates
        self.count_max = count_max
        self.n_cand_max = n_cand_max
        
    def get_population_n_cand(self):
        return len(self.candidates)
    
    def get_gain_n_cand(self):
        return len([c for c in self.candidates if c.get_class().name == "gain"])
    
    def get_loss_n_cand(self):
        return len([c for c in self.candidates if c.get_class().name == "loss"])

    def get_population_count(self):
        return np.sum([c.get_count() for c in self.candidates])
    
    def get_gain_count(self):
        return np.sum([c.get_count() for c in self.candidates if c.get_class().name == "gain"])

    def get_loss_count(self):
        return np.sum([c.get_count() for c in self.candidates if c.get_class().name == "loss"])

class DIPPopulationSimulator():
    def __init__(self, population):
        self.recent_population = population
        self.population_timeseries = [population]

    def get_recent_population(self):
        return self.recent_population

    def run_simulation(self, steps):
        for i in range(steps):
            new_candidates_list = list()
            if i % 100 == 0:
                print(i)
            loss_event_p = 0.1
            loss_event = True if loss_event_p > np.random.random() else False

            count_max = self.get_recent_population().count_max

            for cand in self.get_recent_population().candidates:
                thresh = np.random.random()
                p = 0.1 if loss_event else cand.get_class().get_prob()

                if p > thresh:
                    change = round(cand.get_count() * (1 - self.get_recent_population().get_population_count() / count_max) * cand.get_class().get_rate_increase() * 0.1) # here I have an extra artificial factor included right now
                else:
                    change = round(cand.get_count() * (self.get_recent_population().get_population_count() / count_max) * cand.get_class().get_rate_decrease())

                if change <= 0:
                    change = min(change, -3)

                new_count = cand.get_count() + change
                if new_count > 0:
                    new_candidates_list.append(DIPCandidate(cand.get_class(), new_count))

            de_novo_event_prob = 0.1
            for _ in range(round((count_max - self.get_recent_population().get_population_count())/50)):
                de_novo_event = np.random.random()
                if de_novo_event < de_novo_event_prob:
                    class_p = np.random.random()
                    if class_p < 0.5:
                        cls = "gain"
                    else:
                        cls = "loss"
                    
                    new_candidates_list.append(DIPCandidate(DIPClass(cls), np.random.randint(1, 20)))

            n_cand_max = self.get_recent_population().n_cand_max

            self.recent_population = DIPPopulation(new_candidates_list, count_max, n_cand_max)
            self.population_timeseries.append(copy.deepcopy(self.get_recent_population()))

    def show_results(self):
        fig, axs = plt.subplots(1,2)
        t = list()
        counts = list()
        gain_counts = list()
        loss_counts = list()
        gain_n_cand = list()
        loss_n_cand = list()

        for i, pop_t in enumerate(self.population_timeseries):
            t.append(i)
            counts.append(pop_t.get_population_count())
            gain_counts.append(pop_t.get_gain_count())
            loss_counts.append(pop_t.get_loss_count())
            gain_n_cand.append(pop_t.get_gain_n_cand())
            loss_n_cand.append(pop_t.get_loss_n_cand())

        axs[0].plot(t, gain_n_cand, label="gain")
        axs[0].plot(t, loss_n_cand, label="loss")
        axs[0].set_xlabel("time point")
        axs[0].set_ylabel("number of candidates")
        axs[0].set_ylim(bottom=0)
        axs[0].legend()

        axs[1].plot(t, gain_counts, label="gain")
        axs[1].plot(t, loss_counts, label="loss")
        axs[1].set_xlabel("time point")
        axs[1].set_ylabel("count of candidates")
        axs[1].set_ylim(bottom=0)
        axs[1].legend()

        plt.show()

    def show_candidate_developement(self):
        x = list()
        for i, p in enumerate(self.population_timeseries):
            if i % round(len(self.population_timeseries)/10) == 0:
                x.append([c.get_count() for c in p.candidates])

        fig, axs = plt.subplots(1,1)
        axs.boxplot(x)
        axs.set_xlabel("time point")
        axs.set_ylabel("number of counts")
        axs.legend()        

        plt.show()        


if __name__ == "__main__":
    steps = 1000
       
    candidates = list()

    for _ in range(261):
        candidates.append(DIPCandidate(DIPClass("gain"), np.random.randint(1, 20)))
    for _ in range(373):
        candidates.append(DIPCandidate(DIPClass("loss"), np.random.randint(1, 20)))

    start_population = DIPPopulation(candidates, 50000, 1400)

    Simulator = DIPPopulationSimulator(start_population)
    Simulator.run_simulation(steps)
    Simulator.show_results()
    Simulator.show_candidate_developement()
