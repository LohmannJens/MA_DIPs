'''

'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


sys.path.insert(0, "..")
from utils import RESULTSPATH, SEGMENTS, COLORS, NUCLEOTIDES, QUANT, N_SAMPLES


class DIPClass():
    def __init__(self, name, n_init, prob, r_increase, r_decrease):
        self.name = name
        self.n = n_init
        self.p = prob
        self.r_incr = r_increase
        self.r_decr = r_decrease
        self.candidate_counts = np.array([np.random.randint(1, 100) for _ in range(self.n)])

    def get_size(self):
        return len(self.candidate_counts)

class DIPPopulationSimulator():
    def __init__(self, classes):
        self.classes = classes
    
    def run_simulation(self, steps):
        for i in range(steps):
            print(f"### iteration {i} ###")
            for c in self.classes:
                new_counts = list()
                for cand in c.candidate_counts:
                    thresh = np.random.random()
                    if c.p > thresh:
                        change = round(cand * c.r_incr)
                    else:
                        change = round(cand * c.r_decr)

                    if change < 0:
                        change = min(change, -3)
                    elif change > 0:
                        change = min(change, 10000)

                    cand = cand + change
                    if cand > 0:
                        new_counts.append(cand)
                c.candidate_counts = new_counts

                print(c.get_size())


if __name__ == "__main__":
    # Define rates and initialise starting population
    steps = 100
    Gain = DIPClass("gain", 261, 0.52, 1.0, -0.5)
    Loss = DIPClass("loss", 373, 0.46, 1.0, -0.5)

    Simulator = DIPPopulationSimulator([Gain, Loss])
    Simulator.run_simulation(steps)

