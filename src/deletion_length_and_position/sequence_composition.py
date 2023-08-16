'''
    Loads the Start and End points of the deletion sides from Alnaji 2019 and
    gives insights about the data distribution.

    1.
    Creates a histogram for each line in each strain containing the length of
    the deletion sides multiplied by their occurence.
    
    2.
    Creates a plot where it shows the location of the start and the end points
    of the deletion site in reference to the full length segments

    3.
    Plots length of Start and End part of DI RNA as a scatter plot. Shows if
    they are equally distributed.
'''
import os
import sys
import json

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, STRAINS, NUCLEOTIDES, COLORS
from utils import get_sequence, get_seq_len


def plot_sequence_composition(strain: str):
    '''
    
    '''
    print(f"+++ {strain} +++")
    thresh_dict = dict({
        "PB2": 800, "PB1": 600, "PA": 700, "HA": 700, "NP": 600, "NA": 600, "M": 900, "NS": 500
    })

    overall = dict()
    hotspot = dict()
    for s in ["PB2", "PB1", "PA"]:
        seq = get_sequence(strain, s, full=True)
        seq = seq.seq.transcribe()
        for n in NUCLEOTIDES:
            overall[n] = seq.count(n)
            hotspot[n] = seq[0:thresh_dict[s]].count(n)
                
        overall_n = sum(overall.values())
        hotspot_n = sum(hotspot.values())
        for k in overall.keys():
            print(k)
            print(overall[k]/overall_n)
            print(hotspot[k]/hotspot_n)


if __name__ == "__main__":
    
    plot_sequence_composition("PR8")
    #for st in STRAINS.keys():
     #   plot_sequence_composition(st)
