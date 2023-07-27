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

from scipy import stats

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, STRAINS
from utils import load_pelz_dataset


def long_di_rna_comparision(d: dict):
    '''
    
    '''
    for k, v in d.items():
        df = v[v["NGS_read_count"] >= 30].copy()
        for s in ["PB2", "PB1", "PA"]:
            df_s = df[df["Segment"] == s].copy()
            short_DIs = df_s[(df_s["Start"] <1000) & (df_s["End"] > 1000)]

            
            short_DIs_n = len(short_DIs)
            long_DIs_n = len(df_s) - short_DIs_n

            print(f"{s}\t{long_DIs_n/len(df_s)}")

if __name__ == "__main__":
    pelz_dict = load_pelz_dataset(long_dirna=True)
    long_di_rna_comparision(pelz_dict)
