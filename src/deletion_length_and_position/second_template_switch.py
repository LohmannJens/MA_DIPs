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
from utils import load_pelz_dataset, load_short_reads, load_alnaji_excel, load_WSN_data, load_full_alnaji2021


def long_di_rna_comparision(d: dict):
    '''
    
    '''
    for k, v in d.items():
        df = v[v["NGS_read_count"] >= 30].copy()
        print(k)
        for s in ["PB2", "PB1", "PA"]:
            df_s = df[df["Segment"] == s].copy()
            short_DIs = df_s[(df_s["Start"] <1000) & (df_s["End"] > 1000)]

            long_DIs = df_s[~((df_s["Start"] <1000) & (df_s["End"] > 1000))].copy()
            long_DIs["deletion_len"] = long_DIs["End"] - long_DIs["Start"]



            short_DIs_n = len(short_DIs)
            long_DIs_n = len(df_s) - short_DIs_n

            print(f"{s}\t{long_DIs_n/len(df_s)}")
            print(long_DIs["deletion_len"].mean())
            print(long_DIs["deletion_len"].median())


def run_comparision_all(ds: list, dnames: list):
    '''
    
    '''
    for d, name in zip(ds, dnames):
        print(name)
        long_di_rna_comparision(d)


if __name__ == "__main__":
    ds = list()
    dnames = list()

    ds.append(load_pelz_dataset(long_dirna=True))
    dnames.append("Pelz")

    ds.append(load_short_reads(load_alnaji_excel()))
    dnames.append("Alnaji2019")

    ds.append(load_full_alnaji2021())
    dnames.append("Alnaji2021")

    ds.append(load_WSN_data("Mendes"))
    dnames.append("Mendes")

    ds.append(load_WSN_data("Boussier"))
    dnames.append("Boussier")

    run_comparision_all(ds, dnames)
