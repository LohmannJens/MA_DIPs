'''
Loads the Start and End points of the deletion sides from Alnaji 2019 and gives
insights about the data distribution.

1. Creates a histogram for each line in each strain containing the length of
   the deletion sides multiplied by their occurence.

2. Plots length of Start and End part of DI RNA as a scatter plot. Shows if
   they are equally distributed.
'''

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import load_alnaji_excel, load_short_reads, get_sequence, get_stat_symbol


def plot_deletion_lengths(data: dict)-> None:
    '''
        creates a histogram for each strain, indicating the length of the
        deletions.
        :param data: dictionary with a data frame for each strain

        :return: None
    '''
    for key, value in data.items():
        # create a dict for each segment including the NGS read count
        count_dict = dict()
        for s in SEGMENTS:
            count_dict[s] = dict()

        for i, r in value.iterrows():
            if r["Length"] in count_dict[r["Segment"]]:
                count_dict[r["Segment"]][r["Length"]] += r["NGS_read_count"]
            else:
                count_dict[r["Segment"]][r["Length"]] = r["NGS_read_count"]

        # create a subplot for each key, value pair in count_dict
        fig, axs = plt.subplots(8, 1, figsize=(10, 20), tight_layout=True)
        fig.suptitle(f"absolute occurrences of deletions for {key}", x=0.3)
        for i, s in enumerate(SEGMENTS):
            axs[i].hist(count_dict[s].keys(), weights=count_dict[s].values(), bins=100)
            axs[i].set_title(f"{s}")
            axs[i].set_xlim(left=0)
            axs[i].set_xlabel("deletion length")

        save_path = os.path.join(RESULTSPATH, "deletion_length_and_position", f"{key}_length_del_hist.pdf")
        plt.savefig(save_path)
        plt.close()


def start_vs_end_lengths()-> None:
    '''

    '''


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    plot_deletion_lengths(all_reads_dict)
    

