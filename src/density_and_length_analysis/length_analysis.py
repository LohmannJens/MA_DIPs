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
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import load_alnaji_excel, load_short_reads


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


def start_vs_end_lengths(data: dict)-> None:
    '''
        Plots the length of the start against the length of the end of the DI
        RNA sequences as a scatter plot.
        :param data: dictionary with a data frame (value)  for each strain (key)

        :return: None
    '''
    with open(os.path.join(DATAPATH, "Pelz2021", "packaging_signal.json"), "r") as f:
        packaging_signals = json.load(f)

    for k, v in data.items():
        fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)
        j = 0
        for i, s in enumerate(SEGMENTS):
            v_s = v[v["Segment"] == s]
            if v_s.size != 0:
                start = v_s["Start"]
                end = v_s["Length"] - v_s["Start"]

                axs[i%4,j].scatter(start, end, s=1.0)
                axs[i%4,j].plot([0, 1], [0, 1], transform=axs[i%4,j].transAxes, c="r", linewidth=0.5, linestyle="--")

                max_p = max(start.max(), end.max())
                axs[i%4,j].set_xlim(left=0, right=max_p)
                axs[i%4,j].set_xticks([0,max_p])
                axs[i%4,j].set_yticks([0,max_p])
                axs[i%4,j].set_aspect("equal", "box")

                signals = packaging_signals[s]
                axs[i%4,j].add_patch(plt.Rectangle((0, 0), max_p, signals["incorporation_start"], color="grey", alpha = 0.5, ls="None"))
                axs[i%4,j].add_patch(plt.Rectangle((0, 0), signals["incorporation_end"], max_p, color="grey", alpha = 0.5, ls="None"))
                axs[i%4,j].add_patch(plt.Rectangle((0, 0), max_p, signals["bundling_start"], color="lightgrey", alpha = 0.5, ls="None"))
                axs[i%4,j].add_patch(plt.Rectangle((0, 0), signals["bundling_end"], max_p, color="lightgrey", alpha = 0.5, ls="None"))

            axs[i%4,j].set_title(f"{s} (n={v_s.shape[0]})")
            axs[i%4,j].set_xlabel("start")
            axs[i%4,j].set_ylabel("end")

            if i == 3:
                j = 1

        fig.suptitle(f"relation of length of start and end site {k}", x=0.3)

        save_path = os.path.join(RESULTSPATH, "deletion_length_and_position", f"{k}_length_start_end.pdf")
        plt.savefig(save_path)
        plt.close()


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
#    plot_deletion_lengths(all_reads_dict)
    start_vs_end_lengths(all_reads_dict)   

