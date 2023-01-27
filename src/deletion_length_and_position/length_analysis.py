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
                count_dict[r["Segment"]][r["Length"]] += 1
            else:
                count_dict[r["Segment"]][r["Length"]] = 1

        # create a subplot for each key, value pair in count_dict
        fig, axs = plt.subplots(8, 1, figsize=(10, 15), tight_layout=True)
        for i, s in enumerate(SEGMENTS):
            if len(count_dict[s].keys()) > 1:
                m = round(np.mean(list(count_dict[s].keys())), 2)
                axs[i].hist(count_dict[s].keys(), weights=count_dict[s].values(), bins=100, label=f"{s} (mean={m})")
                axs[i].set_xlim(left=0)
                axs[i].set_xlabel("deletion length")
                axs[i].set_ylabel("# occurrences")
                axs[i].legend()

        axs[0].set_title(f"Deletion length of the eight segments for {key} as histogram")

        plt.rc("font", size=16)
        save_path = os.path.join(RESULTSPATH, "deletion_length_and_position", f"{key}_length_del_hist.png")
        plt.savefig(save_path)
        plt.close()


def plot_start_and_end_positions(data: dict)-> None:
    '''
        Maps the NP density to the start and end position of the deletion
        sites.
        :param data: dict with information about start and end position
        :param density_data: dict with density data (key is segment name)
        :return: counts positions found in NGS data
    '''
    for k, v in data.items():
        # create a dict for each segment using Start and End
        start_dict = dict()
        end_dict = dict()
        for s in SEGMENTS:
            start_dict[s] = dict()
            end_dict[s] = dict()
        for i, r in v.iterrows():
            if r["Start"] in start_dict[r["Segment"]]:
                start_dict[r["Segment"]][r["Start"]] += 1
            else:
                start_dict[r["Segment"]][r["Start"]] = 1
            if r["End"] in end_dict[r["Segment"]]:
                end_dict[r["Segment"]][r["End"]] += 1
            else:
                end_dict[r["Segment"]][r["End"]] = 1
        
        # create a subplot for each key, value pair in count_dict
        fig, axs = plt.subplots(8, 1, figsize=(10, 15), tight_layout=True)
        for i, s in enumerate(SEGMENTS):
            if len(start_dict[s].keys()) > 1:
                axs[i].bar(start_dict[s].keys(), start_dict[s].values(), label=f"{s} start")
                axs[i].bar(end_dict[s].keys(), end_dict[s].values(), label=f"{s} end")

                axs[i].set_xlim(left=0)
                axs[i].set_xlabel("sequence position")
                axs[i].set_ylabel("# occurrences")
                axs[i].legend(bbox_to_anchor=(1.0, 1.0))

        axs[0].set_title(f"Position of deletion site on full sequences for {k}")

        plt.rc("font", size=16)
        save_path = os.path.join(RESULTSPATH, "deletion_length_and_position", f"{k}_start_and_end_positions.png")
        plt.savefig(save_path)
        plt.close()


def start_vs_end_lengths(data: dict,
                         limit: int=0
                         )-> None:
    '''
        Plots the length of the start against the length of the end of the DI
        RNA sequences as a scatter plot.
        :param data: dictionary with a data frame (value)  for each strain (key)
        :param limit: indicates the length that is displayed in the plot
                      if equals 0 full lengths are displayed

        :return: None
    '''
    with open(os.path.join(DATAPATH, "Pelz2021", "packaging_signal.json"), "r") as f:
        packaging_signals = json.load(f)

    for k, v in data.items():
        fig, axs = plt.subplots(2, 4, figsize=(12, 6), tight_layout=True)
        j = 0
        for i, s in enumerate(SEGMENTS):
            v_s = v[v["Segment"] == s].copy()
            if v_s.shape[0] > 1:
                start = v_s["Start"]
                v_s["End_L"] = v_s["Length"] - v_s["Start"]

                axs[j,i%4].scatter(v_s["Start"], v_s["End_L"], s=1.0)
                axs[j,i%4].plot([0, 1], [0, 1], transform=axs[j,i%4].transAxes, c="r", linewidth=0.5, linestyle="--")

                if limit == 0:
                    max_p = max(v_s["Start"].max(), v_s["End_L"].max())
                else:
                    max_p = limit
                    axs[j,i%4].set_xlim(0, max_p)
                    axs[j,i%4].set_ylim(0, max_p)
                    v_s = v_s[(v_s["Start"] <= max_p) & (v_s["End_L"] <= max_p)]

                pearson = stats.pearsonr(v_s["Start"], v_s["End_L"])

                axs[j,i%4].set_xlim(left=0, right=max_p)
                axs[j,i%4].set_xticks([0,max_p])
                axs[j,i%4].set_yticks([0,max_p])
                axs[j,i%4].set_aspect("equal", "box")

                signals = packaging_signals[s]
                axs[j,i%4].add_patch(plt.Rectangle((0, 0), max_p, signals["incorporation_start"], color="grey", alpha = 0.5, ls="None"))
                axs[j,i%4].add_patch(plt.Rectangle((0, 0), signals["incorporation_end"], max_p, color="grey", alpha = 0.5, ls="None"))
                axs[j,i%4].add_patch(plt.Rectangle((0, 0), max_p, signals["bundling_start"], color="lightgrey", alpha = 0.5, ls="None"))
                axs[j,i%4].add_patch(plt.Rectangle((0, 0), signals["bundling_end"], max_p, color="lightgrey", alpha = 0.5, ls="None"))

            axs[j,i%4].set_title(f"{s} (n={v_s.shape[0]}) r={pearson[0]:.2}")
            axs[j,i%4].set_xlabel("start")
            axs[j,i%4].set_ylabel("end")

            if i == 3:
                j = 1

        fig.suptitle(f"Lengths of start and end site {k}", x=0.5)
        if limit == 0:
            filename = f"{k}_length_start_end.png"
        else:
            filename = f"{k}_length_start_end_{limit}.png"

        save_path = os.path.join(RESULTSPATH, "deletion_length_and_position", filename)
        plt.savefig(save_path)
        plt.close()


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    plot_deletion_lengths(all_reads_dict)
    plot_start_and_end_positions(all_reads_dict)
    start_vs_end_lengths(all_reads_dict)
    start_vs_end_lengths(all_reads_dict, 700)
