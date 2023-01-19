'''
    This script test the random sampling approach. It evaluates how big the
    samples size has to be to give reproducable and therefore comparable
    results.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, "..")
from utils import RESULTSPATH, SEGMENTS, QUANT
from utils import load_alnaji_excel, load_short_reads, get_sequence, get_stat_symbol, generate_sampling_data, create_sequence_library, load_kupke, load_alnaji_2021, load_pelz_dataset

sys.path.insert(0, "../direct_repeats")
from search_direct_repeats import count_direct_repeats_overall


def plot_distribution(starts_dict: dict,
                      name: str
                      )-> None:
    '''
        Gets the sampled starting positions for each segment and shows the
        distribution of them on the sequence.
        :param starts_dict: segment as key and list with positions as value
        :param name: prefix for the filename including straing (and author)

        :return: None
    '''
    fig, axs = plt.subplots(4, 2, figsize=(5, 7), tight_layout=True)
    i = 0
    j = 0
    for k, v in starts_dict.items():
        if len(v) == 0 or min(v) == max(v):
            continue
        axs[i,j].hist(v, bins=max(v)-min(v))
        axs[i,j].set_xlabel("position on sequence")
        axs[i,j].set_ylabel("count")
        axs[i,j].set_xlim(left=min(v), right=max(v))
        axs[i,j].set_title(f"{name}")
        
        i = i + 1
        if i == 4:
            j = 1
            i = 0

    fig.suptitle("Distribution of randomly sampled start positions")
    fname = f"{name}_distribution_sampling.png"
    savepath = os.path.join(RESULTSPATH, "general_validation", fname)
    plt.savefig(savepath)
    plt.close()


def test_sampling_approach(seq_dict: dict,
                           author: str=""
                           )-> None:
    '''
        Tests the sampling approach by generating increasingly bigger sets and
        comparing the mean of the start positions. When the difference of the
        mean of the last two batches are smaller than 0.1 % the testing is
        stopped.
        :param seq_dict: dictionary including data as value and strain name as
                         key
        :param author: author of the data set (needed for PR8 strains)

        :return: None
    '''
    plt.rc("font", size=12)
    for k, v in seq_dict.items():
        fig, axs = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)
        starts_dict = dict()
        for s in SEGMENTS:
            v_s = v.loc[(v["Segment"] == s)]
            n = len(v_s.index)
            if n <= 1:
                starts_dict[s] = list()
                continue

            start = (int(v_s.Start.quantile(QUANT)), int(v_s.Start.quantile(1-QUANT)))
            end = (int(v_s.End.quantile(QUANT)), int(v_s.End.quantile(1-QUANT)))
            seq = get_sequence(k, s)

            means = list()
            starts = list()
            n = 100
            i = 0
            thresh = -1

            while True:
                sampling_data = generate_sampling_data(seq, start, end, n)
                starts = starts + sampling_data["Start"].tolist()
                means.append(sum(starts)/len(starts))
                i = i + 1
                if len(means) > 1:
                    # difference is 0.1 %
                    if abs(means[-1] - means[-2]) < np.mean(means) * 0.001 and thresh == -1:
                        thresh = i
                    if i == 20:
                        break

            starts_dict[s] = starts

            rounds = np.arange(n, i*n+1, n)
            axs.scatter(rounds, means, label=s)
            axs.scatter(thresh*n, means[thresh-1], c="black", marker="x")

        axs.set_xlabel("number of samples")
        axs.set_ylabel("mean of start positions")
        name = k if author == "" else f"{author}_{k}"
        axs.set_title(f"{name}")
        axs.legend()

        fname = f"{name}_testing_sampling.png"

        savepath = os.path.join(RESULTSPATH, "general_validation", fname)
        plt.savefig(savepath)
        plt.close()

        plot_distribution(starts_dict, name)


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    sequence_list_dict = create_sequence_library(all_reads_dict)
    test_sampling_approach(sequence_list_dict, author="Alnaji2019")

    for name in ["Pelz", "Kupke", "Alnaji2021"]:
        if name == "Pelz":
            raw_dict = load_pelz_dataset()
        elif name == "Kukpe":
            raw_dict = load_kupke()
        elif name == "Alnaji2021":
            raw_dict = load_alnaji_2021()

        seq_list_dict = create_sequence_library(raw_dict)
        test_sampling_approach(seq_list_dict, author=name)

