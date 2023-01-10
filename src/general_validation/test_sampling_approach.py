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
from utils import RESULTSPATH, SEGMENTS, QUANT, S_ROUNDS
from utils import load_alnaji_excel, load_short_reads, get_sequence, get_stat_symbol, generate_sampling_data, create_sequence_library, load_kupke, load_alnaji_2021, load_pelz_dataset

sys.path.insert(0, "../direct_repeats")
from search_direct_repeats import count_direct_repeats_overall


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
        label = k
        for s in SEGMENTS:
            v_s = v.loc[(v["Segment"] == s)]
            n = len(v_s.index)
            if n <= 1:
                continue

            start = (int(v_s.Start.quantile(QUANT)), int(v_s.Start.quantile(1-QUANT)))
            end = (int(v_s.End.quantile(QUANT)), int(v_s.End.quantile(1-QUANT)))
            seq = get_sequence(k, s)

            means = list()
            starts = list()
            n = 100
            i = 0

            while True:
                sampling_data = generate_sampling_data(seq, start, end, n)
                starts = starts + sampling_data["Start"].tolist()
                means.append(sum(starts)/len(starts))
                i = i + 1
                if len(means) > 1:
                    # difference is 0.1 %
                    if abs(means[-1] - means[-2]) < np.mean(means) * 0.001:
                        break

            rounds = np.arange(n, i*n+1, n)
            axs.scatter(rounds, means, label=s)

        axs.set_xlabel("number of samples")
        axs.set_ylabel("mean of start positions")
        axs.set_title(f"{k}")
        axs.legend()

        fname = f"{k}_testing_samling.png"
        if author != "":
            fname = f"{k}_{author}_testing_samling.png"
        savepath = os.path.join(RESULTSPATH, "general_validation", fname)
        plt.savefig(savepath)
        plt.close()


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    sequence_list_dict = create_sequence_library(all_reads_dict)
    test_sampling_approach(sequence_list_dict)

    for name in ["Pelz", "Kupke", "Alnaji2021"]:
        if name == "Pelz":
            raw_dict = load_pelz_dataset()
        elif name == "Kukpe":
            raw_dict = load_kupke()
        elif name == "Alnaji2021":
            raw_dict = load_alnaji_2021()

        seq_list_dict = create_sequence_library(raw_dict)
        test_sampling_approach(seq_list_dict, author=name)

