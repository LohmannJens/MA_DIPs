'''
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, "..")
from utils import RESULTSPATH, SEGMENTS, QUANT, S_ROUNDS
from utils import load_alnaji_excel, load_short_reads, get_sequence, get_stat_symbol, generate_sampling_data, create_sequence_library

sys.path.insert(0, "../direct_repeats")
from search_direct_repeats import count_direct_repeats_overall


def test_sampling_approach(seq_dict: dict)-> None:
    '''
    '''
    plt.rc("font", size=12)
    for k, v in seq_dict.items():
        fig, axs = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)
        label = k
        SEGMENTS = ["PB1"]
        for s in SEGMENTS:
            v_s = v.loc[(v["Segment"] == s)]
            n = len(v_s.index)
            if n <= 1:
                continue

            start = (int(v_s.Start.quantile(QUANT)), int(v_s.Start.quantile(1-QUANT)))
            end = (int(v_s.End.quantile(QUANT)), int(v_s.End.quantile(1-QUANT)))
            seq = get_sequence(k, s)

            means = list()
            dir_rep = list()
            n = 100
            rounds = np.arange(1, 10000, n)
            
            dir_rep_dict = dict()
            for r in rounds:
                sampling_data = generate_sampling_data(seq, start, end, n)
                new_dict, _ = count_direct_repeats_overall(sampling_data, seq, 1)
                dir_rep_dict = {i: dir_rep_dict.get(i, 0) + new_dict.get(i, 0) for i in set(dir_rep_dict).union(new_dict)}
                num = sum(dir_rep_dict.values())
                means.append(sum([key * value for key, value in dir_rep_dict.items()])/num)

            axs.scatter(rounds, means, label=s)
            axs.set_xlabel("rounds")
            axs.set_ylabel("relative occurrence")
            axs.set_title(f"{s} (n={n})")
            axs.legend()
            plt.ylim([0,1])

        plt.show()

        '''
        fname = f"{k}.png"
        savepath = os.path.join(RESULTSPATH, "general_validation", fname)
        plt.savefig(savepath)
        plt.close()
        '''

if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    sequence_list_dict = create_sequence_library(all_reads_dict)

    test_sampling_approach(sequence_list_dict)
