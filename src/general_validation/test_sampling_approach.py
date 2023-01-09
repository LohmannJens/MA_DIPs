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
            n = 100
            i = 0
            starts = list()
            std = 0

            # using std as stopping criteria
            # is in first round 0.0 thats why i is also checked
      #      while (i < 10) or (std > 1.0):
            while True:
                sampling_data = generate_sampling_data(seq, start, end, n)
                starts = starts + sampling_data["Start"].tolist()
                means.append(sum(starts)/len(starts))
                std = np.std(means)
                i = i + 1
                if len(means) > 1:
                    # difference is 0.001 %
                    if abs(means[-1] - means[-2]) < np.mean(means) * 0.00001:
                        break

            rounds = np.arange(n, i*n+1, n)
            axs.scatter(rounds, means, label=s)
            axs.set_xlabel("samples")
            axs.set_ylabel("relative occurrence")
            axs.set_title(f"{s} (n={n})")
            axs.legend()

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
