"""

"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from scipy import stats

sys.path.insert(0, "..")
sys.path.insert(0, "../density_and_length_analysis")
from utils import RESULTSPATH, SEGMENTS
from utils import load_pelz_dataset, get_stat_symbol
from composition_junction_site import create_sequence_library, nucleotide_occurrence_analysis, count_overlapping_nucleotides_overall, generate_sampling_data


def nuc_overlap_analysis(seq_dict: dict, mode: int)-> None:
    '''
        Calculates the overlapping nucleotides of all sequences of the Pelz
        dataset.
        :param seq_dict: dictionary with the sequences
        :param mode: states which calculation mode is used. Check
                     composition_junction_site.py for more info.

        :return: None
    '''
    print(f"{mode=}")
    fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)
    for k, v in seq_dict.items():
        j = 0
        for i, s in enumerate(SEGMENTS):
            v_s = v.loc[(v["Segment"] == s)]
            nuc_overlap_dict, _ = count_overlapping_nucleotides_overall(v_s, mode)
            n = len(v_s.index)
            if n <= 1:
                continue
    
            seq = v_s.iloc[0]["WholeSequence"]
            q = 0.20
            start = (int(v_s.Start.quantile(q)), int(v_s.Start.quantile(1-q)))
            end = (int(v_s.End.quantile(q)), int(v_s.End.quantile(1-q)))
            m = 5

            sampling_data = generate_sampling_data(seq, start, end, n*m)
            exp, _ = count_overlapping_nucleotides_overall(sampling_data, mode)
            x = list(nuc_overlap_dict.keys())
            h = np.array(list(nuc_overlap_dict.values()))
            h_exp = np.array(list(exp.values()))

            # test statistical significance
            f_obs = list()
            f_exp = list()
            for a in x:
                f_obs.extend([a]*nuc_overlap_dict[a])
                f_exp.extend([a]*exp[a])
            f_obs = np.array(f_obs)
            f_exp = np.array(f_exp)

            res = stats.mannwhitneyu(f_obs, f_exp)
            symbol = get_stat_symbol(res.pvalue)
            print(res.pvalue)

            axs[i%4, j].bar(x=x, height=h/h.sum(), width=-0.4, align="edge", label="observed")
            axs[i%4, j].bar(x=x, height=h_exp/h_exp.sum(), width=0.4, align="edge", label="expected")
            axs[i%4, j].set_xlabel("number of overlapping nucleotides")
            axs[i%4, j].set_ylabel("relative occurrence")
            axs[i%4, j].set_title(f"{s} (n={n}) {symbol}")
            axs[i%4, j].legend(loc="upper right")
            axs[i%4, j].set_ylim(bottom=0.0, top=1.0)

            if i == 3:
                j = 1

    fname = f"overlapping_nucleotides_{k}_mode{mode}.pdf"
    savepath = os.path.join(RESULTSPATH, "control_analysis", fname)
    plt.savefig(savepath)
    plt.close()


if __name__ == "__main__":
    data_dict = load_pelz_dataset()
    seq_list_dict = create_sequence_library(data_dict)

    for s in SEGMENTS:
        nucleotide_occurrence_analysis(seq_list_dict, s)
        src = os.path.join(RESULTSPATH, "relative_occurrence_nucleotides", f"PR_{s}.pdf")
        dst = os.path.join(RESULTSPATH, "control_analysis", f"PR_{s}_nucleotide_occurrence.pdf")
        os.rename(src, dst)

    nuc_overlap_analysis(seq_list_dict, 1)
    nuc_overlap_analysis(seq_list_dict, 2)

