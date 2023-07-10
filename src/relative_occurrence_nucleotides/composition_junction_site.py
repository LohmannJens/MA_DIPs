'''
    Loads the data for the deletion sides (start/end point) from Alnaji 2019.
    Does different analysis with the data.

    Takes a look at the nucleotide distribution around these points.
    Goal is to see if any of the four nucleotides occure more/less often.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import Tuple
from scipy import stats
from decimal import Decimal, ROUND_HALF_UP

sys.path.insert(0, "..")
from utils import RESULTSPATH, SEGMENTS, COLORS, NUCLEOTIDES, QUANT, N_SAMPLES
from utils import load_alnaji_excel, load_short_reads, get_sequence, get_stat_symbol, generate_sampling_data, create_sequence_library


def count_nucleotide_occurrence(seq: str,
                                p: int
                                )-> dict:
    '''
        Counts the number of nucleotides next to a given point.
        Goes 5 steps in both directions.
        :param seq: whole RNA sequence
        :param p: point on the sequence where to count

        :return: returns a counter dict with an entry for each nucleotide. In
                 each entry the counter for each position is given.
    '''
    window = seq[p-5:p+5]
    r_dict = dict({n: np.zeros(10) for n in NUCLEOTIDES})

    for i, char in enumerate(window):
        r_dict[char][i] = 1
    return r_dict


def count_nucleotide_occurrence_overall(df: pd.DataFrame,
                                        seq: str
                                        )-> Tuple[dict, dict]:
    '''
        Counts the occurrence of each nucleotide at different positions around
        the junction site
        :param df: dataframe with sequence and junction site data
        :param seq: rna sequence where to count the occurrence

        :return: tupel with two entries:
                    dict with nucleotide count for start site
                    dict with nucleotide count for end site
    '''

    count_start_dict = dict({n: np.zeros(10) for n in NUCLEOTIDES})
    count_end_dict = dict({n: np.zeros(10) for n in NUCLEOTIDES})
    normalize = 0

    for i, row in df.iterrows():
        seq_start_dict = count_nucleotide_occurrence(seq, row["Start"]) 
        seq_end_dict = count_nucleotide_occurrence(seq, row["End"]-1)
        normalize += 1
        for nuc in count_start_dict.keys():
            count_start_dict[nuc] += seq_start_dict[nuc]
            count_end_dict[nuc] += seq_end_dict[nuc]

    return count_start_dict, count_end_dict


def nucleotide_occurrence_analysis(seq_dict: dict,
                                   seg: str,
                                   author: str=""
                                   )-> None:
    '''
        Gets the sequences for all four strains and calculates the occurrence
        of each nucleotide at the start and end deletion site.
        :param seq_dict: dictionary with the sequences
        :param seg: name of the segment that is analyzed
        :param author: authors name, used to distinguish PR8 datasets

        :return: None
    '''
    for k, v in seq_dict.items():
        v = v.loc[v["Segment"] == seg]
        seq = get_sequence(k, seg)
        count_start_dict, count_end_dict = count_nucleotide_occurrence_overall(v, seq)
        n = len(v.index)
        # only plot results if at least one data point is available
        if n == 0:
            continue

        # get expected values
        s = (int(v.Start.quantile(QUANT)), int(v.Start.quantile(1-QUANT)))
        e = (int(v.End.quantile(QUANT)), int(v.End.quantile(1-QUANT)))

        sampling_data = generate_sampling_data(seq, s, e, N_SAMPLES)
        exp_s, exp_e = count_nucleotide_occurrence_overall(sampling_data, seq)
        
        fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)
        x = np.arange(0.8, 10.8, dtype=np.float64)

        for idx, nuc in enumerate(count_start_dict.keys()):
            h_s = count_start_dict[nuc]/n
            h_e = count_end_dict[nuc]/n
            y_exp_s = exp_s[nuc] / (N_SAMPLES)
            y_exp_e = exp_e[nuc] / (N_SAMPLES)

            axs[idx, 0].bar(x, height=h_s, width=0.3, label=nuc, color=COLORS[nuc])
            axs[idx, 1].bar(x, height=h_e, width=0.3, label=nuc, color=COLORS[nuc])
            axs[idx, 0].bar(x+0.4, height=y_exp_s, width=0.3, label=f"{nuc}_exp", color=COLORS[nuc], alpha=0.5)
            axs[idx, 1].bar(x+0.4, height=y_exp_e, width=0.3, label=f"{nuc}_exp", color=COLORS[nuc], alpha=0.5)

            for j, p in enumerate(y_exp_s):
                result = stats.binomtest(int(count_start_dict[nuc][j]), n, p)
                symbol = get_stat_symbol(result.pvalue)
                axs[idx, 0].annotate(symbol, (j+1, max(h_s[j], p)), fontsize="x-small",
                                     ha="center", stretch="condensed")
            for j, p in enumerate(y_exp_e):
                result = stats.binomtest(int(count_end_dict[nuc][j]), n, p)
                symbol = get_stat_symbol(result.pvalue)
                axs[idx, 1].annotate(symbol, (j+1, max(h_e[j], p)), fontsize="x-small",
                                     ha="center", stretch="condensed")

            for i in range(2):
                axs[idx, i].margins(x=0)
                axs[idx, i].set_xlim(left=0.5, right=10.5)
                axs[idx, i].set_ylim(top=0.8, bottom=0.0)
                axs[idx, i].set_xticks([1,2,3,4,5,6,7,8,9,10])
                axs[idx, i].set_xlabel("position at deletion side")
                axs[idx, i].set_ylabel("relative occurrence")
            axs[idx, 0].add_patch(plt.Rectangle((5.5, 0), 5, 1, color="grey", alpha=0.3))
            axs[idx, 1].add_patch(plt.Rectangle((0.5, 0), 5, 1, color="grey", alpha=0.3))
  
        by_label = dict()
        for ax in axs:
            for a in ax:
                handles, labels = a.get_legend_handles_labels()
                by_label.update(dict(zip(labels, handles)))
        fig.legend(by_label.values(), by_label.keys(), ncol=4, loc="upper center")
        fig.suptitle(f"\n\n\n{seg} ({n})")

        fname = f"{k}_{seg}.png"
        if author != "":
            fname = f"{author}_{fname}"
        savepath = os.path.join(RESULTSPATH, "relative_occurrence_nucleotides", fname)
        plt.savefig(savepath)
        plt.close()


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    sequence_list_dict = create_sequence_library(all_reads_dict)

    # Loop over the different strains and calculate the occurrence of each
    # nucleotide in the sequences
    for s in SEGMENTS:
        nucleotide_occurrence_analysis(sequence_list_dict, s)

