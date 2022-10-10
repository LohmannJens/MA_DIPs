'''

'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats
from matplotlib_venn import venn3
from decimal import Decimal, ROUND_HALF_UP

sys.path.insert(0, "..")
sys.path.insert(0, "../relative_occurrence_nucleotides")
sys.path.insert(0, "../direct_repeats")
from utils import RESULTSPATH, SEGMENTS, COLORS
from utils import load_alnaji_2021, get_sequence, get_stat_symbol, create_sequence_library
from composition_junction_site import count_nucleotide_occurrence_overall
from search_direct_repeats import count_direct_repeats_overall, include_correction


def venn_different_timepoints(data: dict)-> None:
    '''
        Draws a venn diagramm for a given dataset with the groups built by 
        different timepoints (3hpi, 6hpi, 24hpi)
        :param data: dict of the dataset with the strain name as key and a
                     data frame as value
        
        :return: None
    '''
    for k, v in data.items():
        fig, axs = plt.subplots(4, 1, figsize=(8, 10), tight_layout=True)
        for i, t in enumerate(["3hpi", "6hpi", "24hpi", "all"]):
            if t == "all":
                set1 = set(v[v["Timepoint"] == "3hpi"]["DI"])
                set2 = set(v[v["Timepoint"] == "6hpi"]["DI"])
                set3 = set(v[v["Timepoint"] == "24hpi"]["DI"])
                labels = ("3hpi", "6hpi", "24hpi")
            else:
                v_t = v[v["Timepoint"] == t].copy()
                set1 = set(v_t[v_t["Replicate"] == "Rep1"]["DI"])
                set2 = set(v_t[v_t["Replicate"] == "Rep2"]["DI"])
                set3 = set(v_t[v_t["Replicate"] == "Rep3"]["DI"])
                labels = ("Rep1", "Rep2", "Rep3")

            venn3([set1, set2, set3], set_labels=labels, ax=axs[i])
            axs[i].set_title(f"{t}")

        fig.suptitle(f"overlap of replicates at different timepoints for {k}")
        
        save_path = os.path.join(RESULTSPATH, "di_rna_conservation", f"venn_alnaji_timepoints.png")
        plt.savefig(save_path)
        plt.close()


def compare_nucleotide_occurrence(df: object, strain: str)-> None:
    '''
        gets the sequences dataset and calculates the occurrence of each
        nucleotide at the start and end deletion site.
        Compares two groups with each other:
            - DI candidates that occurred in more than [n] replicates
            - DI candidates that occurred in less than [n] replicates
        :param seq_dict: dictionary with the sequences
        :param seg: name of the strain is analyzed

        :return: None
    '''
    for s in SEGMENTS:
        fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)
        seq = get_sequence(strain, s)

        below_df = df.loc[(df["Group"] == "below") & (df["Segment"] == s)]
        above_df = df.loc[(df["Group"] == "above") & (df["Segment"] == s)]

        n_below = len(below_df)
        n_above = len(above_df)

        below_start, below_end = count_nucleotide_occurrence_overall(below_df, seq)
        above_start, above_end = count_nucleotide_occurrence_overall(above_df, seq)

        # only plot results if at least one data point is available
        # for alnaji 2021 and thresh == 2 it is only PB2
        if (n_below == 0 or n_above == 0):
            continue

        x = np.arange(0.8, 9.8, dtype=np.float64)

        for i, nuc in enumerate(below_start.keys()):
            h_below_s = below_start[nuc]/ (n_below)
            h_below_e = below_end[nuc]/ (n_below)
            h_above_s = above_start[nuc] / (n_above)
            h_above_e = above_end[nuc] / (n_above)

            axs[i, 0].bar(x, height=h_below_s, width=0.3,
                          label="below", color=COLORS[nuc])
            axs[i, 1].bar(x, height=h_below_e, width=0.3,
                          label="below", color=COLORS[nuc])
            axs[i, 0].bar(x+0.4, height=h_above_s, width=0.3,
                          label="above", color=COLORS[nuc], alpha=0.5)
            axs[i, 1].bar(x+0.4, height=h_above_e, width=0.3,
                          label="above", color=COLORS[nuc], alpha=0.5)
            axs[i, 0].set_title(nuc)
            axs[i, 1].set_title(nuc)

            # statisitcal testing
            # TODO: maybe delete it, because everything is significant
            for j, (k, p) in enumerate(zip(below_start[nuc], h_above_s)):
                result = stats.binomtest(int(k), n_below, p)
                symbol = get_stat_symbol(result.pvalue)
                axs[i, 0].annotate(symbol, (j+1, max(h_below_s[j], p)),
                                   fontsize="x-small", ha="center", stretch="condensed")
            for j, (k, p) in enumerate(zip(below_end[nuc], h_above_e)):
                result = stats.binomtest(int(k), n_below, p)
                symbol = get_stat_symbol(result.pvalue)
                axs[i, 1].annotate(symbol, (j+1, max(h_below_e[j], p)),
                                   fontsize="x-small", ha="center", stretch="condensed")

            for k in range(2):
                axs[i, k].legend()
                axs[i, k].margins(x=0)
                axs[i, k].set_xlim(left=0.5, right=9.5)
                axs[i, k].set_ylim(top=0.8, bottom=0.0)
                axs[i, k].set_xticks([1,2,3,4,5,6,7,8,9])
                axs[i, k].set_xlabel("position at junction side")
                axs[i, k].set_ylabel("relative occurrence")
            axs[i, 0].add_patch(plt.Rectangle((5.5, 0), 4, 1, color="grey", alpha=0.3))
            axs[i, 1].add_patch(plt.Rectangle((0.5, 0), 4, 1, color="grey", alpha=0.3))
  
        plt.suptitle(f"start (left) and end (right) of {s} of {strain}")
        savepath = os.path.join(RESULTSPATH,
                                "di_rna_conservation",
                                f"{s}_nucleotide_occurrence.png")
        plt.savefig(savepath)
        plt.close()


def compare_direct_repeats(df: dict, strain: str, mode: int, correction: bool)-> None:
    '''
        gets the sequences for all four strains and calculates the overlap of 
        the nucleotides at the junction site.
        :param df: data frame with the sequences
        :param strain: name of the strain
        :param mode: states which calculation mode is used in 
                     calculate_overlapping_nucleotides() check there for info
        :param correction: if True a correction calculation for the counts is
                           made

        :return: None
    '''
    fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)

    i = 0
    j = 0
    for s in SEGMENTS:
        seq = get_sequence(strain, s)
        below_df = df.loc[(df["Group"] == "below") & (df["Segment"] == s)]
        above_df = df.loc[(df["Group"] == "above") & (df["Segment"] == s)]
        n_below = len(below_df)
        n_above = len(above_df)
        below_direct_repeats, _ = count_direct_repeats_overall(below_df, seq, mode)
        above_direct_repeats, _ = count_direct_repeats_overall(above_df, seq, mode)

        # only plot results if at least one data point is available
        if (n_below <= 1 or n_above <= 1):
            continue
        
        if correction:
            below_direct_repeats = include_correction(below_direct_repeats)       
            above_direct_repeats = include_correction(above_direct_repeats)

        x = list(below_direct_repeats.keys())
        below_h = np.array(list(below_direct_repeats.values()))
        above_h = np.array(list(above_direct_repeats.values()))
        
        # test statistical significance
        f_below = list()
        f_above = list()
        for a in x:
            if correction:
                f_below.extend([a]*int(Decimal(below_direct_repeats[a]).to_integral_value(rounding=ROUND_HALF_UP)))
                f_above.extend([a]*int(Decimal(above_direct_repeats[a]).to_integral_value(rounding=ROUND_HALF_UP)))
            else:
                f_below.extend([a]*below_direct_repeats[a])
                f_above.extend([a]*above_direct_repeats[a])

        f_below = np.array(f_below)
        f_above = np.array(f_above)

        # select statistical test here
        #stats_test = "mannwhitneyu"
        stats_test = "ks_2samp"

        if stats_test == "mannwhitneyu":
           res = stats.mannwhitneyu(f_below, f_above)
           symbol = get_stat_symbol(res.pvalue)
        elif stats_test == "ks_2samp":
            stat, pvalue = stats.ks_2samp(f_below, f_above)
            symbol = get_stat_symbol(pvalue)
        
        # plot results as barplot
        axs[i, j].bar(x=x, height=below_h/below_h.sum(), width=-0.4, align="edge", label=f"below (n={n_below})")
        axs[i, j].bar(x=x, height=above_h/above_h.sum(), width=0.4, align="edge", label=f"above (n={n_above})")
        axs[i, j].set_xlabel("number of overlapping nucleotides")
        axs[i, j].set_ylabel("relative occurrence")
        axs[i, j].set_title(f"{s} {symbol}")
        axs[i, j].legend(loc="upper right")
        axs[i, j].set_ylim(bottom=0.0, top=1.0)

        if i == 3:
            i = 0
            j = 1
        else:
            i += 1


    corr = "_corr" if correction else ""    
    fname = f"alnaji_{strain}_mode{mode}_compare_direct_repeats{corr}.pdf"
    savepath = os.path.join(RESULTSPATH, "di_rna_conservation", fname)
    plt.savefig(savepath)
    plt.close()


def slice_by_occurrence(df: object, thresh: int, below: bool)-> object:
    '''
        Allows to slice a dataset by the number of occurrences for each DI
        candidate. Counts the number of occurrences for each DI Name and then
        prints the candidates out, that are above a given threshhold.
        :param data: data frame containing the DI candidates
        :param thresh: occurrences larger or equal to this parameter are
                       included in the written file

        :return: data frame
    '''
    occur_df = pd.DataFrame(df.groupby(["DI"]).size())
    occur_df = occur_df.rename(columns={0: "Occurrences"})

    if below:
        select_list = occur_df[occur_df["Occurrences"] < thresh].index.values.tolist()
    elif not below:
        select_list = occur_df[occur_df["Occurrences"] >= thresh].index.values.tolist()

    selected_df = df[df["DI"].isin(select_list)].copy()
    selected_df = selected_df.groupby("DI").sum() # combine NGS_read_count
    selected_df.reset_index(inplace=True)

    return_df = selected_df[["DI", "NGS_read_count"]].copy()
    return_df[["Segment", "Start", "End"]] = return_df["DI"].str.split("_", expand=True)
    return_df.drop(["DI"], axis=1, inplace=True)
    return_df = return_df[["Segment", "Start", "End", "NGS_read_count"]] # reorder cols

    return_df["Start"] = return_df["Start"].astype(str).astype(int)
    return_df["End"] = return_df["End"].astype(str).astype(int)
    return_df["Segment"] = return_df["Segment"].astype(str)

    return return_df


if __name__ == "__main__":
    data_dict = load_alnaji_2021()
    # check the overlap of the different timepoints
    venn_different_timepoints(data_dict)

    # Compare DI candidates that occur once to those that occur more often
    data_df = data_dict["PR8"]

    below_df = slice_by_occurrence(data_df, 2, below=True)
    above_df = slice_by_occurrence(data_df, 2, below=False)

    below_df = below_df.assign(Group = ["below"] * len(below_df.index))
    above_df = above_df.assign(Group = ["above"] * len(above_df.index))
    concat_df = pd.concat([below_df, above_df])

    sequences_dict = create_sequence_library({"PR8": concat_df})

    compare_nucleotide_occurrence(sequences_dict["PR8"], "PR8")
    compare_direct_repeats(sequences_dict["PR8"], "PR8", mode=1, correction=True)
    compare_direct_repeats(sequences_dict["PR8"], "PR8", mode=1, correction=False)
