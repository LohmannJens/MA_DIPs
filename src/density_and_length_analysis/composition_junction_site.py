'''
    Loads the data for the deletion sides (start/end point) from Alnaji 2019.
    Does different analysis with the data.

    1.
    Takes a look at the nucleotide distribution around these points.
    Goal is to see if any of the four nucleotides occure more/less often.

    2.
    Compares the nucleotides before the junction start with the ones before
    the junction end site. Counts the number of nucleotides, that are the same.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats

sys.path.insert(0, "..")
from utils import RESULTSPATH, SEGMENTS, COLORS, NUCLEOTIDES, QUANT, S_ROUNDS
from utils import load_alnaji_excel, load_short_reads, get_sequence, get_stat_symbol, generate_sampling_data, create_sequence_library


def count_nucleotide_occurrence(seq: str, p: int)-> dict:
    '''
        Counts the number of nucleotides next to a given point.
        Goes 5 steps in both directions.
        :param seq: whole RNA sequence
        :param p: point on the sequence where to count

        :return: returns a counter dict with an entry for each nucleotide. In
                 each entry the counter for each position is given.
    '''
    window = seq[p-5:p+4]
    r_dict = dict({n: np.zeros(9) for n in NUCLEOTIDES})

    for i, char in enumerate(window):
        r_dict[char][i] = 1
    return r_dict


def count_nucleotide_occurrence_overall(df: object, seq: str)-> (dict, dict):
    '''
        Counts the occurrence of each nucleotide at different positions around
        the junction site
        :param df: dataframe with sequence and junction site data
        :param seq: rna sequence where to count the occurrence

        :return: tupel with two entries:
                    dict with nucleotide count for start site
                    dict with nucleotide count for end site
    '''

    count_start_dict = dict({n: np.zeros(9) for n in NUCLEOTIDES})
    count_end_dict = dict({n: np.zeros(9) for n in NUCLEOTIDES})
    normalize = 0

    for i, row in df.iterrows():
        seq_start_dict = count_nucleotide_occurrence(seq, row["Start"]) 
        seq_end_dict = count_nucleotide_occurrence(seq, row["End"])
        normalize += 1
        for nuc in count_start_dict.keys():
            count_start_dict[nuc] += seq_start_dict[nuc]
            count_end_dict[nuc] += seq_end_dict[nuc]

    return count_start_dict, count_end_dict


def nucleotide_occurrence_analysis(seq_dict: dict, seg: str)-> None:
    '''
        gets the sequences for all four strains and calculates the occurrence
        of each nucleotide at the start and end deletion site.
        :param seq_dict: dictionary with the sequences
        :param seg: name of the segment taht is analyzed

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
        n_sampling = n * S_ROUNDS

        sampling_data = generate_sampling_data(seq, s, e, n_sampling)
        exp_s, exp_e = count_nucleotide_occurrence_overall(sampling_data, seq)
        
        fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)
        x = np.arange(0.8, 9.8, dtype=np.float64)

        for idx, nuc in enumerate(count_start_dict.keys()):
            h_s = count_start_dict[nuc]/n
            h_e = count_end_dict[nuc]/n
            y_exp_s = exp_s[nuc] / (n_sampling)
            y_exp_e = exp_e[nuc] / (n_sampling)

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
                axs[idx, i].legend()
                axs[idx, i].margins(x=0)
                axs[idx, i].set_xlim(left=0.5, right=9.5)
                axs[idx, i].set_ylim(top=0.8, bottom=0.0)
                axs[idx, i].set_xticks([1,2,3,4,5,6,7,8,9])
                axs[idx, i].set_xlabel("position at junction side")
                axs[idx, i].set_ylabel("relative occurrence")
            axs[idx, 0].add_patch(plt.Rectangle((5.5, 0), 4, 1, color="grey", alpha=0.3))
            axs[idx, 1].add_patch(plt.Rectangle((0.5, 0), 4, 1, color="grey", alpha=0.3))
  
        plt.suptitle(f"start (left) and end (right) of {seg} of {k} ({n})")
        savepath = os.path.join(RESULTSPATH, "relative_occurrence_nucleotides", f"{k}_{seg}.png")
        plt.savefig(savepath)
        plt.close()


def calculate_overlapping_nucleotides(seq: str, s: int, e: int, w_len: int, m: int)-> (int, str):
    '''
        counts the number of overlapping nucleotides directly before start and
        end of junction site.
        :param seq: nucleotide sequence
        :param w_len: length of window to be searched
        :param s: start point
        :param e: end point
        :param m: mode how the overlap is created
                    full junction site is 1234J6789
                    1: S 1234J | E 1234J (overlap beginning of both sites)
                    2: S 1234J | E J6789 (overlap what stays in deletion seq)

        :return: Tuple with two entries:
                    Integer giving the number of overlapping nucleotides
                    String of the overlapping nucleotides
    '''
    counter = 0

    if m == 1:
        start_window = seq[s-w_len: s]
        end_window = seq[e-1-w_len: e-1]
        for i in range(w_len - 1, -1, -1):
            if start_window[i] == end_window[i]:
                counter += 1
            else:
                break
        overlap_seq = str(start_window[i+1:w_len])

    elif m == 2:
        for i in range(w_len):
            if seq[s-i:s] == seq[e-1:e-1+i]:
                counter = i
                overlap_seq = str(seq[s-i:s])

    assert counter == len(overlap_seq), f"{counter=}, {len(overlap_seq)}"
    if len(overlap_seq) == 0:
        overlap_seq = "_"

    return counter, overlap_seq


def count_overlapping_nucleotides_overall(df: object, seq: str ,mode: int)-> (dict, dict):
    '''
        calculates the number of overlapping nucleotides directly before start
        and end of junction site for each data point.
        :param df: dataframe with sequence and junction site data
        :param seq: RNA sequence of the given segement and strain
        :param mode: states which calculation mode is used in 
                     calculate_overlapping_nucleotides() check there for info

        :return: Tuple including a dict with the count of the length of
                 overlapping sequences and a dict with the overlapping
                 sequences and their count.
    '''
    w_len = 15
    nuc_overlap_dict = dict({i: 0 for i in range(0, w_len+1)})
    overlap_seq_dict = dict()
 
    for i, row in df.iterrows():
        s = row["Start"]
        e = row["End"]
        idx, overlap_seq = calculate_overlapping_nucleotides(seq, s, e, w_len, mode)
        nuc_overlap_dict[idx] += 1
        if overlap_seq in overlap_seq_dict:
            overlap_seq_dict[overlap_seq] += 1
        else:
            overlap_seq_dict[overlap_seq] = 1

    return nuc_overlap_dict, overlap_seq_dict


def count_nucleotide_freq_overlap_seq(seq_dict: dict, strain: str, seg: str)-> object:
    '''
        Counts the number of nucleotides of the overlapping sequences and
        normalizes them. Also gets the expected values, which is the occurrence
        of each nucleotide over the whole seuqence.
        :param seq_dict: dict with overlapping sequences and number of counts
        :param strain: name of strain that is analyzed
        :param seg: name of the segment that is analyzed
    
        :return: Data Frame with one column for each nucleotide
    '''
    w_len = 15
    count_dict = dict({n: np.zeros(w_len) for n in NUCLEOTIDES})
    for k, v in seq_dict.items():
        seq = k
        seq_len = len(seq)
        for n in NUCLEOTIDES:
            count_dict[n][seq_len] += seq.count(n) * v
    df = pd.DataFrame(count_dict)

    full_seq = get_sequence(strain, seg)
    labels = dict()
    for n in NUCLEOTIDES:
        labels[n] = full_seq.count(n)
    exp_df = pd.DataFrame(labels, index=["exp"])
    df.iloc[0] = exp_df.iloc[0]
    sum_df = df.sum(axis=1).astype(int)
    n_df = df.div(df.sum(axis=1), axis=0)
    n_df["Sum"] = sum_df
    n_df = n_df.fillna(0)
    n_df["label"] = n_df.index
    n_df["label"] = n_df["label"].replace([0], "exp").astype(str)
    return n_df


def nucleotide_overlap_analysis(seq_dict: dict, seg: str, mode: int, ngs_thresh: int=0)-> None:
    '''
        gets the sequences for all four strains and calculates the overlap of 
        the nucleotides at the junction site. Also generates the overlapping
        sequences and plots them.
        :param seq_dict: dictionary with the sequences
        :param seg: name of the segment that is analyzed
        :param mode: states which calculation mode is used in 
                     calculate_overlapping_nucleotides() check there for info
        :param ngs_thresh: gives the threshold on which data to include
        :return: None
    '''
    fig, axs = plt.subplots(4, 2, figsize=(10, 10), tight_layout=True,
                            gridspec_kw={"width_ratios": [1, 3]})
    for i, (k, v) in enumerate(seq_dict.items()):
        v = v.loc[(v["Segment"] == seg) & (v["NGS_read_count"] > ngs_thresh)]
        seq = get_sequence(k, seg)
        nuc_overlap_dict, overlap_seq_dict = count_overlapping_nucleotides_overall(v, seq, mode)
        n = len(v.index)
        # only plot results if at least one data point
        if n == 0:
            continue
        
        # get expected values
        s = (int(v.Start.quantile(QUANT)), int(v.Start.quantile(1-QUANT)))
        e = (int(v.End.quantile(QUANT)), int(v.End.quantile(1-QUANT)))
        n_sampling = n * S_ROUNDS

        sampling_data = generate_sampling_data(seq, s, e, n_sampling)
        exp, _ = count_overlapping_nucleotides_overall(sampling_data, seq, mode)
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

        axs[i, 0].bar(x=x, height=h/h.sum(), width=-0.4, align="edge", label="observed")
        axs[i, 0].bar(x=x, height=h_exp/h_exp.sum(), width=0.4, align="edge", label="expected")
        axs[i, 0].set_xlabel("number of overlapping nucleotides")
        axs[i, 0].set_ylabel("relative occurrence")
        axs[i, 0].set_title(f"{k} (n={n}) {symbol}")
        axs[i, 0].legend(loc="upper right")
        axs[i, 0].set_ylim(bottom=0.0, top=1.0)

        plot_df = count_nucleotide_freq_overlap_seq(overlap_seq_dict, k, seg)
        bottom = np.zeros(15)
        for n in NUCLEOTIDES:
            axs[i, 1].bar(x=plot_df["label"], height=plot_df[n], bottom=bottom, label=n, color=COLORS[n])
            bottom += plot_df[n]

        f_exp = plot_df.iloc[0][NUCLEOTIDES]
        for pos in plot_df["label"]:
            symbol = ""
            if pos == "exp":
                pos = 0
            else:
                f_obs = plot_df.iloc[int(pos)][NUCLEOTIDES]
                if f_obs.sum() != 0:
                    res = stats.chisquare(f_obs, f_exp)
                    symbol = get_stat_symbol(res.pvalue)
            text = f"n={str(int(plot_df['Sum'].iloc[int(pos)]))}\n{symbol}"
            axs[i, 1].annotate(text, (pos, 0.0), fontsize="xx-small", ha="center")

        axs[i, 1].legend()
        axs[i, 1].set_title("nucleotide occurrence in overlapping sequence")
        axs[i, 1].set_xlabel("length of overlapping sequence")
        axs[i, 1].set_ylabel("relative nuc. occurrence")

    ngs_thresh = "" if ngs_thresh == 0 else f"NGS{ngs_thresh}_"
    fname = f"{seg}_mode{mode}_{ngs_thresh}sequence_distribution.pdf"
    savepath = os.path.join(RESULTSPATH, "overlapping_nucleotides", fname)
    plt.savefig(savepath)
    plt.close()


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)

    # Create a sequence library for each strain
    sequence_list_dict = create_sequence_library(all_reads_dict)

    # Loop over the different strains and calculate the occurrence of each
    # nucleotide in the sequences
    for s in SEGMENTS:
        nucleotide_occurrence_analysis(sequence_list_dict, s)

    # Check if nucleotides directly before junction site have the same sequence
    # as the ones directly before junction site at the end
    for s in SEGMENTS:
        nucleotide_overlap_analysis(sequence_list_dict, s, mode=1)
        nucleotide_overlap_analysis(sequence_list_dict, s, mode=2)

#        nucleotide_overlap_analysis(sequence_list_dict, s, mode=1, ngs_thresh=1000)
 #       nucleotide_overlap_analysis(sequence_list_dict, s, mode=2, ngs_thresh=1000)

