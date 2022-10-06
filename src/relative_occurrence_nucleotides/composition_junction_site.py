'''
    Loads the data for the deletion sides (start/end point) from Alnaji 2019.
    Does different analysis with the data.

    1.
    Takes a look at the nucleotide distribution around these points.
    Goal is to see if any of the four nucleotides occure more/less often.

    2.
    Compares the nucleotides before the junction start with the ones before
    the junction end site. Counts the number of nucleotides, that are the same.
    --> direct repeats
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats
from decimal import Decimal, ROUND_HALF_UP

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


def count_nuc_freq_overlap(seq_dict: dict, strain: str, seg: str)-> object:
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
    # count nucleotide occurrence split up by overlapping sequence length
    count_dict = dict({n: 0 for n in NUCLEOTIDES})
    for k, v in seq_dict.items():
        seq = k
        for n in NUCLEOTIDES:
            count_dict[n] += seq.count(n) * v
    df = pd.DataFrame(count_dict, index=[0])

    # create expected data by counting over the whole sequence
    full_seq = get_sequence(strain, seg)
    labels = dict()
    for n in NUCLEOTIDES:
        labels[n] = full_seq.count(n)
    exp_df = pd.DataFrame(labels, index=[0])

    df = pd.concat([df, exp_df])
    # Add sum column and normalize data
    sum_df = df.sum(axis=1).astype(int)
    n_df = df.div(df.sum(axis=1), axis=0)
    n_df["Sum"] = sum_df
    n_df["label"] = [f"obs", f"exp"]
    return n_df


def include_correction(nuc_overlap_dict: dict)-> dict:
    '''
        Adds a correction to the counting of the direct repeats. This is due to
        the fact that at these sites the observations get merged towards higher
        lengths of the direct repeat.
        :param nuc_overlap_dict: counting dict of the direct repeats

        :return: corrected counts of the direct repeats (same structure as
                 input
    '''
    new = dict()
    for idx in nuc_overlap_dict.keys():
        org_value = nuc_overlap_dict[idx]
        if org_value != 0:
            divided_value = org_value/(idx+1)
            new[idx] = divided_value
            for idx_2 in range(0, idx):
                new[idx_2] = new[idx_2] + divided_value
        else:
            new[idx] = 0

    return new


def nucleotide_overlap_analysis(seq_dict: dict, seg: str, mode: int, ngs_thresh: int=0, correction: bool=False)-> None:
    '''
        gets the sequences for all four strains and calculates the overlap of 
        the nucleotides at the junction site.
        :param seq_dict: dictionary with the sequences
        :param seg: name of the segment that is analyzed
        :param mode: states which calculation mode is used in 
                     calculate_overlapping_nucleotides() check there for info
        :param ngs_thresh: gives the threshold on which data to include
        :param correction: if True a correction calculation is made

        :return: None
    '''
    overlap_seq_per_strain = dict()
    fig, axs = plt.subplots(4, 1, figsize=(5, 10), tight_layout=True)
    for i, (k, v) in enumerate(seq_dict.items()):
        v = v.loc[(v["Segment"] == seg) & (v["NGS_read_count"] > ngs_thresh)]
        seq = get_sequence(k, seg)
        nuc_overlap_dict, overlap_seq_dict = count_overlapping_nucleotides_overall(v, seq, mode)
        overlap_seq_per_strain[k] = overlap_seq_dict
        n = len(v.index)
        # only plot results if at least one data point
        if n == 0:
            continue
        
        if correction:
            nuc_overlap_dict = include_correction(nuc_overlap_dict)

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
            if correction:
                f_obs.extend([a]*int(Decimal(nuc_overlap_dict[a]).to_integral_value(rounding=ROUND_HALF_UP)))
            else:
                f_obs.extend([a]*nuc_overlap_dict[a])
            f_exp.extend([a]*exp[a])
        f_obs = np.array(f_obs)
        f_exp = np.array(f_exp)

        # select statistical test here
        #stats_test = "mannwhitneyu"
        stats_test = "ks_2samp"

        if stats_test == "mannwhitneyu":
           res = stats.mannwhitneyu(f_obs, f_exp)
           symbol = get_stat_symbol(res.pvalue)
        elif stats_test == "ks_2samp":
            stat, pvalue = stats.ks_2samp(f_obs, f_exp)
            symbol=get_stat_symbol(pvalue)

        # plot results as barplot
        axs[i].bar(x=x, height=h/h.sum(), width=-0.4, align="edge", label="observed")
        axs[i].bar(x=x, height=h_exp/h_exp.sum(), width=0.4, align="edge", label="expected")
        axs[i].set_xlabel("number of overlapping nucleotides")
        axs[i].set_ylabel("relative occurrence")
        axs[i].set_title(f"{k} (n={n}) {symbol}")
        axs[i].legend(loc="upper right")
        axs[i].set_ylim(bottom=0.0, top=1.0)

    ngs_thresh = "" if ngs_thresh == 0 else f"NGS{ngs_thresh}_"
    corr = "_corr" if correction else ""
    fname = f"{seg}_mode{mode}_{ngs_thresh}sequence_distribution{corr}.pdf"
    savepath = os.path.join(RESULTSPATH, "overlapping_nucleotides", fname)
    plt.savefig(savepath)
    plt.close()

    return overlap_seq_per_strain


def overlap_seq_occurrence_analysis(overlap_dict: dict)-> None:
    '''
        Gets a dictionary with all overlapping sequences, split by segments and
        strains. Calculates the relative occurrence of each nucleotide and
        returns the results as a barplot.
        :param overlap_dict: dictionary with the overlapping sequences

        :return: None
    '''
    strains = list(overlap_dict["PB1"].keys())
    fig, axs = plt.subplots(4, 1, figsize=(10, 10), tight_layout=True)
    for i, st in enumerate(strains):
        f_df = pd.DataFrame()
        for s in SEGMENTS:
            # count the occurrence of the four nucleotides in the overlapping sequences
            f_df = pd.concat([f_df, count_nuc_freq_overlap(overlap_dict[s][st], st, s)])
            
        # plot the results as stacked barplot
        x = np.arange(0, 8)
        b_obs = np.zeros(8)
        b_exp = np.zeros(8)
        for n in NUCLEOTIDES:
            obs = f_df[f_df["label"] == "obs"]
            exp = f_df[f_df["label"] == "exp"]
            axs[i].bar(x-0.2, height=obs[n], width=0.3, bottom=b_obs, label=n, color=COLORS[n])
            axs[i].bar(x+0.2, height=exp[n], width=0.3, bottom=b_exp, color=COLORS[n], alpha=0.5)
            b_obs += obs[n]
            b_exp += exp[n]

        # test nucleotide distribution in overlapping sequence for significance
        for pos in range(0, 15, 2):
            f_exp = f_df.iloc[pos + 1][NUCLEOTIDES]
            f_obs = f_df.iloc[pos][NUCLEOTIDES]

            symbol = ""
            if f_obs.sum() != 0:
                res = stats.chisquare(f_obs, f_exp)
                symbol = get_stat_symbol(res.pvalue)
            text = f"n={str(int(f_df['Sum'].iloc[pos]))}\n{symbol}"
            axs[i].annotate(text, (pos/2-0.2, 0.0), fontsize="xx-small", ha="center")

        axs[i].set_xlim(left=-0.5, right=8)
        axs[i].set_xticks(np.arange(0, 8), SEGMENTS)
        axs[i].legend(loc="upper right", fontsize="small")
        axs[i].set_title(f"nucleotide occurrence in overlapping sequence for {st}")
        axs[i].set_xlabel("length of overlapping sequence")
        axs[i].set_ylabel("rel. nuc. occurrence")

    fname = f"nuc_dist_overlap_sequence.pdf"
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
    overlap_seq_per_seg = dict()
    for s in SEGMENTS:
        overlap_per_strain = nucleotide_overlap_analysis(sequence_list_dict, s, mode=1)
        _ = nucleotide_overlap_analysis(sequence_list_dict, s, mode=1, correction=True)
        _ = nucleotide_overlap_analysis(sequence_list_dict, s, mode=2)
        overlap_seq_per_seg[s] = overlap_per_strain

    overlap_seq_occurrence_analysis(overlap_seq_per_seg)

