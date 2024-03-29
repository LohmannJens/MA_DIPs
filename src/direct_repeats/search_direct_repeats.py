'''
    Compares the nucleotides before the junction start with the ones before
    the junction end site. Counts the number of nucleotides, that are the same.
    --> direct repeats
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


def calculate_direct_repeat(seq: str,
                            s: int,
                            e: int,
                            w_len: int,
                            m: int
                            )-> Tuple[int, str]:
    '''
        Counts the number of overlapping nucleotides directly before start and
        end of junction site --> direct repeats
        :param seq: nucleotide sequence
        :param w_len: length of window to be searched
        :param s: start point
        :param e: end point
        :param m: mode how the overlap is created
                    full junction site is 1234J6789
                    1: S 1234J | E 1234J (overlap beginning of both sites)
                    2: S 1234J | E J6789 (overlap what stays in deletion seq)
                    3: S J6789 | E J6789 (overlap in DG-seq (Boussier))

        :return: Tuple with two entries:
                    Integer giving the number of overlapping nucleotides
                    String of the overlapping nucleotides
    '''
    counter = 0

    if m == 1:
        start_window = seq[s-w_len: s]
        end_window = seq[e-1-w_len: e-1]
        
        #if they are the same return directly to avoid off-by-one error
        if start_window == end_window:
            return len(start_window), start_window

        if len(seq) < e:
            return 0, "_"
        
        for i in range(len(end_window) - 1, -1, -1):
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

    elif m == 3:
        start_window = seq[s: s+w_len]
        end_window = seq[e: e+w_len]
        
        #if they are the same return directly to avoid off-by-one error
        if start_window == end_window:
            return len(start_window), start_window

        if len(seq) < e:
            return 0, "_"
     
        for i in range(0, len(end_window)):
            if start_window[i] == end_window[i]:
                counter += 1
            else:
                break
        overlap_seq = str(start_window[:counter])

    assert counter == len(overlap_seq), f"{counter=}, {len(overlap_seq)}"
    if len(overlap_seq) == 0:
        overlap_seq = "_"

    return counter, overlap_seq


def count_direct_repeats_overall(df: pd.DataFrame,
                                 seq: str,
                                 mode: int
                                 )-> Tuple[dict, dict]:
    '''
        Calculates the number of direct repeats for each data point.
        :param df: dataframe with sequence and junction site data
        :param seq: RNA sequence of the given segement and strain
        :param mode: states which calculation mode is used in 
                     calculate_overlapping_nucleotides() check there for info

        :return: Tuple including a dict with the count of the length of
                 overlapping sequences and a dict with the overlapping
                 sequences and their count.
    '''
    w_len = 5
    nuc_overlap_dict = dict({i: 0 for i in range(0, w_len+1)})
    overlap_seq_dict = dict()
 
    for i, row in df.iterrows():
        s = row["Start"]
        e = row["End"]
        idx, overlap_seq = calculate_direct_repeat(seq, s, e, w_len, mode)
        nuc_overlap_dict[idx] += 1
        if overlap_seq in overlap_seq_dict:
            overlap_seq_dict[overlap_seq] += 1
        else:
            overlap_seq_dict[overlap_seq] = 1

    return nuc_overlap_dict, overlap_seq_dict


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
        orig_value = nuc_overlap_dict[idx]
        if orig_value != 0:
            divided_value = orig_value/(idx+1)
            new[idx] = divided_value
            for idx_2 in range(0, idx):
                new[idx_2] = new[idx_2] + divided_value
        else:
            new[idx] = 0

    return new


def direct_repeats_analysis(seq_dict: dict,
                            mode: int,
                            top: bool=False,
                            correction: bool=False,
                            author: str=""
                            )-> None:
    '''
        Calculates the direct repeats of all sequences of the Pelz dataset.
        :param seq_dict: dictionary with the sequences
        :param mode: states which calculation mode is used. Check
                     composition_junction_site.py for more info.
        :param top: states if the whole dataset or just the top DI RNA are used
        :param correction: if True a correction calculation is made
        :param author: name of author data published the given data
        :return: None
    '''
    plt.rc("font", size=12)
    for k, v in seq_dict.items():
        fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)
        j = 0
        for i, s in enumerate(SEGMENTS):
            v_s = v.loc[(v["Segment"] == s)]
            seq = get_sequence(k, s)
            nuc_overlap_dict, _ = count_direct_repeats_overall(v_s, seq, mode)  
            n = len(v_s.index)
            if n <= 1:
                if k != "NC":
                    axs[i%4, j].set_visible(False)
                continue

            if correction:
                nuc_overlap_dict = include_correction(nuc_overlap_dict)

            start = (int(v_s.Start.quantile(QUANT)), int(v_s.Start.quantile(1-QUANT)))
            end = (int(v_s.End.quantile(QUANT)), int(v_s.End.quantile(1-QUANT)))
            sampling_data = generate_sampling_data(seq, start, end, N_SAMPLES)
            exp, _ = count_direct_repeats_overall(sampling_data, seq, mode)

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

            stat, pvalue = stats.ks_2samp(f_obs, f_exp)
            symbol=get_stat_symbol(pvalue)

            axs[i%4, j].bar(x=x, height=h/h.sum(), width=-0.4, align="edge", label="observed")
            axs[i%4, j].bar(x=x, height=h_exp/h_exp.sum(), width=0.4, align="edge", label="expected")
            axs[i%4, j].set_xlabel("length of direct repeat")
            axs[i%4, j].set_ylabel("relative occurrence")
            axs[i%4, j].set_title(f"{s} (n={n}) {symbol}")
            axs[i%4, j].set_xticks([0, 5, 10, 15])
            axs[i%4, j].set_yticks([1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
            axs[i%4, j].legend(loc="upper right")
            axs[i%4, j].set_ylim(bottom=0.0, top=1.0)

            if i == 3:
                j = 1

        fname = f"{k}_mode{mode}"
        if author != "":
            fname = f"{author}_{fname}"
        if top:
            fname = f"{fname}_top"
        if correction:
            fname = f"{fname}_corr"
        fname = f"{fname}.png"

        savepath = os.path.join(RESULTSPATH, "direct_repeats", fname)
        plt.savefig(savepath)
        plt.close()


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)

    # Create a sequence library for each strain
    sequence_list_dict = create_sequence_library(all_reads_dict)

    # Check for direct repeats
    direct_repeats_analysis(sequence_list_dict, mode=1)
    direct_repeats_analysis(sequence_list_dict, mode=1, correction=True)
    direct_repeats_analysis(sequence_list_dict, mode=2)

