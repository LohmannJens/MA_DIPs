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
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio.Seq import Seq

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, load_excel, load_short_reads, get_sequence


COLORS = dict({"A": "blue", "C": "orange", "G": "green", "U": "red"})
NUCLEOTIDES = list(["A", "C", "G", "U"])


def create_sequence(s: int, e: int, strain: str, seg: str, crop: bool = False)-> str:
    '''
        Loads a DNA sequence and removes a part of it. The deletion is given by
        a start and an end. Those positions are part of the remaining sequence.
        :param s: start of the deletion (included in the remaining sequence)
        :param e: end of the deletion (included in the remaining sequence)
        :param strain: name of the virus strain
        :param seg: segment where the sequenz is coming from
        :param crop: Indicates if the sequence should be cropped (True) to
                     exclude the deletion side or take whole sequence (False)

        :return: RNA sequence with or without deletion side
    '''
    full_dna_record = get_sequence(strain, seg)
    full_dna_seq = full_dna_record.seq
    full_rna_seq = full_dna_seq.transcribe()

    if crop:
        return full_rna_seq[:s] + full_rna_seq[e-1:]
    else:
        return full_rna_seq


def create_sequence_library(data_dict: dict)-> dict:
    '''
        gets the raw loaded sequence data, which is a dict over all strains.
        In each dict the value is a data frame with the rows and columns from
        the loaded excel file.
        Creates the deletion sequence and saves it with other features as 
        sequence length, ... in a pandas data frame.
        :param data_dict: dictionary of the loaded excel

        :return: dictionary with key for each strain. Value is a pandas df.
    '''
    for k, v in data_dict.items():
        del_sequence_list = list()
        whole_sequence_list = list()
        # loop over rows of value dataframe
        for i, row in v.iterrows():
            start = row["Start"]
            end = row["End"]
            segment = row["Segment"]
            del_sequence = create_sequence(start, end, k, segment, crop=True)
            whole_sequence = create_sequence(start, end, k, segment, crop=False)
            del_sequence_list.append(del_sequence)
            whole_sequence_list.append(whole_sequence)
        
        data_dict[k]["DelSequence"] = del_sequence_list
        data_dict[k]["WholeSequence"] = whole_sequence_list

    return data_dict


def count_nucleotide_occurrence_overall(df):
    '''
        Counts the occurrence of each nucleotide at different positions around
        the junction site
        :param df: dataframe with sequence and junction site data

        :return: tupel with three entries:
                    dict with nucleotide count for start site
                    dict with nucleotide count for end site
    '''
    def count_nucleotide_occurrence(seq: str, p: int)-> dict: 
        window = seq[p-4:p+5]
        r_dict = dict({n: np.zeros(9) for n in NUCLEOTIDES})
 
        for i, char in enumerate(window):
            r_dict[char][i] = 1
        return r_dict

    count_start_dict = dict({n: np.zeros(9) for n in NUCLEOTIDES})
    count_end_dict = dict({n: np.zeros(9) for n in NUCLEOTIDES})
    normalize = 0

    for i, row in df.iterrows():
        seq_start_dict = count_nucleotide_occurrence(row["WholeSequence"], row["Start"]) 
        seq_end_dict = count_nucleotide_occurrence(row["WholeSequence"], row["End"])
        normalize += 1
        for nuc in count_start_dict.keys():
            count_start_dict[nuc] += seq_start_dict[nuc]
            count_end_dict[nuc] += seq_end_dict[nuc]

    return count_start_dict, count_end_dict


def count_overlapping_nucleotides_overall(df)-> (dict, dict):
    '''
        calculates the number of overlapping nucleotides directly before start
        and end of junction site for each data point.
        :param df: dataframe with sequence and junction site data

        :return: Tuple including a dict with the count of the length of
                 overlapping sequences and a dict with the overlapping
                 sequences and their count.
    '''
    def calculate_overlapping_nucleotides(seq: str, s: int, e: int)-> (int, str):
        window_len = 10
        start_window = seq[s-window_len: s]
        end_window = seq[e-1-window_len: e-1]
        counter = 0

        for i in range(window_len - 1, -1, -1):
            if start_window[i] == end_window[i]:
                counter += 1
            else:
                break
        return counter, str(start_window[i:window_len])

    nuc_overlap_dict = dict({i: 0 for i in range(0, 11)})
    overlap_seq_dict = dict()
 
    for i, row in df.iterrows():
        sequence = row["WholeSequence"]
        idx, overlap_seq = calculate_overlapping_nucleotides(sequence, row["Start"], row["End"])
        nuc_overlap_dict[idx] += 1
        if overlap_seq in overlap_seq_dict:
            overlap_seq_dict[overlap_seq] += 1
        else:
            overlap_seq_dict[overlap_seq] = 1

    return nuc_overlap_dict, overlap_seq_dict


def generate_sampling_data(seq: str, s: (int, int), e: (int, int),  n: int):
    '''
        generates sampling data by creating random start and end points for
        artificial junction sites. Generated data is used to calculate the
        expected values. Sample set is 3 times the size of the observation set.
        :param seq: sequence of the segment
        :param s: tuple with start and end point of the range for the artifical
                  start point of the junction
        :param e: tuple with start and end point of the range for the artifical
                  end point of the junction
        :param n: size of the observation data set

        :return: dataframe with the artifical data set
    '''
    sampling = dict({"WholeSequence": [], "Start": [], "End": []})
    for _ in range(n*3):
        sampling["WholeSequence"].append(seq)
        sampling["Start"].append(random.randint(s[0], s[1]))
        sampling["End"].append(random.randint(e[0], e[1]))
    return pd.DataFrame(data=sampling)


def nucleotide_occurrence_analysis(seq_dict: dict, seg: str)-> None:
    '''
        gets the sequences for all four strains and calculates the occurrence
        of each nucleotide at the start and end deletion site.
        :param seq_dict: dictionary with the sequences
        :param seg: name of the segment that is analyzed; 'all' if all 8 are used

        :return: None
    '''
    for k, v in seq_dict.items():
        v = v.loc[v["Segment"] == seg]
        
        count_start_dict, count_end_dict = count_nucleotide_occurrence_overall(v)
        n = len(v.index)

        # only plot results if more then 10 data points
        if n <= 10:
            continue

        # get expected values
        seq = v.iloc[0]["WholeSequence"]
        q = 0.20
        s = (int(v.Start.quantile(q)), int(v.Start.quantile(1-q)))
        e = (int(v.End.quantile(q)), int(v.End.quantile(1-q)))
        sampling_data = generate_sampling_data(seq, s, e, n)
        expected_start, expected_end = count_nucleotide_occurrence_overall(sampling_data)

        fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)
        x = np.arange(0.7, 9.7, dtype=np.float64)

        for idx, nuc in enumerate(count_start_dict.keys()):
            h_s = count_start_dict[nuc]/n
            h_e = count_end_dict[nuc]/n
            h_s_exp = height=expected_start[nuc]/(n*3)
            h_e_exp = height=expected_end[nuc]/(n*3)
            axs[idx, 0].bar(x, height=h_s, width=0.3, label=nuc, color=COLORS[nuc])
            axs[idx, 1].bar(x, height=h_e, width=0.3, label=nuc, color=COLORS[nuc])
            axs[idx, 0].bar(x+0.4, height=h_s_exp, width=0.3, label=f"{nuc}_exp", color=COLORS[nuc], alpha=0.5)
            axs[idx, 1].bar(x+0.4, height=h_e_exp, width=0.3, label=f"{nuc}_exp", color=COLORS[nuc], alpha=0.5)

            for i in range(2):
                axs[idx, i].legend()
                axs[idx, i].margins(x=0)
                axs[idx, i].set_xlim(left=0.5)
                axs[idx, i].set_ylim(top=0.8)
                axs[idx, i].set_xticks([1,2,3,4,5,6,7,8,9])
                axs[idx, i].set_xlabel("position at junction side")
                axs[idx, i].set_ylabel("relative occurrence")
            axs[idx, 0].add_patch(plt.Rectangle((5.5, 0), 4, 1, color="grey", alpha=0.3))
            axs[idx, 1].add_patch(plt.Rectangle((0.5, 0), 4, 1, color="grey", alpha=0.3))
  
        plt.suptitle(f"start (left) and end (right) of deletion site of {seg} of {k} ({n})")
        savepath = os.path.join(RESULTSPATH, "relative_occurrence_nucleotides", f"{k}_{seg}.pdf")
        plt.savefig(savepath)
        plt.close()


def nucleotide_overlap_analysis(seq_dict: dict, seg: str)-> None:
    '''
        gets the sequences for all four strains and calculates the overlap of 
        the nucleotides at the junction site. Also generates the overlapping
        sequences and plots them.
        :param seq_dict: dictionary with the sequences
        :param seg: name of the segment that is analyzed

        :return: None
    '''
    fig, axs = plt.subplots(4, 2, figsize=(10, 10), tight_layout=True, gridspec_kw={"width_ratios": [1, 3]})
    for i, (k, v) in enumerate(seq_dict.items()):
        v = v.loc[v["Segment"] == seg]
        
        nuc_overlap_dict, overlap_seq_dict = count_overlapping_nucleotides_overall(v)
        n = len(v.index)
        # only plot results if more then 10 data points
        if n <= 10:
            continue
        
        # get expected values
        seq = v.iloc[0]["WholeSequence"]
        q = 0.20
        s = (int(v.Start.quantile(q)), int(v.Start.quantile(1-q)))
        e = (int(v.End.quantile(q)), int(v.End.quantile(1-q)))

        sampling_data = generate_sampling_data(seq, s, e, n)
        expected, _ = count_overlapping_nucleotides_overall(sampling_data)

        x = list(nuc_overlap_dict.keys())
        h = np.array(list(nuc_overlap_dict.values()))
        x_exp = list(expected.keys())
        y_exp = np.array(list(expected.values()))
        
        axs[i, 0].bar(x=x, height=h/h.sum(), width=-0.4, align="edge", label="observed")
        axs[i, 0].bar(x=x_exp, height=y_exp/y_exp.sum(), width=0.4, align="edge", label="expected")

        axs[i, 0].set_xlabel("numbaer of overlapping nucleotides")
        axs[i, 0].set_ylabel("relative occurrence")
        axs[i, 0].set_title(f"{k}")
        axs[i, 0].legend("upper left")

        plot_dict = dict()
        for key, value in overlap_seq_dict.items():
            if value >= 3:
                plot_dict[key] = value

        # order elements by decreasing value
        plot_dict = {k: v for k, v in sorted(plot_dict.items(), key=lambda item: item[1], reverse=True)}

        x = list(plot_dict.keys())
        h = list(plot_dict.values())
        axs[i, 1].bar(x=x, height=h, width=0.5)

        axs[i, 1].set_xlabel("nucleotide sequence")
        axs[i, 1].set_ylabel("absolute occurrence")
        axs[i, 1].set_title(f"occurrence of sequences that overlap at start and end")

    savepath = os.path.join(RESULTSPATH, "overlapping_nucleotides", f"{seg}_sequence_distribution.pdf")
    plt.savefig(savepath)
    plt.close()


if __name__ == "__main__":
    filepath = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    cleaned_data_dict = load_excel(filepath)

    short_reads_filepath = os.path.join(DATAPATH, "alnaji2019", "Small_deletionSize_FA.xlsx")
    all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)

    # Create a sequence library for each strain
    sequence_list_dict = create_sequence_library(all_reads_dict)

    # Loop over the different strains and calculate the occurrence of each
    # nucleotide in the sequences
    for s in SEGMENTS:
        nucleotide_occurrence_analysis(sequence_list_dict, s)

    # Check if nucleotides directly before junction site have the same sequence
    # as the ones directly before junction site at the end
    for s in SEGMENTS:
        nucleotide_overlap_analysis(sequence_list_dict, s)
