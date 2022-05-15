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
import argparse

import numpy as np
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


def count_nucleotide_occurrence(seq: str, p: int)-> dict: 
    '''
        Calculates the occurrence of the four nucleotides around a point p.
        It takes the nucleotide at point p and four nucleotides up- and
        downstream the sequence.
        :param seq: full RNA sequence to count the nucleotide occurrence
        :param p: middle point of the window of size 9 where the occurrence
                  should be counted

        :return: dict with the nucleotide as key and a list of the absolute
                 occurrences as value
    '''
    window = seq[p-4:p+5]

    r_dict = {"A": np.zeros(9),
              "C": np.zeros(9),
              "G": np.zeros(9),
              "U": np.zeros(9)}

    for i, char in enumerate(window):
        r_dict[char][i] = 1

    return r_dict


def calculate_overlapping_nucleotides(seq: str, s: int, e: int)-> (int, str):
    '''
        calculates the number of overlapping nucleotides directly before start
        and end of junction site.
        :param seq: sequence where the overlap will be calculated on
        :param s: start point of the junction site
        :param e: end point of the junction site

        :return: integer indicating the length of the overlapping sequence
        :return: string of the overlapping sequence
    '''
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

def get_expected_by_sampling(seq: str, s: (int, int), e: (int, int), n: int)-> (dict, dict):
    '''
        Does a sampling to generate a expected value. Selects a random start
        and end point to generate a random deletion. Then counts the occurrence
        of nucleotides around this site. This is done 3 * n while n is the
        number of observed deletions in the data set.
        :param seq: sequence of the segment
        :param s: Tupel with the lower and upper bound for the start of the
                  deletion site
        :param e: Tupel with the lower and upper bound for the end of the
                  deletion site
        :param n: number of samples for this strain and segment observerd in
                  the data set

        :return: Tupel of two dicts including the counts for each nucleotide
                 at each position. One for start and one for the end point
    '''

    start_expected = dict({"A": np.zeros(9), "C": np.zeros(9), "G": np.zeros(9), "U": np.zeros(9)})
    end_expected = dict({"A": np.zeros(9), "C": np.zeros(9), "G": np.zeros(9), "U": np.zeros(9)})

    for _ in range(n*3):
        start_count = count_nucleotide_occurrence(seq, random.randint(s[0], s[1]))
        end_count = count_nucleotide_occurrence(seq, random.randint(e[0], e[1]))
        for nuc in NUCLEOTIDES:
            start_expected[nuc] += start_count[nuc]
            end_expected[nuc] += end_count[nuc]

    return start_expected, end_expected


def nucleotide_occurrence_deletion_site(seq_dict: dict, seg: str, weighted: bool)-> None:
    '''
        gets the sequences for all four strains and calculates the occurrence
        of each nucleotide at the start and end deletion site.
        :param seq_dict: 
        :param seg: name of the segment that is analyzed; 'all' if all 8 are used
        :param weighted: True if the NGS count should be included

        :return: None
    '''
    for k, v in seq_dict.items():
        v = v.loc[v["Segment"] == seg]
        count_start_dict = dict({"A": np.zeros(9), "C": np.zeros(9), "G": np.zeros(9), "U": np.zeros(9)})
        count_end_dict = dict({"A": np.zeros(9), "C": np.zeros(9), "G": np.zeros(9), "U": np.zeros(9)})
        normalize = 0
        
        for i, row in v.iterrows():
            seq_start_dict = count_nucleotide_occurrence(row["WholeSequence"], row["Start"]) 
            seq_end_dict = count_nucleotide_occurrence(row["WholeSequence"], row["End"])
            weight = row["NGS_read_count"] if weighted else 1
            normalize += weight
            for nuc in count_start_dict.keys():
                count_start_dict[nuc] += seq_start_dict[nuc] * weight
                count_end_dict[nuc] += seq_end_dict[nuc] * weight
        # only plot results if more than 10 DI RNA samples 
        if normalize <= 10:
            continue

        # get expected values
        seq = v.iloc[1]["WholeSequence"]
        n = len(v.index)
        q = 0.25
        s = (int(v.Start.quantile(0.25)), int(v.Start.quantile(0.75)))
        e = (int(v.End.quantile(0.25)), int(v.End.quantile(0.75)))

        expected_start, expected_end = get_expected_by_sampling(seq, s, e, n)


        fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)
        x = np.arange(0.7, 9.7, dtype=np.float64)
        for idx, nuc in enumerate(count_start_dict.keys()):
            axs[idx, 0].bar(x, height=count_start_dict[nuc]/normalize, width=0.3, label=nuc, color=COLORS[nuc])
            axs[idx, 1].bar(x, height=count_end_dict[nuc]/normalize, width=0.3, label=nuc, color=COLORS[nuc])
            #plot expected values
            axs[idx, 0].bar(x+0.4, height=expected_start[nuc]/(n*3), width=0.3, label=f"{nuc}_exp", color=COLORS[nuc], alpha=0.5)
            axs[idx, 1].bar(x+0.4, height=expected_end[nuc]/(n*3), width=0.3, label=f"{nuc}_exp", color=COLORS[nuc], alpha=0.5)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze NP density and junction site of deletions")
    parser.add_argument("-w", "--weighted", action="store_true")
    args = parser.parse_args()
    weighted = args.weighted


    filepath = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    cleaned_data_dict = load_excel(filepath)

    short_reads_filepath = os.path.join(DATAPATH, "alnaji2019", "Small_deletionSize_FA.xlsx")
    all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)

    # Create a sequence library for each strain
    sequence_list_dict = create_sequence_library(all_reads_dict)

    # Loop over the different strains and calculate the occurrence of each
    # nucleotide in the sequences
    for s in SEGMENTS:
        nucleotide_occurrence_deletion_site(sequence_list_dict, s, weighted)


    # Check if nucleotides directly before junction site have the same sequence
    # as the ones directly before junction site at the end
    overlap_seq_dict = dict()
    for k, v in sequence_list_dict.items():
        nuc_overlap_dict = dict({0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10:0})

        for i, row in v.iterrows():
            sequence = row["WholeSequence"]
            idx, overlap_seq = calculate_overlapping_nucleotides(sequence, row["Start"], row["End"])

            weight = row["NGS_read_count"] if weighted else 1
            nuc_overlap_dict[idx] += weight
            if overlap_seq in overlap_seq_dict:
                overlap_seq_dict[overlap_seq] += weight
            else:
                overlap_seq_dict[overlap_seq] = weight

        y = np.array([*nuc_overlap_dict.values()])

        def expected(x):
            return pow(0.25, x) * 0.75
        expected_values = np.fromfunction(expected, (11,), dtype=float)
        plt.figure()
        plt.bar(x=nuc_overlap_dict.keys(), height=y/np.sum(y), width=-0.4, align="edge")
        plt.bar(x=list(nuc_overlap_dict.keys()), height=expected_values, width=0.4, align="edge")

        plt.xlabel("number of overlapping nucleotides")
        plt.ylabel("relative occurrence")
        plt.title(f"nucleotide overlap at junction site for {k}")

        savepath = os.path.join(RESULTSPATH, "overlapping_nucleotides", f"{k}.pdf")
        plt.savefig(savepath)
        
    plot_dict = dict()
    for key, value in overlap_seq_dict.items():
        if value >= 10:
            plot_dict[key] = value

    plot_dict = {k: v for k, v in sorted(plot_dict.items(), key=lambda item: item[1], reverse=True)}

    plt.figure(figsize=(15, 5))
    plt.bar(x=plot_dict.keys(), height=plot_dict.values(), width=0.5)

    plt.xlabel("nucleotide sequence")
    plt.ylabel("absolute occurrence (without weight)")
    plt.title(f"occurrence of sequences that overlap at start and end")

    savepath = os.path.join(RESULTSPATH, "overlapping_nucleotides", f"sequence_distribution.pdf")
    plt.savefig(savepath)

