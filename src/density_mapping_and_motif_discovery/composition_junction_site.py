'''
Loads the data for the deletion sides (start/end point) from Alnaji 2019.
Takes a look at the nucleotide distribution around these points.
Goal is to see if any of the four nucleotides occure more/less often.
'''
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from Bio.Seq import Seq

sys.path.insert(0, "..")
from utils import SEGMENTS, load_excel, load_short_reads, get_sequence


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


data_folder = os.path.join("..", "..", "data", "alnaji2019")
filepath = os.path.join(data_folder, "DI_Influenza_FA_JVI.xlsx")
cleaned_data_dict = load_excel(filepath)

short_reads_filepath = os.path.join(data_folder, "Small_deletionSize_FA.xlsx")
all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)

# Create a sequence library for each strain
sequence_list_dict = dict()
nuc_count = dict({"A": 0, "C": 0, "G": 0, "U": 0})

for key, value in all_reads_dict.items():
    sequence_list = list()
    # loop over rows of value dataframe
    for i, row in value.iterrows():
        length = int(row["Length"])
        count = int(row["NGS_read_count"])
        start = int(row["Start"])
        end = int(row["End"])
        segment = row["Segment"]
        del_sequence = create_sequence(start, end, key, segment, crop=True)
        whole_sequence = create_sequence(start, end, key, segment, crop=False)
    
        entry_dict = dict({"Length": length,
                           "Count": count,
                           "Start": start,
                           "End": end,
                           "Segment": segment,
                           "DelSequence": del_sequence,
                           "WholeSequence": whole_sequence,
                           "NucleotideCount": nuc_count})
        sequence_list.append(entry_dict)
    
        for n in NUCLEOTIDES:
            nuc_count[n] = nuc_count[n] + whole_sequence.count(n)
    sequence_list_dict[key] = sequence_list

'''
# Loop over the different strains
for k, v in sequence_list_dict.items():

    count_start_dict = dict({"A": np.zeros(9), "C": np.zeros(9), "G": np.zeros(9), "U": np.zeros(9)})
    count_end_dict = dict({"A": np.zeros(9), "C": np.zeros(9), "G": np.zeros(9), "U": np.zeros(9)})

    normalize = 0

    for entry in v:
        seq_start_dict = count_nucleotide_occurrence(entry["WholeSequence"], entry["Start"]) 
        seq_end_dict = count_nucleotide_occurrence(entry["WholeSequence"], entry["End"])
        normalize += entry["Count"]
        for nuc in count_start_dict.keys():
            count_start_dict[nuc] += seq_start_dict[nuc] * entry["Count"]
            count_end_dict[nuc] += seq_end_dict[nuc] * entry["Count"]

    fig, axs = plt.subplots(2, 1, figsize=(5, 10), tight_layout=True)
    x = np.arange(0.7, 9.7, dtype=np.float64)

    for key in count_start_dict.keys():
        axs[0].bar(x, height=count_start_dict[key]/normalize, width=0.2, label=key)
        axs[1].bar(x, height=count_end_dict[key]/normalize, width=0.2, label=key)
        x += 0.2

    colors = dict({"A": "blue", "C": "yellow", "G": "green", "U": "red"})

    for i in range(2):
        for n in NUCLEOTIDES:
            axs[i].axhline(y=float(nuc_count[n]/sum(nuc_count.values())), c=colors[n], ls="--", lw=0.5)
        axs[i].legend()
        axs[i].set_ylim(top=0.8)
        axs[i].set_xticks([1,2,3,4,5,6,7,8,9])
        axs[i].set_xlabel("position at junction side")
        axs[i].set_ylabel("relative occurrence")
        pos = "start" if i == 0 else "end"
        axs[i].set_title(f"relative occurrence of nucleotides at {pos} of deletion site")

    savepath = os.path.join("results", f"{k}_relative_occurrence_nucleotides.pdf")
    plt.savefig(savepath)
'''

def calculate_overlapping_nucleotides(seq: str, s: int, e: int)-> int:
    '''
    calculates the number of overlapping nucleotides directly before start
        and end of junction site.
        :param seq: sequence where the overlap will be calculated on
        :param s: start point of the junction site
        :param e: end point of the junction site

        :return: integer indicating the length of the overlapping sequence
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
    return counter

####
#TODO
#   Check if nucleotides directly before junction site have the same sequence
#   as the ones directly before junction site at the end
for k, v in sequence_list_dict.items():
    nuc_overlap_dict = dict({0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10:0})

    for seq_dict in v:
        sequence = seq_dict["WholeSequence"]
        i = calculate_overlapping_nucleotides(sequence, seq_dict["Start"], seq_dict["End"])
        nuc_overlap_dict[i] += seq_dict["Count"]

    y = np.array([*nuc_overlap_dict.values()])

    def expected(x):
        return pow(0.25, x) * 0.75

    expected_values = np.fromfunction(expected, (11,), dtype=float)
    print(expected_values)
    print(nuc_overlap_dict.keys())
    plt.figure()
    plt.bar(x=nuc_overlap_dict.keys(), height=y/np.sum(y), width=-0.4, align="edge")
    plt.bar(x=list(nuc_overlap_dict.keys()), height=expected_values, width=0.4, align="edge")

    plt.xlabel("number of overlapping nucleotides")
    plt.ylabel("relative occurrence")
    plt.title(f"nucleotide overlap at junction site for {k}")

    savepath = os.path.join("results", f"{k}_overlapping_nucleotides.pdf")
    plt.savefig(savepath)



