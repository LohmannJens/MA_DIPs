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


def create_sequence(s: int, e: int, strain: str, seg: str, crop: bool = False)-> str:
    '''
        Loads a sequence and removes a part of it. The deletion is given by
        a start and an end. Those positions are part of the remaining sequence.
        :param s: start of the deletion (included in the remaining sequence)
        :param e: end of the deletion (included in the remaining sequence)
        :param strain: name of the virus strain
        :param seg: segment where the sequenz is coming from
        :param crop: Indicates if the sequence should be cropped (True) to
                     exclude the deletion side or take whole sequence (False)

        :return: RNA sequence without deletion side
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
                           "WholeSequence": whole_sequence})
        sequence_list.append(entry_dict)
    sequence_list_dict[key] = sequence_list


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
# TODO: rethink if weighting by NGS count makes sense
            count_start_dict[nuc] += seq_start_dict[nuc] * entry["Count"]
            count_end_dict[nuc] += seq_end_dict[nuc] * entry["Count"]
    
    fig, axs = plt.subplots(2, 1, figsize=(5, 10), tight_layout=True)
    x = np.arange(0.7, 9.7, dtype=np.float64)

    for key in count_start_dict.keys():
        axs[0].bar(x, height=count_start_dict[key]/normalize, width=0.2, label=key)
        axs[1].bar(x, height=count_end_dict[key]/normalize, width=0.2, label=key)
        x += 0.2

    for i in range(2):
        axs[i].set_xticks([1,2,3,4,5,6,7,8,9])
        axs[i].legend()
        axs[i].set_xlabel("position at junction side")
        axs[i].set_ylabel("relative occurrence")
        pos = "start" if i == 0 else "end"
        axs[i].set_title(f"relative occurrence of nucleotides at {pos} of deletion site")
    
    savepath = os.path.join("results", f"{k}_relative_occurrence_nucleotides.pdf")
    plt.savefig(savepath)







