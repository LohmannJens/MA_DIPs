'''
Loads the Start and End points of the deletion sides from Alnaji 2019 and gives
insights about the data distribution.

1. Creates a histogram for each line in each strain containing the length of
   the deletion sides multiplied by their occurence.

2. Creates a plot where the start and end point of the deletion sides are
   plotted onto the sequence.
'''

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from sklearn.linear_model import LinearRegression

sys.path.insert(0, "..")
from utils import SEGMENTS, get_seq_len, load_excel, load_short_reads


data_folder = os.path.join("..", "..", "data", "alnaji2019")
filepath = os.path.join(data_folder, "DI_Influenza_FA_JVI.xlsx")
cleaned_data_dict = load_excel(filepath)

short_reads_filepath = os.path.join(data_folder, "Small_deletionSize_FA.xlsx")
all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)


# is used to check if the length of the deletions is very often modulo 3
# would be an indication that keeping the codons is important
seq_len_modulo_dict = dict({'0': 0, '1': 0, '2': 0})

# create a histogram for each line, indicating the length of the deletions
# just to get an overview about the data
for key, value in all_reads_dict.items():
    # create a dict for each segment including the NGS read count
    count_dict = dict()
    for s in SEGMENTS:
        count_dict[s] = dict()
    for i, r in value.iterrows():
        if r["Length"] in count_dict[r["Segment"]]:
            count_dict[r["Segment"]][r["Length"]] += r["NGS_read_count"]
        else:
            count_dict[r["Segment"]][r["Length"]] = r["NGS_read_count"]

    for seg_dict in count_dict.values():
        for k, v in seg_dict.items():
            seq_len_modulo_dict[str(int(k % 3))] += v

    # create a subplot for each key, value pair in count_dict
    fig, axs = plt.subplots(8, 1, figsize=(10, 20), tight_layout=True)
    fig.suptitle(f"absolute occurrences of deletions for {key}", x=0.3)
    for i, s in enumerate(SEGMENTS):
        axs[i].hist(count_dict[s].keys(), weights=count_dict[s].values(), bins=50)
        axs[i].set_title(f"{s}")
        axs[i].set_xlim(left=0)
        axs[i].set_xlabel("deletion length")
    save_path = os.path.join("results", f"{key}_length_del_hist.pdf")
    plt.savefig(save_path)

print(seq_len_modulo_dict)

# uses 'Start' and 'End' to indicate where the deletions happen on the sequence
for key, value in all_reads_dict.items():
    # create a dict for each segment using Start and End
    count_dict = dict()
    for s in SEGMENTS:
        count_dict[s] = dict()
    for i, r in value.iterrows():
        count_dict[r["Segment"]][r["Start"]] = r["NGS_read_count"]
        count_dict[r["Segment"]][r["End"]] = r["NGS_read_count"]

    # create a subplot for each key, value pair in count_dict
    fig, axs = plt.subplots(8, 1, figsize=(10, 20), tight_layout=True)
    fig.suptitle(f"position of deletions on sequence for {key}", x=0.3)
    for i, s in enumerate(SEGMENTS):
        axs[i].hist(count_dict[s].keys(), weights=count_dict[s].values(), bins=50)
        axs[i].set_title(f"{s}")
        axs[i].set_xlim(left=0)
        axs[i].set_xlabel("sequence position")
    save_path = os.path.join("results", f"{key}_del_position.pdf")
    plt.savefig(save_path)
