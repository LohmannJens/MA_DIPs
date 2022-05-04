import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO


SEGMENTS = ["PB2", "PB1", "PA", "NP", "HA", "NA", "M", "NS"]


def get_seq_len(strain: str, seg:str)-> int:
    '''
        calculates the length of a specific segment.
        uses a fasta file for that.
        :param strain: name of the strain (including suffix "_l1" or "_l2")
        :param seg: name of the segment of the strain

        :return: returns the length of the sequence as int
    '''
    fasta_file = os.path.join("data", strain[:-3], f"{seg}.fasta")
    return len(SeqIO.read(fasta_file, "fasta"))


filepath = os.path.join("data", "DI_Influenza_FA_JVI.xlsx")
data_dict = pd.read_excel(io=filepath,
                          sheet_name=None,
                          header=1,
                          na_values=["", "None"],
                          keep_default_na=False)

# Splitting up the two lines in new data frames and cleaning NaN data
# For l2 the columns get renamed, they get the same names as in l1
# Cleaned data is stored in a dict, can be accessed by [datasetname]_[l1/l2]
# dataset names are "Cal07", "NC", "Perth", "B_LEE"
cleaned_data_dict = dict()
for key, value in data_dict.items():
    cleaned_data_dict[f"{key}_l1"] = data_dict[key].iloc[:, 0:4]
    cleaned_data_dict[f"{key}_l2"] = data_dict[key].iloc[:, 5:9]

    cleaned_data_dict[f"{key}_l1"].dropna(how="all", inplace=True)
    cleaned_data_dict[f"{key}_l2"].dropna(how="all", inplace=True)

    cleaned_data_dict[f"{key}_l2"].columns = cleaned_data_dict[f"{key}_l1"].columns 


# calculate length for all seqments and add as extra column
for key, value in cleaned_data_dict.items():
    value["Length"] = np.nan
    for s in SEGMENTS:
        seq_length = get_seq_len(key, s)
        value.loc[value["Segment"] == s, "Length"] = value["Start"] + (seq_length - value["End"] + 1)

# is used to check if the length of the deletions is very often modulo 3
# would be an indication that keeping the codons is important
seq_len_modulo_dict = dict({'0': 0, '1': 0, '2': 0})

# create a histogram for each line, indicating the length of the deletions
# just to get an overview about the data
for key, value in cleaned_data_dict.items():
    # create a dict for each segment using Start and End
    count_dict = dict()
    for s in SEGMENTS:
        count_dict[s] = dict()
    for i, r in value.iterrows():
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
        axs[i].set_xlim(0, get_seq_len(key, s))
        axs[i].set_xlabel("deletion length")
    save_path = os.path.join("results", f"{key}_length_del_hist.pdf")
    plt.savefig(save_path)

print(seq_len_modulo_dict)

# uses 'Start' and 'End' to indicate where the deletions happen on the sequence
for key, value in cleaned_data_dict.items():
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
        axs[i].set_xlim(0, get_seq_len(key, s))
        axs[i].set_xlabel("sequence position")
    save_path = os.path.join("results", f"{key}_del_position.pdf")
    plt.savefig(save_path)




