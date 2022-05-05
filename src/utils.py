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

        :return: length of the sequence as int
    '''
    fasta_file = os.path.join("..", "..", "data", "alnaji2019", strain[:-3], f"{seg}.fasta")
    return len(SeqIO.read(fasta_file, "fasta"))

def load_excel(path: str)-> dict:
    '''
        loads the excel file containing the start and end positions of 
        the deletion sides.
        Cleans up nan data fields and splits up the strains and the lines.
        Also calculates the length of each deletion and adds it as new column.
        :param path: path to the file location

        :return: dictionary with 8 key, value pairs (4 strains * 2 lines each)
    '''
    data_dict = pd.read_excel(io=path,
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
            length = value["Start"] + (seq_length - value["End"] + 1)
            value.loc[value["Segment"] == s, "Length"] = length

    return cleaned_data_dict

def load_short_reads(data_dict: dict, path: str)-> dict:
    '''

    '''

    # combine lineage 1 and 2 for each of the 4 strains


    # load short read excel as dict of dataframes


    # merge every dataframe in each of the four key, value pairs


