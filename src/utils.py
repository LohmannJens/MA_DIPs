import os

import numpy as np
import pandas as pd

from Bio import SeqIO

# needs to be changed when cloning the repo
REPOPATH = "/home/jens/Masterarbeit/MA_DIPs"

DATAPATH = os.path.join(REPOPATH, "data")
RESULTSPATH = os.path.join(REPOPATH, "results")

SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]


def get_sequence(strain: str, seg: str)-> object:
    '''
        loads a DNA sequence by strain and segment.
        :param strain: name of the strain
        :param seg: name of the segment of the strain
    '''
    fasta_file = os.path.join(DATAPATH, "alnaji2019", strain, f"{seg}.fasta")
    return SeqIO.read(fasta_file, "fasta")

def get_seq_len(strain: str, seg: str)-> int:
    '''
        calculates the length of a specific segment.
        uses a fasta file for that.
        :param strain: name of the strain
        :param seg: name of the segment of the strain

        :return: length of the sequence as int
    '''
    return len(get_sequence(strain, seg))

def load_alnaji_excel()-> dict:
    '''
        loads the excel file containing the start and end positions of 
        the deletion sides.
        Cleans up nan data fields and splits up the strains and the lines.
        Also calculates the length of each deletion and adds it as new column.

        :return: dictionary with 8 key, value pairs (4 strains * 2 lines each)
    '''
    file_path = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    data_dict = pd.read_excel(io=file_path,
                              sheet_name=None,
                              header=1,
                              na_values=["", "None"],
                              keep_default_na=False,
                              converters={"Start": int,"End": int, "NGS_read_count": int,
                                          "Start.1": int,"End.1": int, "NGS_read_count.1": int})
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
            seq_length = get_seq_len(key[:-3], s)
            length = value["Start"] + (seq_length - value["End"] + 1)
            value.loc[value["Segment"] == s, "Length"] = length

    return cleaned_data_dict

def join_lineages(data: dict)-> dict:
    '''
        gets the loaded excel file of Alnaji 2019 and joins the lineages by
        strains. Dict with eight key-value pairs gets reduced to four pairs
        :param data: loaded dict with strain as key and data frame as value

        :return: dict with dataframes joined by strain name
    '''
    merged_dict = dict()
    for k, v in data.items():
        if k in merged_dict:
            merged_dict[k[:-3]] = pd.concat([merged_dict[k[:-3]], v])
        else:
            merged_dict[k[:-3]] = v
    return merged_dict

def load_short_reads(data_dict: dict)-> dict:
    '''
        loads the short reads from extra excel file and adds them to an
        existing dictionary.
        :param data_dict: dictionary with longer deletions

        :return: gives a combined dictionary of both sources
    '''
    file_path = os.path.join(DATAPATH, "alnaji2019", "Small_deletionSize_FA.xlsx")
    short_data_dict = pd.read_excel(io=file_path,
                                    sheet_name=None,
                                    header=0,
                                    na_values=["", "None"],
                                    keep_default_na=False,
                                    converters={"Start": int,"End": int, "NGS_read_count": int})

    # calculate length for all seqments and add as extra column
    for k, v in short_data_dict.items():
        v["Length"] = np.nan
        for s in SEGMENTS:
            seq_length = get_seq_len(k, s)
            length = v["Start"] + (seq_length - v["End"] + 1)
            v.loc[v["Segment"] == s, "Length"] = length
    
    # merge each dataframe to the loaded dictionary with the short deletions
    for k, v in data_dict.items():
        short_data_dict[k[:-3]] = pd.concat([short_data_dict[k[:-3]], v])

    return short_data_dict

def get_stat_symbol(p: float)-> str:
    '''
        Indicates the statistical significance by letters. Is used for plots.
        :param p: p-value of the test

        :return: letter indicating the significance
    '''
    if p == 0.0:
        return ""
    elif p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

def load_pelz_dataset()-> object:
    '''
        Loads the data from Pelz et al 2019 publication.
        Is structured the same way as data from alnaji 2019.

        :return: dictionary with one key, value pair
    '''
    file_path = os.path.join(DATAPATH, "Pelz2021", "ShortDeletions_AbsoluteValues.xlsx")
    data_dict = pd.read_excel(io=file_path,
                              sheet_name=None,
                              header=0,
                              na_values=["", "None"])

    return data_dict

