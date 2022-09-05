'''
    General functions and global parameters, that are used in different scripts
    Functions include: loading of sequences, loading of datasets, ...
    Parameters include: paths to results and data, Segments names, ...
'''
import os
import random

import numpy as np
import pandas as pd

from Bio import SeqIO

# needs to be changed when cloning the repo
REPOPATH = "/home/jens/Masterarbeit/MA_DIPs"

# paths to the data and results folder
DATAPATH = os.path.join(REPOPATH, "data")
RESULTSPATH = os.path.join(REPOPATH, "results")

SEGMENTS = list(["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"])
COLORS = dict({"A": "deepskyblue", "C": "gold", "G": "springgreen", "U": "salmon"})
NUCLEOTIDES = list(["A", "C", "G", "U"])

# parameters for the sampling
QUANT = 0.1
S_ROUNDS = 5

def get_sequence(strain: str, seg: str, full: bool=False)-> object:
    '''
        loads a DNA sequence by strain and segment.
        :param strain: name of the strain
        :param seg: name of the segment of the strain
        :param full: if true the whole Biopython Seq Object is returned

        :return: Biopython Seq Object including the sequence. To get raw string
                 use: str(RETURN.seq)
    '''
    fasta_file = os.path.join(DATAPATH, "strain_segment_fastas", strain, f"{seg}.fasta")
    seq_obj = SeqIO.read(fasta_file, "fasta")
    if full:
        return seq_obj
    else:
        return str(seq_obj.seq.transcribe())

    return

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
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

def load_pelz_dataset()-> dict:
    '''
        Loads the data from Pelz et al 2019 publication.
        Is structured the same way as data from alnaji 2019.

        :return: dictionary with one key, value pair
    '''
    file_path = os.path.join(DATAPATH, "Pelz2021", "ShortDeletions_AbsoluteValues.xlsx")
    data_dict = pd.read_excel(io=file_path,
                              sheet_name=None,
                              header=0,
                              na_values=["", "None"],
                              keep_default_na=False)

    return data_dict

def generate_sampling_data(seq: str, s: (int, int), e: (int, int),  n: int) -> object:
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
    sampling = dict({"Start": [], "End": []})
    for _ in range(n):
        sampling["Start"].append(random.randint(s[0], s[1]))
        sampling["End"].append(random.randint(e[0], e[1]))
    return pd.DataFrame(data=sampling)

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
        del_seq_list = list()
        for i, row in v.iterrows():
            full_seq = get_sequence(k, row["Segment"])
            del_seq = full_seq[:row["Start"]] + full_seq[row["End"]-1:]
            del_seq_list.append(del_seq)

        data_dict[k]["DelSequence"] = del_seq_list

    return data_dict

def load_kupke()-> dict:
    '''
        Loads the data set of Kupke et al. 2020. Does not load unused columns.
        Returns a dictionary with the data.

        :return: dictionary with strain name as key and data frame as value
    '''
    path = os.path.join(DATAPATH, "Kupke2020", "Table_S5_chim_junctions_h1n1.csv")
    data = pd.read_csv(path, usecols=["junctions", "infection", "num_reads"])

    replace_dict = dict({"H1N1_seg1": "PB2", "H1N1_seg2": "PB1", "H1N1_seg3": "PA"})

    data[["Segment", "DI"]] = data["junctions"].str.split(":", expand=True)
    data.replace(to_replace=replace_dict, inplace=True)
    data.drop(["junctions"], axis=1, inplace=True)

    dic = dict({"PR8": data})
    return dic

