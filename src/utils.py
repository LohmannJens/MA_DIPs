'''
    General functions and global parameters, that are used in different scripts
'''
import os
import random
import json

import numpy as np
import pandas as pd

from Bio import SeqIO


# load config and assign values to global variables
DATAPATH = json.load(open("../../.config.json"))["DATAPATH"]
RESULTSPATH = json.load(open("../../.config.json"))["RESULTSPATH"]

# segments, nuclotides, and strains
SEGMENTS = list(["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"])
NUCLEOTIDES = list(["A", "C", "G", "U"])
STRAINS = dict({"Cal07": "A/California/07/2009",
                "NC": "A/New Caledonia/20-JY2/1999",
                "Perth": "A/Perth/16/2009",
                "BLEE": "B/Lee/1940",
                "PR8": "A/Puerto Rico/8/1934",
                "WSN": "A/WSN/1933"
                })

# global colors for plotting
COLORS = dict({"A": "deepskyblue", "C": "gold", "G": "springgreen", "U": "salmon"})

# parameters for the sampling
QUANT = 0.1
N_SAMPLES = 2000

def get_sequence(strain: str, seg: str, full: bool=False)-> object:
    '''
        Loads a DNA sequence given the strain and segment.
        :param strain: name of the strain
        :param seg: name of the segment
        :param full: if True the whole Biopython Seq Object is returned
                     if False a string object is returned

        :return: Biopython Seq Object or str() of the sequence
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
        Calculates the length of a specific sequence given the strain and
        segment.
        :param strain: name of the strain
        :param seg: name of the segment

        :return: length of the sequence as int
    '''
    return len(get_sequence(strain, seg))

def load_alnaji_excel()-> dict:
    '''
        Loads the excel file of Alnaji2019 containing the start and end
        positions of the deletion sites.
        Cleans up NaN data fields and splits up the strains and the lines.
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
    # dataset names are "Cal07", "NC", "Perth", "BLEE"
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
        Gets the loaded excel file of Alnaji 2019 and joins the lineages by
        strain. Dict with eight key-value pairs gets reduced to four pairs
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
        Loads the short reads of Alnaji 2019 from extra excel file and adds
        them to an existing dictionary.
        :param data_dict: dictionary with longer deletions

        :return: combined dictionary of both sources
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

        :return: letter indicating the significance level
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

def load_pelz_dataset(de_novo: bool=False,
                      long_dirna: bool=False,
                      by_time: bool=False)-> dict:
    '''
        Loads the data from Pelz et al 2021 publication.
        Is structured the same way as data from Alnaji 2019.
        :param de_novo: if True only de novo candidates are taken
        :param long_dirna: if True loads data set that includes long DI RNA
                           candidates
        :param by_time: if True loads the dataset split up by timepoints

        :return: dictionary with one key, value pair
    '''
    if long_dirna:
        filename = "NGS_SC_3_cutoffMPI_numbers_long_DI_RNAs.xlsx"
    elif by_time:
        filename = "ShortDeletions_by_timepoints.xlsx"
    else:
        filename = "ShortDeletions_AbsoluteValues.xlsx"
    file_path = os.path.join(DATAPATH, "Pelz2021", filename)
    data_dict = pd.read_excel(io=file_path,
                              sheet_name=None,
                              header=0,
                              na_values=["", "None"],
                              keep_default_na=False)

    if de_novo:
        d = data_dict["PR8"]
        d = d[d["class"].isin(["de_novo_loss", "de_novo_gain"])]
        data_dict["PR8"] = d

    return data_dict


def generate_sampling_data(seq: str, s: (int, int), e: (int, int),  n: int) -> object:
    '''
        Generates sampling data by creating random start and end points for
        artificial deletion sites. Generated data is used to calculate the
        expected values.
        :param seq: RNA sequence
        :param s: tuple with start and end point of the range for the artifical
                  start point of the deletion site
        :param e: tuple with start and end point of the range for the artifical
                  end point of the deletion site
        :param n: number of samples to generate

        :return: dataframe with the artifical data set
    '''
    sampling = dict({"Start": [], "End": []})
    for _ in range(n):
        sampling["Start"].append(random.randint(s[0], s[1]))
        sampling["End"].append(random.randint(e[0], e[1]))
    return pd.DataFrame(data=sampling)

def create_sequence_library(data_dict: dict)-> dict:
    '''
        Gets the raw loaded sequence data, which is a dict over all strains.
        In each dict the value is a data frame including DI RNA candidates.
        Creates the DI RNA sequence and adds it to the data frame.
        :param data_dict: dictionary key is strain names, value is df of DI RNA
                          candiates

        :return: dictionary with key for each strain. Value is a pandas df.
    '''
    for k, v in data_dict.items():
        del_seq_list = list()
        for i, row in v.iterrows():
            full_seq = get_sequence(k, row["Segment"])
            del_seq = full_seq[:row["Start"]] + full_seq[row["End"]-1:]
            del_seq_list.append(del_seq)

        data_dict[k]["DIRNASequence"] = del_seq_list

    return data_dict

def load_kupke(corrected: bool)-> dict:
    '''
        Loads the data set of Kupke et al. 2020. Does not load unused columns.
        Returns a dictionary with the data.
        :param corrected: indicates if the corrected dataset should be loaded

        :return: dictionary with strain name as key and data frame as value
    '''
    if corrected:
        path = os.path.join(DATAPATH, "Kupke2020", "corrected_data.xlsx")
        data = pd.read_excel(path)
    else:
        path = os.path.join(DATAPATH, "Kupke2020", "Table_S5_chim_junctions_h1n1.csv")
        data = pd.read_csv(path, usecols=["junctions", "infection"])

        replace_dict = dict({"H1N1_seg1": "PB2", "H1N1_seg2": "PB1", "H1N1_seg3": "PA"})

        data[["Segment", "DI"]] = data["junctions"].str.split(":", expand=True)
        data.replace(to_replace=replace_dict, inplace=True)
        data.drop(["junctions"], axis=1, inplace=True)
    
    dic = dict({"PR8": data})
    return dic

def load_alnaji_2021()-> dict:
    '''
        Loads the data set of Alnaji et al. 2021. Returns a dictionary with the
        data.

        :return: dictionary with strain name as key and data frame as value
    '''
    path = os.path.join(DATAPATH, "Alnaji2021", "Early_DIs_mbio.xlsx")
    data = pd.read_excel(path, na_values=["", "None"], keep_default_na=False)
    dic = dict({"PR8": data})
    return dic

def load_hpi14_alnaji()-> dict:
    '''
        Loads the hpi 14 data set of Alnaji et al. 2021. Returns a dictionary
        with the data.

        :return: dictionary with strain name as key and data frame as value
    '''
    path = os.path.join(DATAPATH, "Alnaji2021", "14hpi_inter_exter.xlsx")
    data = pd.read_excel(path, na_values=["", "None"], keep_default_na=False)
    dic = dict({"PR8": data})
    return dic

def load_full_alnaji2021()-> dict:
    '''
        Loads the data for 3, 6 and 24 hpi and the data for both 14 hpis and
        merges them together.

        :return: dictionary with strain name as key and data frame as value
    '''
    data_dict = load_alnaji_2021()
    hpi14_dict = load_hpi14_alnaji()
    data_df = pd.concat([data_dict["PR8"], hpi14_dict["PR8"]]).reset_index()
    data_df.drop(columns=["index"], inplace=True)

    data_df.loc[(data_df["Timepoint"] == "14hpi") & (data_df["Class"] == "internal"), "Timepoint"] = "14hpi_internal"
    data_df.loc[(data_df["Timepoint"] == "14hpi") & (data_df["Class"] == "external"), "Timepoint"] = "14hpi_external"

    return data_df

