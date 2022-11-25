'''
    General functions and global parameters, that are used in different scripts
    Functions include: loading of sequences, loading of datasets, ...
    Parameters include: paths to results and data, Segments names, ...
'''
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "..")
from utils import load_pelz_dataset, load_kupke, load_full_alnaji2021, load_short_reads, load_alnaji_excel

def load_all_sets()-> object:
    '''
        Loads all data sets together in one data frame. Provides the columns
        Segment, Start, End, NGS_read_count and dataset_name.

        :return: data frame
    '''
    def log_and_norm(df):
        df["NGS_read_count"] = df["NGS_read_count"].astype(float)
        df = df[df["NGS_read_count"] > 0]
        df["NGS_log"] = np.log(df["NGS_read_count"]).astype(float)
        df["NGS_norm"] = df["NGS_read_count"]/max(df["NGS_read_count"])
        df["NGS_log_norm"] = df["NGS_log"]/max(df["NGS_log"])

        return df

    def merge_duplicates(df):
        df = df.groupby(["Segment", "Start", "End"]).sum(["NGS_read_count"]).reset_index()
        return df

    # load pelz dataset
    df = load_pelz_dataset()["PR8"]
    df["dataset_name"] = "Pelz"
    df = log_and_norm(df)

    # load kupke dataset
    kupke = load_kupke(corrected=True)["PR8"]
    kupke.drop(["DI", "Length", "Infection", "Num_sample", "Correction"], axis=1, inplace=True)
    kupke = merge_duplicates(kupke)
    kupke["dataset_name"] = "Kupke"
    kupke = log_and_norm(kupke)
    df = pd.concat([df, kupke])

    # load alnaji 2021 dataset
    alnaji2021 = load_full_alnaji2021()
    alnaji2021.drop(["DI", "Replicate", "Timepoint", "Class"], axis=1, inplace=True)
    alnaji2021 = merge_duplicates(alnaji2021)
    alnaji2021["dataset_name"] = "Alnaji2021"
    alnaji2021 = log_and_norm(alnaji2021)
    df = pd.concat([df, alnaji2021])

    # load four datasets of alnaji 2019
    alnaji2019 = load_short_reads(load_alnaji_excel())
    for k, v in alnaji2019.items():
        v.drop(["Length"], axis=1, inplace=True)
        v["NGS_read_count"] = v["NGS_read_count"].astype(int)
        v = merge_duplicates(v)
        v["dataset_name"] = f"Alnaji2019_{k}"
        v = log_and_norm(v)
        df = pd.concat([df, v])

    df.reset_index(inplace=True)
    df.drop(["index"], axis=1, inplace=True)

    return df

