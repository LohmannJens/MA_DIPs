'''
    General functions and global parameters, that are used in different scripts
    Functions include: loading of sequences, loading of datasets, ...
    Parameters include: paths to results and data, Segments names, ...
'''
import os
import sys

import numpy as np
import pandas as pd

from sklearn.preprocessing import OneHotEncoder

sys.path.insert(0, "..")
from utils import load_pelz_dataset, load_kupke, load_full_alnaji2021, load_short_reads, load_alnaji_excel, get_sequence, get_seq_len

sys.path.insert(0, "../direct_repeats")
from search_direct_repeats import calculate_direct_repeat

### feature generation ###

def segment_ohe(df: object)-> (object, list):
    '''
        Converts the column with segment names into an one hot encoding.
        :param df: data frame including a row called 'Segment'
        :return: Tuple with two entries:
                    data frame including original data and OHE data
                    list with the column names of the OHE
    '''
    ohe = OneHotEncoder()
    segment_df = pd.DataFrame(ohe.fit_transform(df[["Segment"]]).toarray())
    ohe_cols = ohe.get_feature_names_out().tolist()
    segment_df.columns = ohe_cols
    df = df.join(segment_df)
    return df, ohe_cols

def junction_site_ohe(df: object, position: str)-> (object, list):
    '''
        Gets the sequence around the start or end of a given junction site and
        converts the sequence into an one hot encoding.
        :param df: data frame including Start, End, Strain, and Segment
        :param position: is either 'Start' or 'End' to indicate which site
        :return: Tuple with two entries:
                    data frame including original data and OHE data
                    list with the column names of the OHE
    '''
    # defining static parameters
    CHARS = 'ACGU'
    CHARS_COUNT = len(CHARS)
    n = df.shape[0]
    res = np.zeros((n, CHARS_COUNT * 10), dtype=np.uint8)

    # getting sequence window for each row and do OHE
    for i, r in df.iterrows():
        s = r[position]
        seq = get_sequence(r["Strain"], r["Segment"])
        seq = seq[s-5:s+5]
        # Write down OHE in numpy array
        for j, char in enumerate(seq):
            pos = CHARS.rfind(char)
            res[i][j*CHARS_COUNT+pos] = 1

    # format as data frame and create columns names of OHE
    encoded_df = pd.DataFrame(res)
    col_names = [f"{position}_{i}_{ch}" for i in range(1, 11) for ch in CHARS]
    encoded_df.columns = col_names
    df = df.join(encoded_df)

    return df, col_names

def full_sequence_ohe(df: object)-> (object, list):
    '''
        Gets the whole sequence as an one hot encoding. Sequences get
        normalized to the longest sequence length by adding * at the end
        :param df: data frame including Start, End, Strain, and Segment
        :return: Tuple with two entries:
                    data frame including original data and OHE data
                    list with the column names of the OHE
    '''
    # defining static parameters
    CHARS = 'ACGU'
    CHARS_COUNT = len(CHARS)
    MAX_LEN = 2361 # B Lee
    n = df.shape[0]
    res = np.zeros((n, CHARS_COUNT * MAX_LEN), dtype=np.uint8)

    # getting sequence window for each row and do OHE
    for i, r in df.iterrows():
        seq = get_sequence(r["Strain"], r["Segment"])
        seq = seq + "*" * (MAX_LEN - len(seq))
        # Write down OHE in numpy array
        for j, char in enumerate(seq):
            pos = CHARS.rfind(char)
            res[i][j*CHARS_COUNT+pos] = 1

    # format as data frame and create columns names of OHE
    encoded_df = pd.DataFrame(res)
    col_names = [f"{i}_{ch}" for i in range(1, MAX_LEN+1) for ch in CHARS]
    encoded_df.columns = col_names
    df = df.join(encoded_df)

    return df, col_names

def get_dirna_length(row: list)-> int:
    '''
        Calculates the length of the DI RNA sequence given a row of a data
        frame with the necessary data.
        :param row: data frame row including Strain, Segment, Start, and End
        :return: length of DI RNA sequence
    '''
    seq_len = get_seq_len(row["Strain"], row["Segment"])
    return row["Start"] + (seq_len - row["End"] + 1)

def get_direct_repeat_length(row)-> int:
    '''
        Calculates the length of the direct repeat given a row of a data frame
        with the necessary data.
        :param row: data frame row including Strain, Segment, Start, and End
        :return: length of direct repeat
    '''
    seq = get_sequence(row["Strain"], row["Segment"])
    s = row["Start"]
    e = row["End"]
    n, _ = calculate_direct_repeat(seq, s, e, 15, 1)
    return n

def get_3_to_5_ratio(row)-> float:
    '''
        Calculates the proportion of the 3' sequence to the 5' sequence given
        a row of a data frame.
        :param row: data frame row including Strain, Segment, Start, and End
        :return: ratio of 3' to 5' sequence length
    '''
    seq_len = get_seq_len(row["Strain"], row["Segment"])
    len3 = row["Start"]
    len5 = seq_len - row["End"] + 1
    return len3/len5

def get_length_proportion(row)-> float:
    '''
        Calculates the proportion of the length of the DI RNA sequence to the
        full length sequence given a row of a data frame.
        :param row: data frame row including Strain, Segment, Start, and End
        :return: ratio of DI RNA lenght to full length sequence
    '''
    seq_len = get_seq_len(row["Strain"], row["Segment"])
    dirna_len = row["Start"] + (seq_len - row["End"] + 1)
    return dirna_len/seq_len

### others ###

def load_all_sets()-> object:
    '''
        Loads all data sets together in one data frame. Provides the columns
        Segment, Start, End, NGS_read_count and dataset_name.

        :return: pandas data frame including all available data sets
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
    df["Strain"] = "PR8"
    df = log_and_norm(df)

    # load kupke dataset
    kupke = load_kupke(corrected=True)["PR8"]
    kupke.drop(["DI", "Length", "Infection", "Num_sample", "Correction"], axis=1, inplace=True)
    kupke = merge_duplicates(kupke)
    kupke["dataset_name"] = "Kupke"
    kupke["Strain"] = "PR8"
    kupke = log_and_norm(kupke)
    df = pd.concat([df, kupke])

    # load alnaji 2021 dataset
    alnaji2021 = load_full_alnaji2021()
    alnaji2021.drop(["DI", "Replicate", "Timepoint", "Class"], axis=1, inplace=True)
    alnaji2021 = merge_duplicates(alnaji2021)
    alnaji2021["dataset_name"] = "Alnaji2021"
    alnaji2021["Strain"] = "PR8"
    alnaji2021 = log_and_norm(alnaji2021)
    df = pd.concat([df, alnaji2021])

    # load four datasets of alnaji 2019
    alnaji2019 = load_short_reads(load_alnaji_excel())
    for k, v in alnaji2019.items():
        v.drop(["Length"], axis=1, inplace=True)
        v["NGS_read_count"] = v["NGS_read_count"].astype(int)
        v = merge_duplicates(v)
        v["dataset_name"] = f"Alnaji2019_{k}"
        v["Strain"] = k
        v = log_and_norm(v)
        df = pd.concat([df, v])

    df.reset_index(inplace=True)
    df.drop(["index"], axis=1, inplace=True)

    return df

def set_labels(df: object, n_bins: int, style: str, labels: list=[])-> object:
    '''
        Sets the labels for the classifer. Can be done by using pd.cut() or by
        using the median/33-percentil as split.
        :param df: data frame including the data
        :param n_bins: number of bins to use
        :param style: declares how to create the labels
        :param labels: list with labels to use, using 'pd.cut()'

        :return: pandas Series including the labels
    '''
    if style == "pd.cut":
        y = pd.cut(df["NGS_log_norm"], bins=n_bins, labels=labels, ordered=False)
    elif style == "median":
        y = list()
        if n_bins == 2:
            median = df["NGS_log_norm"].median()
            for row in df.iterrows():
                r = row[1]
                y.append("low" if r["NGS_log_norm"] < median else "high")
        elif n_bins == 3:
            perc1 = df["NGS_log_norm"].quantile(q=0.33)
            perc2 = df["NGS_log_norm"].quantile(q=0.66)
            for row in df.iterrows():
                r = row[1]
                if r["NGS_log_norm"] < perc1:
                    y.append("low")
                elif r["NGS_log_norm"] > perc2:
                    y.append("high")
                else:
                    y.append("mid")
        y = pd.Series(y)

    return y

def select_datasets(df, dataset_name: str, features: list, n_bins: int, label_style: str)-> (object, object, object, object):
    '''
        Selects training a test data by a given name.
        :param df: pandas data frame including all data sets and features
        :param dataset_name: string indicating which data sets to include
        :param features: list with all features, that should be selected
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes

        :return: tuple with 4 entries, where each is a pandas data frame
                    X:     input data for training
                    y:     True labels for training
                    X_val: input data for validation
                    y_val: True labels for validation
    '''

    if dataset_name == "Alnaji2019":
        train = ["Alnaji2019_Cal07", "Alnaji2019_NC", "Alnaji2019_Perth"]
        val = ["Alnaji2019_BLEE"]
    elif dataset_name == "PR8":
        train = ["Pelz", "Alnaji2021"]
        val = ["Kupke"]
    else:
        train = ["Alnaji2019_Cal07", "Alnaji2019_NC", "Alnaji2019_Perth",
                 "Alnaji2019_BLEE", "Pelz", "Alnaji2021", "Kupke"]
        val = list()

    labels = ["low", "high"]
    if n_bins == 3:
        labels.insert(1, "mid")

    t_df = df.loc[df["dataset_name"].isin(train)].copy().reset_index()
    X = t_df[features]
    y = set_labels(t_df, n_bins, label_style, labels)

    if len(val) != 0:
        v_df = df.loc[df["dataset_name"].isin(val)].copy().reset_index()
        X_val = v_df[features]
        y_val = set_labels(v_df, n_bins, label_style, labels)
    else:
        X_val = pd.DataFrame()
        y_val = pd.DataFrame()

    return X, y, X_val, y_val

