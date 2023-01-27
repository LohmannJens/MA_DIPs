'''
    General functions and global parameters, that are used in different scripts
    of the ML part. Includes functions to load data and generate features.
'''
import os
import sys
import RNA

import numpy as np
import pandas as pd

from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split

sys.path.insert(0, "..")
from utils import DATAPATH, SEGMENTS
from utils import load_pelz_dataset, load_kupke, load_full_alnaji2021, load_short_reads, load_alnaji_excel, get_sequence, get_seq_len

sys.path.insert(0, "../direct_repeats")
from search_direct_repeats import calculate_direct_repeat

### global paramters ###
CHARS = 'ACGU'
CHARS_COUNT = len(CHARS)
MAX_LEN = 2361 # B Lee


### feature generation ###
def generate_features(df: pd.DataFrame,
                      features: list,
                      load_precalc: bool
                      )-> (pd.DataFrame, list):
    '''
        Main function to generate/load the features.
        :param df: data frame including the deletion site data
        :param features: list indicating which features to include
        :param load_precalc: if True features are loaded from precalculated
                             file

        :return: data frame including the features,
                 list with the names of the columns for the features
    '''
    if load_precalc:
        # load precalculated data from file
        path = os.path.join(DATAPATH, "ML", "precalc.csv")
        df = pd.read_csv(path, na_values=["", "None"], keep_default_na=False)

        # select correct columns from a file?
        feature_cols = ["Start", "End"]
        if "Segment" in features:
            segment_cols = [f"Segment_{s}" for s in SEGMENTS]
            feature_cols = feature_cols + segment_cols
        if "DI_Length" in features:
            feature_cols.append("DI_Length")
        if "Direct_repeat" in features:
            feature_cols.append("Direct_repeat")
        if "Junction" in features:
            start_cols = [f"Start_{i}_{ch}" for i in range(1, 11) for ch in CHARS]
            end_cols = [f"End_{i}_{ch}" for i in range(1, 11) for ch in CHARS]
            feature_cols = feature_cols + start_cols + end_cols
        if "3_5_ratio" in features:
            feature_cols.append("3_5_ratio")
        if "length_proportion" in features:
            feature_cols.append("length_proportion")
        if "full_sequence" in features:
            sequence_cols = [f"{i}_{ch}" for i in range(1, MAX_LEN+1) for ch in CHARS]
            feature_cols = feature_cols + sequence_cols
        if "delta_G" in features:
            feature_cols.append("delta_G")

    else:
        feature_cols = ["Start", "End"]
        if "Segment" in features:
            df, segment_cols = segment_ohe(df)
            feature_cols = feature_cols + segment_cols
        if "DI_Length" in features:
            df["DI_Length"] = df.apply(get_dirna_length, axis=1)
            feature_cols.append("DI_Length")
        if "Direct_repeat" in features:
            df["Direct_repeat"] = df.apply(get_direct_repeat_length, axis=1)
            feature_cols.append("Direct_repeat")
        if "Junction" in features:
            df, junction_start_cols = junction_site_ohe(df, "Start")
            df, junction_end_cols = junction_site_ohe(df, "End")
            feature_cols = feature_cols + junction_start_cols + junction_end_cols
        if "3_5_ratio" in features:
            df["3_5_ratio"] = df.apply(get_3_to_5_ratio, axis=1)
            feature_cols.append("3_5_ratio")
        if "length_proportion" in features:
            df["length_proportion"] = df.apply(get_length_proportion, axis=1)
            feature_cols.append("length_proportion")
        if "full_sequence" in features:
            df, sequence_cols = full_sequence_ohe(df)
            feature_cols = feature_cols + sequence_cols
        if "delta_G" in features:
            df["delta_G"] = df.apply(get_delta_G, axis=1)
            feature_cols.append("delta_G")

    return df, feature_cols

def segment_ohe(df: pd.DataFrame)-> (pd.DataFrame, list):
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

def get_dirna_length(row: pd.Series)-> int:
    '''
        Calculates the length of the DI RNA sequence given a row of a data
        frame with the necessary data.
        :param row: data frame row including Strain, Segment, Start, and End
        
        :return: length of DI RNA sequence
    '''
    seq_len = get_seq_len(row["Strain"], row["Segment"])
    return row["Start"] + (seq_len - row["End"] + 1)

def get_direct_repeat_length(row: pd.Series)-> int:
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

def junction_site_ohe(df: pd.DataFrame,
                      position: str
                      )-> (pd.DataFrame, list):
    '''
        Gets the sequence around the start or end of a given junction site and
        converts the sequence into an one hot encoding.
        :param df: data frame including Start, End, Strain, and Segment
        :param position: is either 'Start' or 'End' to indicate which site
        
        :return: Tuple with two entries:
                    data frame including original data and OHE data
                    list with the column names of the OHE
    '''
    # initializing matrix
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

def get_3_to_5_ratio(row: pd.Series)-> float:
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

def get_length_proportion(row: pd.Series)-> float:
    '''
        Calculates the proportion of the length of the DI RNA sequence to the
        full length sequence given a row of a data frame.
        :param row: data frame row including Strain, Segment, Start, and End
        
        :return: ratio of DI RNA lenght to full length sequence
    '''
    seq_len = get_seq_len(row["Strain"], row["Segment"])
    dirna_len = row["Start"] + (seq_len - row["End"] + 1)
    return dirna_len/seq_len

def full_sequence_ohe(df: pd.DataFrame)-> (pd.DataFrame, list):
    '''
        Gets the whole sequence as an one hot encoding. Sequences get
        normalized to the longest sequence length by adding * at the end
        :param df: data frame including Start, End, Strain, and Segment
        
        :return: Tuple with two entries:
                    data frame including original data and OHE data
                    list with the column names of the OHE
    '''
    # defining initializing matrix
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

def get_delta_G(row: pd.Series)-> float:
    '''
        :param row: data frame row including Strain, Segment, Start, and End
        
        :return: ratio of DI RNA lenght to full length sequence
    '''
    seq = get_sequence(row["Strain"], row["Segment"])
    del_seq = seq[:row["Start"]] + seq[row["End"]-1:]
    mfe = RNA.fold_compound(del_seq).mfe()[1]
    return mfe/len(del_seq)

### others ###

def get_duplicate_info(df: pd.DataFrame)-> pd.DataFrame:
    '''
        Adds a new column to a given data frame, that shows if the DI
        candidate is a duplicate (1) or not (0). 
        :param df: data frame

        :return: Data frame with a new column called 'Duplicate'
    '''
    df["DI"] = df["Segment"] + "_" + df["Start"].astype(str) + "_" + df["End"].astype(str)
    t_df = pd.DataFrame(df.groupby(["DI"]).size())
    t_df = t_df.rename(columns={0: "Occurrences"})
    is_dupl_list = t_df[t_df["Occurrences"] > 1].index.values.tolist()

    df["Duplicate"] = 0
    df.loc[df["DI"].isin(is_dupl_list), "Duplicate"] = 1
    df.drop(["DI"], axis=1, inplace=True)
    df["comb_dup"] = ((df["Duplicate"] == 1) | (df["int_dup"] == 1))

    return df

def load_all_sets()-> pd.DataFrame:
    '''
        Loads all data sets together in one data frame. Provides the columns
        Segment, Start, End, NGS_read_count and dataset_name.

        :return: Data frame including all available data sets
    '''
    def log_and_norm(df):
        df["NGS_read_count"] = df["NGS_read_count"].astype(float)
        df = df[df["NGS_read_count"] > 0].copy()
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
    # only load post samples to avoid duplicates
    kupke = kupke[kupke.Infection == "post"]

    kupke.drop(["DI", "Length", "Infection", "Num_sample", "Correction"], axis=1, inplace=True)
    kupke = merge_duplicates(kupke)
    kupke["dataset_name"] = "Kupke"
    kupke["Strain"] = "PR8"
    kupke = log_and_norm(kupke)
    df = pd.concat([df, kupke])

    df["int_dup"] = 0

    # load alnaji 2021 dataset
    alnaji2021 = load_full_alnaji2021()

    t_df = pd.DataFrame(alnaji2021.groupby(["DI"]).size())
    t_df = t_df.rename(columns={0: "Occurrences"})
    dupl_list = t_df[t_df["Occurrences"] > 1].index.values.tolist()

    alnaji2021 = merge_duplicates(alnaji2021)
    alnaji2021["DI"] = alnaji2021["Segment"] + "_" + alnaji2021["Start"].astype(str) + "_" + alnaji2021["End"].astype(str)
    alnaji2021["int_dup"] = 0
    alnaji2021.loc[alnaji2021["DI"].isin(dupl_list), "int_dup"] = 1
    alnaji2021.drop(["DI"], axis=1, inplace=True)

    alnaji2021["dataset_name"] = "Alnaji2021"
    alnaji2021["Strain"] = "PR8"
    alnaji2021 = log_and_norm(alnaji2021)
    df = pd.concat([df, alnaji2021])

    # load four datasets of alnaji 2019
    alnaji2019 = load_short_reads(load_alnaji_excel())
    for k, v in alnaji2019.items():
        v["DI"] = v["Segment"] + "_" + v["Start"].astype(str) + "_" + v["End"].astype(str)
        t_df = pd.DataFrame(v.groupby(["DI"]).size())
        t_df = t_df.rename(columns={0: "Occurrences"})
        dupl_list = t_df[t_df["Occurrences"] > 1].index.values.tolist()

        v["NGS_read_count"] = v["NGS_read_count"].astype(int)
        v = merge_duplicates(v)
        v["DI"] = v["Segment"] + "_" + v["Start"].astype(str) + "_" + v["End"].astype(str)
        v["int_dup"] = 0
        v.loc[v["DI"].isin(dupl_list), "int_dup"] = 1

        v.drop(["DI"], axis=1, inplace=True)
        v["dataset_name"] = f"Alnaji2019_{k}"
        v["Strain"] = k
        v = log_and_norm(v)
        df = pd.concat([df, v])

    df.reset_index(inplace=True)
    df.drop(["index"], axis=1, inplace=True)

    df = get_duplicate_info(df)

    return df

def duplicates_set_labels(df: pd.DataFrame,
                          col: str
                          )-> pd.Series:
    '''
        Sets the labels for the classifier if the duplicates should be
        predicted.
        :param df: data frame
        :param col: name of the column to use as y value

        :return: pandas Series containing the labels (0 and 1)
    '''
    return df[col]

def ngs_set_labels(df: pd.DataFrame,
                   n_bins: int,
                   style: str,
                   labels: list=[]
                   )-> pd.Series:
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

def select_datasets(df: pd.DataFrame,
                    t_datasets: list,
                    v_datasets: list,
                    features: list,
                    n_bins: int,
                    label_style: str,
                    y_column: str
                    )-> (pd.DataFrame, pd.Series, pd.DataFrame, pd.Series):
    '''
        Selects training a test data by a given name.
        :param df: pandas data frame including all data sets and features
        :param t_datasets: list of datasets to use for training
        :param v_datasets: list of datasets to use for valdation
        :param features: list with all features, that should be selected
        :param n_bins: number of classes to create
        :param label_style: declares how to create the labels/classes
        :param y_column: indicates the columne where to take the y values from

        :return: tuple with 4 entries, where each is a pandas data frame
                    X:     input data for training
                    y:     True labels for training
                    X_val: input data for validation
                    y_val: True labels for validation
    '''
    labels = ["low", "high"]
    if n_bins == 3:
        labels.insert(1, "mid")

    t_df = df.loc[df["dataset_name"].isin(t_datasets)].copy().reset_index()
    X = t_df[features]
    if y_column == "NGS_log_norm":
        y = ngs_set_labels(t_df, n_bins, label_style, labels)
    else:
        y = duplicates_set_labels(t_df, y_column)

    if len(v_datasets) != 0:
        v_df = df.loc[df["dataset_name"].isin(v_datasets)].copy().reset_index()
        X_val = v_df[features]
        if y_column == "NGS_log_norm":
            y_val = ngs_set_labels(v_df, n_bins, label_style, labels)
        else:
            y_val = duplicates_set_labels(v_df, y_column)
    else:
        X, X_val, y, y_val = train_test_split(X, y, test_size=0.2, random_state=42)

    return X, y, X_val, y_val

