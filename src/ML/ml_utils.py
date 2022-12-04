'''
    General functions and global parameters, that are used in different scripts
    Functions include: loading of sequences, loading of datasets, ...
    Parameters include: paths to results and data, Segments names, ...
'''
import os
import sys

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

sys.path.insert(0, "..")
from utils import load_pelz_dataset, load_kupke, load_full_alnaji2021, load_short_reads, load_alnaji_excel

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

def select_classifier(clf_name: str)-> object:
    '''
        Selects a scikit-learn classifier by a given name. Is implemented in an
        extra function to use the same parameters in each usage of one of the 
        classifiers.
        :param clf_name: name of the classifier

        :return: Selected classifier as class implemented in scikit-learn
    '''
    if clf_name == "logistic_regression":
        clf = LogisticRegression(max_iter=4000)
    elif clf_name == "svc":
       clf = SVC(gamma=2, C=1)
    elif clf_name == "random_forest":
        clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
    else:
        print(f"classifier {clf_name} unknown!")
        exit()
    return clf

def set_labels(df: object, style: str, n_bins: int, labels: list=[])-> object:
    '''
        Sets the labels for the classifer. Can be done by using pd.cut() or by
        using the median/33-percentil as split.
        :param df: data frame including the data
        :param style: declares how to create the labels
        :param n_bins: number of bins to use
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

def select_datasets(df, dataset_name: str, features: list, n_bins: int)-> (object, object, object, object):
    '''
        Selects training a test data by a given name.
        :param df: pandas data frame including all data sets and features
        :param dataset_name: string indicating which data sets to include
        :param features: list with all features, that should be selected

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

    labeling_style = "median"

    t_df = df.loc[df["dataset_name"].isin(train)].copy().reset_index()
    X = t_df[features]
    y = set_labels(t_df, labeling_style, n_bins, labels)

    if len(val) != 0:
        v_df = df.loc[df["dataset_name"].isin(val)].copy().reset_index()
        X_val = v_df[features]
        y_val = set_labels(v_df, labeling_style, n_bins, labels)
    else:
        X_val = pd.DataFrame()
        y_val = pd.DataFrame()

    return X, y, X_val, y_val

