'''
    Analyzing the datasets to check their distribution etc
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH
from ml_utils import load_all_sets, set_labels

def check_distributions(df: object)-> None:
    '''
        Checks the distribution of the NGS counts for all given datasets.
        Creates a histogram for that. A second one is using the log() function.
        :param df: data frame including all data sets

        :return: None
    '''
    dataset_names = df["dataset_name"].unique().tolist()
    dataset_names.append("all")
    for name in dataset_names:
        if name == "all":
            n_df = df
        else:
            n_df = df[df["dataset_name"] == name]

        fig, axs = plt.subplots(2, 1, figsize=(4, 7), tight_layout=True)

        axs[0].hist(n_df["NGS_read_count"], bins=100)
        axs[0].set_xlabel("NGS count")
        axs[0].set_ylabel("# occurrences")
        axs[0].set_title(f"NGS distribution for {name}")

        axs[1].hist(n_df["NGS_log_norm"], bins=100)
        axs[1].set_xlabel("log(NGS count)")
        axs[1].set_ylabel("# occurrences")
        axs[1].set_title(f"log(NGS) distribution for {name}")

        save_path = os.path.join(RESULTSPATH, "ML", f"distribution_{name}.png")
        plt.savefig(save_path)
        plt.close()


def check_stat_parameters(df: object)-> None:
    '''
        Checks statistics for all data sets (like mean, max, min, etc). Creates
        a LaTeX table as a result and saves it.
        :param df: data frame including all data sets

        :return: None
    '''
    d = dict()
    g = df.groupby("dataset_name")
    d["Datasetname"] = g.size().index
    d["N_samples"] = g.size().tolist()
    d["Max"] = g.max("NGS_read_count")["NGS_read_count"].tolist()
    d["Min"] = g.min("NGS_read_count")["NGS_read_count"].tolist()
    d["Mean"] = g.mean("NGS_read_count")["NGS_read_count"].tolist()
    d["Median"] = g.median("NGS_read_count")["NGS_read_count"].tolist()
    r_df = pd.DataFrame(d)

    all_l = ["all", len(df.index),
             df["NGS_read_count"].max(), df["NGS_read_count"].min(),
             df["NGS_read_count"].mean(), df["NGS_read_count"].median()]
    r_df.loc[len(r_df.index)] = all_l
    r_df = r_df.set_index("Datasetname")
    
    path = os.path.join(RESULTSPATH, "ML", "stat_param_datasets.tex")
    r_df.to_latex(path, index=False, float_format="%.2f", longtable=True)


def check_duplicates(df: object)-> None:
    '''
        For the Segment, Start, End combinations that occur more than once this
        function checks how many of them have the same label
        :param df: data frame including all data points

        :return: None
    '''
    df["DI"] = df["Segment"] + "_" + df["Start"].map(str) + "_" + df["End"].map(str)
    entries = df["DI"].tolist()
    dupl = [e for e in entries if entries.count(e) > 1]

    df["class"] = set_labels(df, "median", 3)

    eq = 0
    a = 0
    for d in dupl:
        l_list = df.loc[df["DI"] == d]["class"].tolist()
        if len(l_list) == 2:
            a = a + 1
            if l_list[0] == l_list[1]:
                eq = eq + 1

    print(eq/a)


if __name__ == "__main__":
    all_df = load_all_sets()
#    check_distributions(all_df)
 #   check_stat_parameters(all_df)
    check_duplicates(all_df)
