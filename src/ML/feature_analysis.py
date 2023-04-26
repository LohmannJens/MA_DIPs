'''
    Analyzing the different features used in the classification.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH
from ml_utils import load_all_sets, junction_site_ohe, get_direct_repeat_length, get_3_to_5_ratio, get_length_proportion

def split_data(df: pd.DataFrame,
               feature: str,
               thresh: str="",
               col: str=""
               )-> (pd.DataFrame, pd.DataFrame):
    '''
        Splits a data frame by a given feature and threshold.
        :param df: pandas data frame
        :param feature: name of the feature to split by
        :param thresh: threshold indicating the point where to split the df
        :param col: name of the column to split by

        return: tuple containing two data frames:
                    data points above or equal to threshold
                    data points below threshold
    '''
    if col == "": col = feature

    if feature == "junction_site":
        pos = col.split("_")[0]
        df, cols = junction_site_ohe(df, pos)
    elif feature == "Direct_repeat":
        df[col] = df.apply(get_direct_repeat_length, axis=1)
    elif feature == "3_5_ratio":
        df[col] = df.apply(get_3_to_5_ratio, axis=1)
        thresh = df[col].median()
    elif feature == "length_proportion":
        df[col] = df.apply(get_length_proportion, axis=1)
        thresh = df[col].median()

    df1 = df[df[col] >= thresh]
    df2 = df[df[col] < thresh]

    return df1, df2


def test_feature(df: pd.DataFrame,
                 feature: str,
                 thresh: str="",
                 col: str=""
                 )-> None:
    '''
        Test for a given feature how the data is split up applying a given
        threshold.
        :param df: pandas data frame
        :param feature: name of the feature to split by
        :param thresh: threshold indicating the point where to split the df
        :param col: name of the column to split by

        :return: None
    '''
    df1, df2 = split_data(df, feature, thresh, col)

    if feature in ["3_5_ratio", "length_proportion"]:
        col = pd.concat([df1, df2])[feature].median()
        col = f"{col:.2}"

    fig, axs = plt.subplots(1, 1, figsize=(7, 7), tight_layout=True)

    axs.hist(df1["NGS_log_norm"], bins=100, alpha=0.3, label=col, color="b")
    axs.hist(df2["NGS_log_norm"], bins=100, alpha=0.3, label="other", color="y")
    axs.axvline(df1["NGS_log_norm"].mean(), color="b")
    axs.axvline(df2["NGS_log_norm"].mean(), color="y")
    axs.set_xlabel("log(NGS count)")
    axs.set_ylabel("# occurrences")
    axs.set_title(f"log(NGS) distribution for {feature}")
    axs.legend()

    save_path = os.path.join(RESULTSPATH, "ML", f"check_{feature}_feature.png")
    plt.rc("font", size=14)
    plt.savefig(save_path)
    plt.close()


if __name__ == "__main__":
    all_df = load_all_sets()
    features = [("junction_site", 1, "Start_5_A"),
                ("Direct_repeat", 2, ""),
                ("3_5_ratio", "", ""),
                ("length_proportion", "", "")]

    for f in features:
        test_feature(all_df, f[0], f[1], f[2])

