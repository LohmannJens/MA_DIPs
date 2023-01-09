'''
    Does a linear regression for each timepoint of the pelz data set and
    compares the change of the slope over time. Includes only the de novo
    candidates
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import get_seq_len
sys.path.insert(0, "../regression_length_vs_occurrence")
from regression_length_occurrence import format_dataset_for_plotting


def assign_label(row)-> str:
    '''
        Give the different de novo and gain/loss labels to a row of the pelz
        data set.
        :param row: a single row of the data set

        :return: string indicating the label/class
    '''
    start = row["VB3-Saat"]
    end = row["VB3-48"]

    if (start == 0 and end == 0):
        return "de_novo_loss"
    elif (start == 0 and end !=0):
        return "de_novo_gain"
    elif (start != 0 and start >= end):
        return "loss"
    elif (start != 0 and start < end):
        return "gain"
    else:
        print(start, end)
        return "unknown"


if __name__ == "__main__":
    path = os.path.join(DATAPATH, "Pelz2021", "pelz_split_by_timepoints.xlsx")
    df = pd.read_excel(path, na_values="", keep_default_na=False)
    df["class"] = df.apply(assign_label, axis=1)

    t_x = [4.1, 4.6, 5.0, 7.0, 7.6, 8.1, 8.6, 9.1, 11.5, 12.6, 13.0, 16.0, 16.6, 17.1, 19.6, 20.6, 21.0, 21.6, 23.1, 23.6, 24.0, 24.6]
    m = list()
    for i in np.arange(1, 23):
        t_df = df.iloc[:, [0, 1, 2, i+3, -1]].copy()
        t_df = t_df.loc[t_df["class"].isin(["de_novo_loss", "de_novo_gain"])]
        t_df.columns = ["Segment", "Start", "End", "NGS_read_count", "class"]

        x, y, _, _ = format_dataset_for_plotting(t_df, "PR8")
        model = LinearRegression().fit(x.reshape((-1, 1)), y)
        m.append(model.coef_[0])

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)
    ax.plot(t_x, m)
    ax.set_xlabel("number of measurement")
    ax.set_ylabel("slope of regression")

    save_path = os.path.join(RESULTSPATH, "PR8", "Pelz_slope_over_time.png")
    plt.savefig(save_path)
    plt.close()
