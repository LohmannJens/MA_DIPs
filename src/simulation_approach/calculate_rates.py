'''

'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import Tuple

sys.path.insert(0, "..")
from utils import DATAPATH


def load_excel()-> pd.DataFrame:
    '''
    
    '''
    file_path = os.path.join(DATAPATH, "Pelz2021", "rates_for_simulation_frac.xlsx")
    data_dict = pd.read_excel(io=file_path,
                              sheet_name=None,
                              header=0,
                              na_values=["", "None"],
                              keep_default_na=False,
                              converters={"Start": int,"End": int})

    return data_dict["PR8"]


def calculate_probabilities(df: pd.DataFrame)-> float:
    '''
    
    '''
    fig, ax = plt.subplots(1,1)

    for c in ["gain", "loss"]:
        increase = list()
        decrease = list()
        p_df = df[df["class"] == c]

        for i in range(3, p_df.shape[1]-2):
            series = p_df.iloc[:,i] - p_df.iloc[:,i+1]

            increase.append(series[series < 0].count())
            decrease.append(series[series > 0].count())

        p_list = [a/(a+b) for a, b in zip(increase, decrease)]
        ax.plot(p_list, label = c)

        probability = sum(increase) / (sum(increase) + sum(decrease))

        print(f"{c}:/t{probability}")
    
    ax.set_xlabel("timepoint")
    ax.set_ylabel("prob. of increase")
    ax.set_ylim(bottom=0, top=1.0)
    ax.legend()

    plt.show()
    
    return sum(increase) / (sum(increase) + sum(decrease))


def starting_conditions(df)-> Tuple[int, int]:
    '''
    
    '''
    g = df.groupby(["class"]).count()

    n_gain = g.iloc[2,2]
    n_loss = g.iloc[3,3]

    return n_gain, n_loss


def check_label_developement_to_gt(df):
    '''
    
    '''
    differ_list = list()

    ground_truth_class = df["class"]
    
    start = df.iloc[:,3]

    for i in range(4, df.shape[1]-1):
        class_series = pd.Series(np.where((start == 0.0), np.where((start < df.iloc[:,i]), "de novo gain", "de novo loss"), np.where((start < df.iloc[:,i]), "gain", "loss")))
        differ_list.append(ground_truth_class.compare(class_series).shape[0])

    fig, ax = plt.subplots(1,1)
    ax.plot(differ_list)
    ax.set_xlabel("timepoint")
    ax.set_ylabel("no. of different labels")
    ax.set_ylim(bottom=0)
    
    plt.show()
    

def check_label_distribution_over_time(df):
    '''
    
    '''
    only_top = False
    if only_top: # filter for top 5 candidates
        top_candidates = ["PB1_139_2056", "PA_164_2028", "PA_138_1948", "PA_244_2074", "PA_163_1990",
            "PA_137_1916", "PB1_113_1897", "PB2_206_2152", "PB2_269_2202", "PA_124_1940", "PB1_218_2091"]
        df["DI"] = df["Segment"] + "_" + df["Start"].map(str) + "_" + df["End"].map(str)
        df = df[df["DI"].isin(top_candidates)]
        df.drop(["DI"], axis=1, inplace=True)

    data_dict = dict({"gain": list(), "loss": list(), "de novo gain" : list(), "de novo loss": list()})
    start = df.iloc[:,3]
    for i in range(4, df.shape[1]-1):
        labels_t = np.where((start == 0), np.where((start < df.iloc[:,i]), "de novo gain", "de novo loss"), np.where((start < df.iloc[:,i]), "gain", "loss"))
        for l in ["gain", "loss", "de novo gain", "de novo loss"]:
            data_dict[l].append(np.count_nonzero(labels_t == l))

    fig, ax = plt.subplots(1,1)
    for l in ["gain", "loss", "de novo gain", "de novo loss"]:
        ax.plot(data_dict[l], label=l)
    ax.set_xlabel("timepoint")
    ax.set_ylabel("count of labels")
    ax.set_ylim(bottom=0)
    ax.legend()

    plt.show()


if __name__ == "__main__":
    df = load_excel()

    n_gain, n_loss = starting_conditions(df)

    calculate_probabilities(df)
    check_label_distribution_over_time(df)
    check_label_developement_to_gt(df)