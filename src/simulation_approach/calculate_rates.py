'''

'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from typing import Tuple

sys.path.insert(0, "..")
from utils import DATAPATH, SEGMENTS, COLORS, NUCLEOTIDES, QUANT, N_SAMPLES


def load_excel()-> pd.DataFrame:
    '''
    
    '''
    file_path = os.path.join(DATAPATH, "Pelz2021", "rates_for_simulation.xlsx")
    data_dict = pd.read_excel(io=file_path,
                              sheet_name=None,
                              header=0,
                              na_values=["", "None"],
                              keep_default_na=False,
                              converters={"Start": int,"End": int})

    return data_dict["PR8"]


def calculate_probabilities(df: pd.DataFrame)-> Tuple[float, float]:
    '''
    
    '''
    increase = list()
    decrease = list()

    for i in range(3, df.shape[1]-2):
        series = df.iloc[:,i] - df.iloc[:,i+1]

        increase.append(series[series < 0].count())
        decrease.append(series[series > 0].count())

    """
    p_list = [a/(a+b) for a, b in zip(increase, decrease)]
    plt.plot(p_list)
    plt.xlabel("timepoint")
    plt.ylabel("prob. of increase")
    plt.show()
    """
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

  #  df.drop(df[df["class"] == "de novo gain"].index, inplace=True)
   # df.drop(df[df["class"] == "de novo loss"].index, inplace=True)
    #df.reset_index(drop=True, inplace=True)

    ground_truth_class = df["class"]
    ground_truth_class.replace({"de novo gain": "gain", "de novo loss": "loss"}, inplace=True)
    
    start_series = df.iloc[:,3]
    denovo_series = np.where(start_series == 0, "denovo", "no")

    for i in range(4, df.shape[1]-1):
        class_series = pd.Series(np.where((df.iloc[:,3] < df.iloc[:,i]), "gain", "loss"))
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

    """
    n_gain, n_loss = starting_conditions(df)

    prob_gain = calculate_probabilities(df[df["class"] == "gain"])
    print(prob_gain)
    print("#######################################")    
    prob_loss = calculate_probabilities(df[df["class"] == "loss"])
    print(prob_loss)
    """

    check_label_distribution_over_time(df)

 #   check_label_developement_to_gt(df)