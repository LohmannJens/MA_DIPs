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
    use_abs = True
    if use_abs:
        file_path = os.path.join(DATAPATH, "Pelz2021", "rates_for_simulation.xlsx")
    else:
        file_path = os.path.join(DATAPATH, "Pelz2021", "rates_for_simulation_frac.xlsx")
    
    data_dict = pd.read_excel(io=file_path,
                              sheet_name=None,
                              header=0,
                              na_values=["", "None"],
                              keep_default_na=False,
                              converters={"Start": int,"End": int})

    return data_dict["PR8"]


def starting_conditions(df: pd.DataFrame)-> Tuple[int, int]:
    '''
    
    '''
    #print(df.shape[0]) this is max number of DI RNA candidates

    g = df.groupby(["class"]).count()

    n_gain = g.iloc[2,2]
    n_loss = g.iloc[3,3]

    return n_gain, n_loss


def calculate_probabilities(df: pd.DataFrame)-> None:
    '''
    
    '''
    fig, ax = plt.subplots(1,1)
    prob = dict()

    for c in ["gain", "loss", "de novo gain", "de novo loss"]:
        increase = list()
        decrease = list()
        p_df = df[df["class"] == c]

        for i in range(3, p_df.shape[1]-2):
            series = p_df.iloc[:,i] - p_df.iloc[:,i+1]

            increase.append(series[series < 0].count())
            decrease.append(series[series > 0].count())

        p_list = [a/(a+b) for a, b in zip(increase, decrease)]
        ax.plot(p_list, label = c)

        prob[c] = sum(increase) / (sum(increase) + sum(decrease))
        print(f"{c}:\t{prob[c]}")

    print(f"gain\t{(prob['gain'] + prob['de novo gain'])/2}")
    print(f"loss\t{(prob['loss'] + prob['de novo loss'])/2}")
    
    ax.set_xlabel("timepoint")
    ax.set_ylabel("prob. of increase")
    ax.set_ylim(bottom=0, top=1.0)
    ax.legend()

    plt.show()


def calculate_rates(df: pd.DataFrame)-> None:
    '''
    
    '''
    fig, ax = plt.subplots(1,1)
    inc_d = dict()
    dec_d = dict()

    for c in ["gain", "loss", "de novo gain", "de novo loss"]:
        increase = list()
        decrease = list()
        p_df = df[df["class"] == c]

        for i in range(3, p_df.shape[1]-2):
            series = (p_df.iloc[:,i] - p_df.iloc[:,i+1]) / p_df.iloc[:,i] # calculate percentage

            increase.append(series[series < 0].replace(-np.inf, np.nan).mean())
            decrease.append(series[series > 0].mean())

        increase = [x for x in increase if str(x) != "nan"]
        decrease = [x for x in decrease if str(x) != "nan"]

        ax.plot(increase, label = f"{c} inc.")
        ax.plot(decrease, label = f"{c} dec.")

        inc_d[c] = sum(increase) / len(increase)
        dec_d[c] = sum(decrease) / len(decrease)

        print(f"{c} increase:\t{inc_d[c]}")
        print(f"{c} decrease:\t{dec_d[c]}")

    print(f"gain inc\t{(inc_d['gain'] + inc_d['de novo gain'])/2}")
    print(f"gain dec\t{(dec_d['gain'] + dec_d['de novo gain'])/2}")

    print(f"loss inc\t{(inc_d['loss'] + inc_d['de novo loss'])/2}")
    print(f"loss dec\t{(dec_d['loss'] + dec_d['de novo loss'])/2}")

    ax.set_xlabel("timepoint")
    ax.set_ylabel("rel. increase")
    ax.legend()

    plt.show()


def calculate_de_novo_events(df: pd.DataFrame):
    '''
    
    '''
    for c in ["gain", "loss"]:
        gain_ratios = list()
        loss_ratios = list()
        p_df = df[df["class"].isin([c, f"de novo {c}"])]

        for i in range(3, p_df.shape[1]-2):
            start = p_df.iloc[:,i]
            s_candidates = start[start != 0].count()

            end = p_df.iloc[:,i+1]
            e_candidates = end[end != 0].count()

            d_candidates = e_candidates - s_candidates

            if d_candidates > 0:
                gain_ratios.append(d_candidates/s_candidates)
            else:
                loss_ratios.append(d_candidates/s_candidates)

        print(f"{c}\t{sum(gain_ratios[1:])/(len(gain_ratios) + len(loss_ratios) - 1)}")
        print(f"{c}\t{sum(loss_ratios)/(len(gain_ratios) + len(loss_ratios) - 1)}")


def calculate_initialising_range(df: pd.DataFrame):
    '''
    
    '''
    for i in range(3, df.shape[1]-2):
        a_df = df[(df.iloc[:,i] == 0) & (df.iloc[:,i+1] != 0)]
        print(a_df.iloc[:,i+1].mean())


def check_label_developement_to_gt(df: pd.DataFrame):
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
    

def check_label_distribution_over_time(df: pd.DataFrame):
    '''
    
    '''
    only_top = False
    if only_top: # filter for top candidates
        file_path = os.path.join(DATAPATH, "Pelz2021", "top_candidates_ranked.xlsx")
        top_df = pd.read_excel(io=file_path,
                                  sheet_name=None,
                                  header=0,
                                  na_values=["", "None"],
                                  keep_default_na=False)["PR8"]

        all_candidates = top_df["DI"].tolist()

        top_candidates = ["PB1_139_2056", "PA_164_2028", "PA_138_1948", "PA_244_2074", "PA_163_1990",
            "PA_137_1916", "PB1_113_1897", "PB2_206_2152", "PB2_269_2202", "PA_124_1940", "PB1_218_2091"]

        n = 200
        top_candidates = all_candidates[0:n]
        print(top_candidates)
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


def ground_truth_population(df):
    '''
    
    '''
    d = dict({"gain": list(), "loss": list()})
    
    for i in range(3, df.shape[1]-1):
        d["gain"].append(df[df["class"].isin(["gain", "de novo gain"]) & df.iloc[:,i] != 0].count()[0])
        d["loss"].append(df[df["class"].isin(["loss", "de novo loss"]) & df.iloc[:,i] != 0].count()[0])

    fig, ax = plt.subplots(1,1)
    for k, v in d.items():
        ax.plot(v, label=k)
    ax.set_xlabel("time point")
    ax.set_ylabel("number of candidates")
    ax.set_ylim(bottom=0)
    plt.legend()
    plt.show()
    

if __name__ == "__main__":
    df = load_excel()
    '''
    n_gain, n_loss = starting_conditions(df)
    
    calculate_probabilities(df)
    calculate_rates(df)

    calculate_de_novo_events(df)
    '''
    #check_label_distribution_over_time(df)
    #check_label_developement_to_gt(df)

    calculate_initialising_range(df)

    #ground_truth_population(df)