'''
    Plots the results of the sliding window approach conducted in 
    "run_sliding_window_approach.py". Plots the estimated delta G together with
    the start and end positions of the junction sites to show if there is a
    relation between them.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import load_alnaji_excel, load_short_reads, get_stat_symbol


def slice_dataset(seg_df, energy_df, seg):
    '''

    '''
    slice_index = dict({"PB1": 700, "PB2": 500,
                  "PA": 700, "HA": 500,
                  "NP": 400, "NA": 600,
                  "M": 400, "NS": 500})
    seg_df = seg_df[seg_df.index < slice_index[seg]]
    energy_df = energy_df[energy_df["position"] < slice_index[seg]]

    return seg_df, energy_df


def plot_deletions_with_delta_G(d: dict, w_s: int, s_s: int)-> None:
    '''
        Loads the deletion start and end points and the estimated delta G and
        Plots them together in one plot.
        :param d: dictionary with the start and end point data split by strains
        :param w_s: size of the window in the sliding window approach
        :param s_s: step size of the sliding window approach

        :return: None
    '''
    energy_path = os.path.join(DATAPATH, "energy_calculation", "sliding_window")
    for k, v in d.items():
        fig, axs = plt.subplots(8, 1, figsize=(10, 15), tight_layout=True)

        for i, s in enumerate(SEGMENTS):
            seg_df = v.loc[v["Segment"] == s]
            if seg_df.size == 0:
                continue

            start_seg_df = seg_df.loc[:, ["Start", "NGS_read_count"]]
            start_seg_df.rename(columns={"Start": "Position"}, inplace=True)
            end_seg_df = seg_df.loc[:, ["End", "NGS_read_count"]]
            end_seg_df.rename(columns={"End": "Position"}, inplace=True)
            concat_seg_df = pd.concat([start_seg_df, end_seg_df])
            concat_seg_df = concat_seg_df.groupby("Position").sum()
            
            energy_file = os.path.join(energy_path, f"{k}_{s}_{w_s}_{s_s}.csv")
            energy_df = pd.read_csv(energy_file)
            
            concat_seg_df, energy_df = slice_dataset(concat_seg_df, energy_df, s)

            l1 = axs[i].twinx().bar(concat_seg_df.index, concat_seg_df["NGS_read_count"])
            l2, = axs[i].plot(energy_df["position"], energy_df["delta_G"], color="green")
            l3 = axs[i].axhline(y=energy_df["delta_G"].mean(), c="r", linestyle="--")

            axs[i].set_title(s)
            axs[i].set_xlim(left=0.0, right=max(energy_df["position"]))
            axs[i].set_ylim(bottom=0.0, top=min(energy_df["delta_G"]))
            axs[i].set_xlabel("Sequence position")
            axs[i].set_ylabel("\u0394 G")

        fig.legend([l1, l2, l3], ["NGS count", "\u0394 G", "mean of \u0394 G"])
        fig.suptitle(k, ha="left")

        save_path = os.path.join(RESULTSPATH, "free_energy_estimations")
        save_file = os.path.join(save_path, f"{k}_sliding_window_{w_s}_{s_s}_sliced.pdf")
        fig.savefig(save_file)
        plt.close()


def create_boxplots(d: dict, w_s: int, s_s: int)-> None:
    '''
        Loads the deletion start and end points and the estimated delta G and
        Plots them together in one plot.
        :param d: dictionary with the start and end point data split by strains
        :param w_s: size of the window in the sliding window approach
        :param s_s: step size of the sliding window approach

        :return: None
    '''
    energy_path = os.path.join(DATAPATH, "energy_calculation", "sliding_window")
    fig, axs = plt.subplots(4, 1, figsize=(10, 15), tight_layout=True)
    for i, (k, v) in enumerate(d.items()):
        data = list()
        annotations = list()
        y_min = 0.0

        for idx, s in enumerate(SEGMENTS):
            seg_df = v.loc[v["Segment"] == s]
            if seg_df.size != 0:
                start_seg_df = seg_df.loc[:, ["Start", "NGS_read_count"]]
                start_seg_df.rename(columns={"Start": "Position"}, inplace=True)
                end_seg_df = seg_df.loc[:, ["End", "NGS_read_count"]]
                end_seg_df.rename(columns={"End": "Position"}, inplace=True)
                concat_seg_df = pd.concat([start_seg_df, end_seg_df])
                concat_seg_df = concat_seg_df.groupby("Position").sum()

                energy_file = os.path.join(energy_path, f"{k}_{s}_{w_s}_{s_s}.csv")
                energy_df = pd.read_csv(energy_file)

                concat_seg_df, energy_df = slice_dataset(concat_seg_df, energy_df, s)

                if y_min > min(energy_df["delta_G"]):
                    y_min = min(energy_df["delta_G"]) + 1
                
                df1 = energy_df[energy_df["position"].isin(list(concat_seg_df.index))]
                df2 = energy_df[~energy_df["position"].isin(list(concat_seg_df.index))]
                data.extend([list(df1["delta_G"]), list(df2["delta_G"])])

                # do statistics
                # t-test is not suitable, because no normal deviation is given
                res = stats.mannwhitneyu(df1["delta_G"], df2["delta_G"])
                symbol = get_stat_symbol(res.pvalue)
            else:
                symbol = ""
                data.extend([0, 0])

            annotations.append(axs[i].annotate(f"{s} {symbol}", (idx*2 +1.5, 0.0), ha="center", size="small"))

        for idx, text in enumerate(annotations):
            text.set_position((idx*2+1.5, y_min))

        box = axs[i].boxplot(data, labels=["has del.", "no del."]*8)
        axs[i].set_title(k)
        axs[i].set_xlabel("Segments")
        axs[i].set_ylabel("\u0394 G")
        axs[i].set_ylim(bottom=y_min)
        axs[i].set_ylim(top=20.0)

    save_path = os.path.join(RESULTSPATH, "free_energy_estimations")
    save_file = os.path.join(save_path, f"boxplot_sliding_window_{w_s}_{s_s}_sliced.pdf")
    fig.savefig(save_file)
    plt.close()


def create_scatterplots(d: dict, w_s: int, s_s: int)-> None:
    '''
        Loads the deletion start and end points and the estimated delta G and
        Plots the delta G against the NGS count in a scatter plot. Also does a
        linear regression to detect a possible correlation.
        :param d: dictionary with the start and end point data split by strains
        :param w_s: size of the window in the sliding window approach
        :param s_s: step size of the sliding window approach

        :return: None
    '''
    energy_path = os.path.join(DATAPATH, "energy_calculation", "sliding_window")
    fig, axs = plt.subplots(4, 1, figsize=(10, 15), tight_layout=True)
    for i, (k, v) in enumerate(d.items()):
        data = list()

        for s in SEGMENTS:
            seg_df = v.loc[v["Segment"] == s]
            if seg_df.size != 0:
                start_seg_df = seg_df.loc[:, ["Start", "NGS_read_count"]]
                start_seg_df.rename(columns={"Start": "Position"}, inplace=True)
                end_seg_df = seg_df.loc[:, ["End", "NGS_read_count"]]
                end_seg_df.rename(columns={"End": "Position"}, inplace=True)
                concat_seg_df = pd.concat([start_seg_df, end_seg_df])
                concat_seg_df = concat_seg_df.groupby("Position").sum()

                energy_file = os.path.join(energy_path, f"{k}_{s}_{w_s}_{s_s}.csv")
                energy_df = pd.read_csv(energy_file)

                concat_seg_df["position"] = concat_seg_df.index
                merged_df = pd.merge(concat_seg_df, energy_df)
                data.append(merged_df)

        data = pd.concat(data, ignore_index=True)

        axs[i].scatter(data["NGS_read_count"], data["delta_G"])

        b, a, r, p, std_err = stats.linregress(data["NGS_read_count"], data["delta_G"], )
        axs[i].plot(data["NGS_read_count"], a + b * data["delta_G"])

        axs[i].set_title(f"{k} R²:{r:.2}")
        axs[i].set_xlabel("Segments")
        axs[i].set_ylabel("\u0394 G")
        axs[i].invert_yaxis()


    save_path = os.path.join(RESULTSPATH, "free_energy_estimations")
    save_file = os.path.join(save_path, f"boxplot_sliding_window_scatter.pdf")
    fig.savefig(save_file)
    plt.close()


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)

    window_size = 1
    step_size = 1
    plot_deletions_with_delta_G(all_reads_dict, window_size, step_size)
    create_boxplots(all_reads_dict, window_size, step_size)
    create_scatterplots(all_reads_dict, window_size, step_size)

