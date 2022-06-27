'''
    Plots the results of the sliding window approach conducted in 
    "run_sliding_window_approach.py". Plots the estimated delta G together with
    the start and end positions of the junction sites to show if there is a
    relation between them.
'''
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, load_excel, load_short_reads


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
        save_file = os.path.join(save_path, f"{k}_sliding_window_{w_s}_{s_s}.pdf")
        fig.savefig(save_file)
        plt.close()


if __name__ == "__main__":
    filepath = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    short_reads_filepath = os.path.join(DATAPATH, "alnaji2019", "Small_deletionSize_FA.xlsx")    
    cleaned_data_dict = load_excel(filepath)
    all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)


    window_size = 20
    step_size = 1
    plot_deletions_with_delta_G(all_reads_dict, window_size, step_size)

