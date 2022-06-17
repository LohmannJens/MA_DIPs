"""
This script analyses the results that were created by Vienna RNA package.
It loads the .fold files and analyses the delta G in different ways.
"""
import os
import re
import sys
import shutil

import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats
from mpl_toolkits.mplot3d import Axes3D

sys.path.insert(0, "..")
sys.path.insert(0, "../density_and_length_analysis")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, get_seq_len


def extract_data_from_file(f: str)-> (int, float):
    '''
        Opens a .fold file form the free energy calculations and gets the 
        length of the sequence and the delta G
        :param f: filename

        :return: Tuple with length of sequence and delta G
    '''
    # open file
    with open(f, "r") as handle:
        # read header from first line
        l1 = handle.readline()

        # get sequence length from second line (is length of line)
        l2 = handle.readline()
        length = len(l2) - 1

        # get delta G from third line
        l3 = handle.readline()
        match = re.findall("\((\S*)\)$", l3)
        delta_G = float(match[0])

    return length, delta_G


def structure_dataframe(path: str)-> object:
    '''
        Gets a path to a folder with .fold files. Opens each of them and
        extracts the lenght and delta G of the sequences. Then a data frame
        is created including this and more information.
        :param path: directory with .fold files inside

        :return: pandas data frame
    '''
    strains = list()
    segments = list()
    lengths = list()
    delta_Gs = list()
    delta_Gs_shuffled = list()
    delta_Gs_random = list()
    starts = list()
    ends = list()
    NGSs = list()

    # loop over .fold files
    for f in os.listdir(path):
        fname = os.fsdecode(f)
        if fname.endswith(".fold"):
            # remove ".fold" and split up filename to get: strain, segment, (start), (end)
            f_split = f[:-5].split(sep="_")
            strains.append(f_split[0])
            segments.append(f_split[1])
            if len(f_split) > 2:
                starts.append(int(f_split[2]))
                ends.append(int(f_split[3]))
                NGSs.append(int(f_split[4]))

                path_shuffled = os.path.join(f"{path}_shuffled",
                                             f"{f[:-5]}_shuffled.fold")
                _, delta_G_shuffled = extract_data_from_file(path_shuffled)
                delta_Gs_shuffled.append(delta_G_shuffled)

                path_random = os.path.join(f"{path}_randomcrop",
                                           f"{f[:-5]}_randomcrop.fold")
                _, delta_G_random = extract_data_from_file(path_random)
                delta_Gs_random.append(delta_G_random)

            length, delta_G = extract_data_from_file(os.path.join(path, f))
            lengths.append(length)
            delta_Gs.append(delta_G)


    # convert dict to data frame
    d = dict()
    d["strain"] = strains
    d["segment"] = segments
    d["length"] = lengths
    d["delta_G"] = delta_Gs
    if (len(starts) != 0 and len(ends) != 0):
        d["start"] = starts
        d["end"] = ends
        d["NGS_read_count"] = NGSs
        d["delta_G_shuffled"] = delta_Gs_shuffled
        d["delta_G_random"] = delta_Gs_random

    df = pd.DataFrame.from_dict(d)
    return df


def plot_deltaG_length(df: object, path: str, d_set: str)-> None:
    '''
        Plots delta G against the sequence length as a scatter plot
        :param df: data frame with the data
        :param path: path to the results folder
        :param d_set: is "full" or "cropped" depending on which set was used

        :return: None
    '''
    r, p = stats.pearsonr(df["delta_G"], df["length"])

    fig, ax = plt.subplots(1, 1, figsize=(10,10), tight_layout=True)
    for s in SEGMENTS:
        s_df = df[df["segment"] == s]
        ax.scatter(s_df["delta_G"], s_df["length"], label=s, alpha=0.3)

    ax.set_title(f"correlation of delta G to sequence length for {d_set} data set (r = {r:.2})")
    ax.set_xlabel("delta G")
    ax.set_ylabel("sequence length")
    ax.legend()

    save_path = os.path.join(path, f"deltaG_length_{d_set}.pdf")
    plt.savefig(save_path)
    plt.close()


def plot_deltaG_NGS(df: object, path: str, normalize: bool)-> None:
    '''
        Plots delta G against NGS count as a scatter plot.
        :param df: data frame with the data
        :param path: path to the results folder
        :param normalize: indicates if delta G should be normalized by sequence
                          length

        :return: None
    '''
    fig, ax = plt.subplots(1, 1, figsize=(10,10), tight_layout=True)
    for s in SEGMENTS:
        s_df = df[cropped_df["segment"] == s]
        x = s_df["delta_G"]/s_df["length"] if normalize else s_df["delta_G"]
        ax.scatter(x, s_df["NGS_read_count"], label=s, alpha=0.3)

    plt.locator_params(axis="y", nbins=10)
    n = "(normalized by sequence length)" if normalize else ""
    ax.set_title(f"correlation of delta G to NGS count {n}")
    ax.set_xlabel("delta G")
    ax.set_ylabel("NGS count")
    ax.legend()

    n = "_normalized" if normalize else ""
    save_path = os.path.join(path, f"deltaG_NGScount{n}.pdf")
    plt.savefig(save_path)
    plt.close()


def plot_delta_G_observed_expected(df: object, path: str, mode: str)-> None:
    '''
        Plots delta G against delta G of randomly shuffled or randomly cut
        sequences.
        :param df: data frame with the data
        :param path: path to the results folder
        :param mode: is either 'shuffled' or 'random'

        :return: None

    '''
    fig, axs = plt.subplots(2, 2, figsize=(10,10), tight_layout=True)
    j = 0
    for i, strain in enumerate(["Cal07", "NC", "Perth", "BLEE"]):
        strain_df = df[df["strain"] == strain]
        if i == 2:
            j = 1
        for s in SEGMENTS:
            seg_df = strain_df[strain_df["segment"] == s]
            axs[i%2][j].scatter(seg_df["delta_G"], seg_df[f"delta_G_{mode}"], label=s, alpha=0.3)

        r, p = stats.pearsonr(strain_df["delta_G"], strain_df["length"])
        
        axs[i%2][j].set_title(f"delta G vs delta G {mode} {strain} (r = {r:.2})")
        axs[i%2][j].set_xlabel("delta G")
        axs[i%2][j].set_ylabel(f"delta G {mode}")
        axs[i%2][j].set_xlim([-650, 0])
        axs[i%2][j].set_ylim([-650, 0])
        axs[i%2][j].plot([0,1], [0,1], transform=axs[i%2][j].transAxes, c="grey", linestyle="--")
        axs[i%2][j].legend()

    plt.locator_params(axis="y", nbins=10)

    save_path = os.path.join(path, f"deltaG_observed_{mode}.pdf")
    plt.savefig(save_path)
    plt.close()


def create_difference_boxplots(df: object, path: str)-> None:
    '''
        Creates boxplots for the difference of the delta G against the expected
        values. Is done for the random and shuffled approach and saved in one
        figure. The boxes of the boxplot are split up by the segments.
        :param df: data frame with the values
        :param path: path to the results folder
        
        :return: None
    '''
    # calculate differences of delta G to random cropped data
    df["random_diff"] = df["delta_G"] - df["delta_G_random"]
 #   data = [df.loc[df["segment"] == s, "random_diff"] for s in SEGMENTS]

    fig, axs = plt.subplots(2, 2, figsize=(10,10), tight_layout=True)
    j = 0
    for i, strain in enumerate(["Cal07", "NC", "Perth", "BLEE"]):
        strain_df = df[df["strain"] == strain]
        data = [strain_df.loc[strain_df["segment"] == s, "random_diff"] for s in SEGMENTS]

        if i == 2:
            j = 1

        axs[i%2][j].boxplot(data, labels=SEGMENTS)
        axs[i%2][j].set_title(f"random approach {strain}")
        axs[i%2][j].set_xlabel("segments")
        axs[i%2][j].set_ylabel("\u0394\u0394G (\u0394G - \u0394G random approach)")
        axs[i%2][j].set_ylim([-60, 50])
        axs[i%2][j].axhline(y=0, c="r", linestyle="--", linewidth=0.5)

    save_path = os.path.join(path, "boxplot_delta_delta_G.pdf")
    fig.savefig(save_path)
    plt.close()


if __name__ == "__main__":
    path = os.path.join(DATAPATH, "energy_calculation")

    path_full = os.path.join(path, "full_sequences")
    path_cropped = os.path.join(path, "cropped_sequences")

    full_df = structure_dataframe(path_full)
    cropped_df = structure_dataframe(path_cropped)

    results_path = os.path.join(RESULTSPATH, "free_energy_estimations")
    
    plot_deltaG_length(full_df, results_path, "full")
    plot_deltaG_length(cropped_df, results_path, "cropped")

    plot_deltaG_NGS(cropped_df, results_path, False)
    plot_deltaG_NGS(cropped_df, results_path, True)
    
    plot_delta_G_observed_expected(cropped_df, results_path, "random")
    plot_delta_G_observed_expected(cropped_df, results_path, "shuffled")
    
    create_difference_boxplots(cropped_df, results_path)

    # This part creates a 3D plot of delta G, length and NGS count.
    # It basically just combines the two plots above.
    # Not sure if I will use it, but it will stay here for the moment
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    for s in SEGMENTS:
        df = cropped_df[cropped_df["segment"] == s]
        ax.scatter(df["delta_G"], df["length"], df["NGS_read_count"], label=s)

    ax.set_xlabel("delta G")
    ax.set_ylabel("length")
    ax.set_zlabel("NGS_read_count")
    plt.show()
    '''
