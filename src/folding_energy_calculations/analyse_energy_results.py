"""
This script analyses the results that were created by Vienna RNA package.
It loads the .fold files and compares them.
"""
import os
import re
import sys
import shutil

import pandas as pd
import matplotlib.pyplot as plt

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
    fig, ax = plt.subplots(1, 1, figsize=(10,10), tight_layout=True)
    for s in SEGMENTS:
        s_df = df[df["segment"] == s]
        ax.scatter(s_df["delta_G"], s_df["length"], label=s)

    ax.set_title(f"correlation of delta G to sequence length for {d_set} data set")
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
        ax.scatter(x, s_df["NGS_read_count"], label=s)

    plt.locator_params(axis="y", nbins=10)
    n = "(normalized by sequence length)" if normalize else ""
    ax.set_title(f"correlation of delta G to NGS count {n}")
    ax.set_xlabel("delta G")
    ax.set_ylabel("NGS count")
    ax.legend()

    n = "_normalized" if normalize else ""
    save_path = os.path.join(path, f"deltaG_NGScount_{n}.pdf")
    plt.savefig(save_path)
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
    ax.set_ylabel("sequence length")
    ax.set_zlabel("NGS count")
    plt.show()
    '''
