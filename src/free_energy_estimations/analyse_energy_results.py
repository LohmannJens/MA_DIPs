'''
    This script analyses the results that were created by Vienna RNA package.
    It loads the .fold files and analyses the delta G in different ways.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, QUANT, N_SAMPLES
from utils import get_seq_len, load_alnaji_excel, load_short_reads, get_stat_symbol, generate_sampling_data


def extract_data_from_file(f: str)-> (int, float, str):
    '''
        Opens a .fold file form the free energy calculations and gets the 
        length of the sequence and the delta G
        :param f: filename

        :return: Tuple with length of sequence, delta G and secondary structure
                 in Dot-Bracket notation
    '''
    # open file
    with open(f, "r") as handle:
        # read header from first line
        l1 = handle.readline()

        # get sequence length from second line (is length of line)
        l2 = handle.readline()
        length = len(l2) - 1

        # get delta G and secondary structure from third line
        l3 = handle.readline()
        l3_splitted = l3.rstrip().split(" ")
        sec_structure = l3_splitted[0]
        delta_G = float(l3_splitted[1][1:-1])

    return length, delta_G, sec_structure


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
    sec_structures = list()
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
                _, delta_G_shuffled, _ = extract_data_from_file(path_shuffled)
                delta_Gs_shuffled.append(delta_G_shuffled)

                path_random = os.path.join(f"{path}_randomcrop",
                                           f"{f[:-5]}_randomcrop.fold")
                _, delta_G_random, _ = extract_data_from_file(path_random)
                delta_Gs_random.append(delta_G_random)

            length, delta_G, sec_struct = extract_data_from_file(os.path.join(path, f))
            lengths.append(length)
            delta_Gs.append(delta_G)
            sec_structures.append(sec_struct)

    # convert dict to data frame
    d = dict()
    d["strain"] = strains
    d["segment"] = segments
    d["length"] = lengths
    d["delta_G"] = delta_Gs
    d["sec_structure"] = sec_structures
    if (len(starts) != 0 and len(ends) != 0):
        d["start"] = starts
        d["end"] = ends
        d["NGS_read_count"] = NGSs
        d["delta_G_shuffled"] = delta_Gs_shuffled
        d["delta_G_random"] = delta_Gs_random

    df = pd.DataFrame.from_dict(d)
    return df


def plot_deltaG_length(df: pd.DataFrame,
                       path: str,
                       d_set: str
                       )-> None:
    '''
        Plots delta G against the sequence length as a scatter plot
        :param df: data frame with the data
        :param path: path to the results folder
        :param d_set: is "full" or "cropped" depending on which set was used

        :return: None
    '''
    plt.rc("font", size=14)
    r, p = stats.pearsonr(df["delta_G"], df["length"])

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)
    for s in SEGMENTS:
        s_df = df[df["segment"] == s]
        ax.scatter(s_df["delta_G"], s_df["length"], label=s, alpha=0.3)

    ax.set_title(f"{d_set} sequences (r = {r:.2})")
    ax.set_xlabel("\u0394 G")
    ax.set_ylabel("sequence length")
    ax.legend()

    save_path = os.path.join(path, f"deltaG_length_{d_set}.png")
    plt.savefig(save_path)
    plt.close()


def plot_deltaG_NGS(df: pd.DataFrame,
                    path: str,
                    normalize: bool
                    )-> None:
    '''
        Plots delta G against NGS count as a scatter plot.
        :param df: data frame with the data
        :param path: path to the results folder
        :param normalize: indicates if delta G should be normalized by sequence
                          length

        :return: None
    '''
    plt.rc("font", size=14)
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)
    for s in SEGMENTS:
        s_df = df[cropped_df["segment"] == s]
        x = s_df["delta_G"]/s_df["length"] if normalize else s_df["delta_G"]
        ax.scatter(x, s_df["NGS_read_count"], label=s, alpha=0.3)

    plt.locator_params(axis="y", nbins=10)
    n = "normalized" if normalize else ""
    ax.set_title("correlation of \u0394 G to NGS count")
    ax.set_xlabel(f"{n} \u0394 G")
    ax.set_ylabel("NGS count")
    ax.legend()

    n = "_normalized" if normalize else ""
    save_path = os.path.join(path, f"deltaG_NGScount{n}.png")
    plt.savefig(save_path)
    plt.close()


def plot_delta_G_observed_expected(df: pd.DataFrame,
                                   path: str,
                                   mode: str
                                   )-> None:
    '''
        Plots delta G against delta G of randomly shuffled or randomly cut
        sequences.
        :param df: data frame with the data
        :param path: path to the results folder
        :param mode: is either 'shuffled' or 'random'

        :return: None
    '''
    plt.rc("font", size=14)
    fig, axs = plt.subplots(2, 2, figsize=(7, 5), tight_layout=True)
    j = 0
    for i, strain in enumerate(["Cal07", "NC", "Perth", "BLEE"]):
        strain_df = df[df["strain"] == strain]
        if i == 2:
            j = 1
        for s in SEGMENTS:
            seg_df = strain_df[strain_df["segment"] == s]
            axs[i%2][j].scatter(seg_df["delta_G"], seg_df[f"delta_G_{mode}"], label=s, alpha=0.3)

        r, p = stats.pearsonr(strain_df["delta_G"], strain_df[f"delta_G_{mode}"])
        
        axs[i%2][j].set_title(f"{strain} (r = {r:.2})")
        axs[i%2][j].set_xlabel("\u0394 G")
        axs[i%2][j].set_ylabel(f"\u0394 G {mode}")
        axs[i%2][j].set_xlim([-650, 0])
        axs[i%2][j].set_xticks([-600, -300, 0])
        axs[i%2][j].set_ylim([-650, 0])
        axs[i%2][j].set_yticks([-600, -300, 0])
        axs[i%2][j].plot([0,1], [0,1], transform=axs[i%2][j].transAxes, c="grey", linestyle="--")
    
    axs[0][1].legend(bbox_to_anchor=(1.0,1.0))

    plt.locator_params(axis="y", nbins=10)

    save_path = os.path.join(path, f"deltaG_observed_{mode}.png")
    plt.savefig(save_path)
    plt.close()


def create_difference_boxplots(df: pd.DataFrame,
                               path: str,
                               mode: str
                               )-> None:
    '''
        Creates boxplots for the difference of the delta G against the expected
        values. Is done for the random and shuffled approach and saved in one
        figure. The boxes of the boxplot are split up by the segments.
        :param df: data frame with the values
        :param path: path to the results folder
        :param mode: is either 'shuffled' or 'random'
        
        :return: None
    '''
    # calculate differences of delta G to random cropped data
    df["delta_G_diff"] = df["delta_G"] - df[f"delta_G_{mode}"]

    plt.rc("font", size=14)
    fig, axs = plt.subplots(2, 2, figsize=(5, 5), tight_layout=True)
    j = 0
    for i, strain in enumerate(["Cal07", "NC", "Perth", "BLEE"]):
        strain_df = df[df["strain"] == strain]
        data = [strain_df.loc[strain_df["segment"] == s, "delta_G_diff"] for s in SEGMENTS]

        if i == 2:
            j = 1

        axs[i%2][j].boxplot(data, labels=SEGMENTS)
        axs[i%2][j].set_title(f"{mode} approach {strain}")
        axs[i%2][j].set_xlabel("segments")
        axs[i%2][j].set_ylabel(f"\u0394\u0394G (\u0394G - \u0394G {mode} approach)")
        axs[i%2][j].set_ylim([-60, 50])
        axs[i%2][j].axhline(y=0, c="r", linestyle="--", linewidth=0.5)

        for k, s in enumerate(SEGMENTS):
            if len(data[k]) > 1:
                statistic, p = stats.ttest_1samp(data[k], 0, alternative="less")
                symbol = get_stat_symbol(p)
                axs[i%2][j].annotate(symbol, (k+1, max(data[k])), horizontalalignment="center")

    save_path = os.path.join(path, f"boxplot_delta_delta_G_{mode}.png")
    fig.savefig(save_path)
    plt.close()


def count_bound_bases(df: pd.DataFrame,
                      sec_struct: str
                      )-> int:
    '''
        Checks if a nucleotide is bound or not and counts the overall bound
        bases over the whole data set.
        :param df: Data frame including the start and end point of the junction
                   sites
        :param sec_struct: secondary structure in Dot-Bracket notation

        :return: number of bound bases
    '''
    def is_bound_base(sec_struct: str, p: int)-> bool:
        if sec_struct[p] in ["(", ")"]:
            return True
        elif sec_struct[p] == ".":
            return False
        else:
            exit(f"Error: Unknown symbol ({sec_struct[p]})!")

    bound_bases = 0
    for r in df.iterrows():
        r = r[1]
        s = r["Start"]
        e = r["End"]
        if is_bound_base(sec_struct, s-1):
            bound_bases += 1
        if is_bound_base(sec_struct, e):
            bound_bases += 1

    return bound_bases


def check_secondary_structures(all_reads_dict: dict,
                               df: pd.DataFrame,
                               path: str
                               )-> None:
    '''
        Gets the junction site data and the secondary structures of the
        segments. Counts the number of bound nucleotides for each strain and
        segment and plots the results.
        :param all_reads_dict: Data including Start and End of junction site
        :param df: Secondary structure of the segments
        :param path: path to the results folder

        :return: None
    '''
    fig, axs = plt.subplots(4, 1, figsize=(20, 10), tight_layout=True)

    for i, (k, v) in enumerate(all_reads_dict.items()):
        for idx, seg in enumerate(SEGMENTS):
            sec_struct = df.loc[(df["segment"] == seg) & (df["strain"] == k)]["sec_structure"].tolist()[0]
            seg_df = v[v["Segment"] == seg]
            n = len(seg_df.index)
            if n == 0:
                obs_bound_ratio = 0.0
                exp_bound_ratio = 0.0
                symbol = ""
            else:
                obs_bound = count_bound_bases(seg_df, sec_struct)
                obs_bound_ratio = obs_bound/(n*2)

                s = (int(seg_df.Start.quantile(QUANT)), int(seg_df.Start.quantile(1-QUANT)))
                e = (int(seg_df.End.quantile(QUANT)), int(seg_df.End.quantile(1-QUANT)))

                sampling_df = generate_sampling_data(sec_struct, s, e, N_SAMPLES)
                exp_bound = count_bound_bases(sampling_df, sec_struct)
                exp_bound_ratio = exp_bound / (n_sampling*2)

                # n * 2 because Start and End is counted together
                result = stats.binomtest(obs_bound, n*2, exp_bound_ratio)
                symbol = get_stat_symbol(result.pvalue)

            axs[i].bar([f"{seg} obs", f"{seg} exp"], [obs_bound_ratio, exp_bound_ratio])
            axs[i].annotate(f"(n={n}) {symbol}", (idx*2+0.5,
                            max(obs_bound_ratio, exp_bound_ratio)),
                            horizontalalignment="center")

        axs[i].set_ylim(top=1.0)
        axs[i].set_xlabel("Segments")
        axs[i].set_ylabel("Ratio unbound/bound bases")
        axs[i].set_title(f"normalized bound bases of {k} (n={len(v.index)})", pad=10.0)
        axs[i].set_xticks(ticks=np.arange(0,16), labels=["obs", "exp"]*8)
    
    plt.legend(SEGMENTS)

    save_path = os.path.join(path, "bound_bases_occurrence.png")
    fig.savefig(save_path)
    plt.close()


if __name__ == "__main__":
    path = os.path.join(DATAPATH, "energy_calculation")
    path_full = os.path.join(path, "full_sequences")
    path_cropped = os.path.join(path, "cropped_sequences")
    results_path = os.path.join(RESULTSPATH, "free_energy_estimations")

    full_df = structure_dataframe(path_full)
    cropped_df = structure_dataframe(path_cropped)

    plot_deltaG_length(full_df, results_path, "full")
    plot_deltaG_length(cropped_df, results_path, "cropped")

    plot_deltaG_NGS(cropped_df, results_path, False)
    plot_deltaG_NGS(cropped_df, results_path, True)
    
    plot_delta_G_observed_expected(cropped_df, results_path, "random")
    plot_delta_G_observed_expected(cropped_df, results_path, "shuffled")
    
    create_difference_boxplots(cropped_df, results_path, "random")
    create_difference_boxplots(cropped_df, results_path, "shuffled")
    
    # Check junction sites and secondary structure
#    cleaned_data_dict = load_alnaji_excel()
 #   all_reads_dict = load_short_reads(cleaned_data_dict)
    
#    check_secondary_structures(all_reads_dict, full_df, results_path)

