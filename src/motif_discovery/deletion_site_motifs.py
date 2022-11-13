'''
    Takes the fimo results from window xstreme run and does analysis with that.
    Before using this script the mode of the generation of the input sequences
    for MEME suite has to be defined (directly at beginning of main()).
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib_venn import venn2

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH
from plot_fimo_motifs_against_sequence import get_fimo_folders, load_fimo_files


W_SIZE = 50


def get_motif_overlap(row: object)-> str:
    '''
        Gets the row of a fimo.tsv file. Checks if the motif is fully on the
        start or end part of the sequence or if it is on both parts.
        :param row: row of a data frame, that was loaded from a fimo.tsv file

        :return: position of the motif in reference to the deletion site
                 'start', 'end' or 'overlap'
    '''
    if row["stop"] <= min(W_SIZE*2, row["del_start_pos"] + W_SIZE):
        return "start"
    elif row["start"] <= min(W_SIZE*2, row["del_start_pos"] + W_SIZE) and \
         row["stop"] > min(W_SIZE*2, row["del_start_pos"] + W_SIZE):
        return "overlap"
    elif row["start"] > min(W_SIZE*2, row["del_start_pos"] + W_SIZE):
        return "end"


def get_motif_overlap_remains(row: object)-> str:
    '''
        Gets the row of a fimo.tsv file. Checks if the motif is fully on the
        start or end part of the sequence or if it is on both parts. This
        function has to be used, if the input sequences where merged directly
        at the deletion site [only_remain == True] in preparation script.
        :param row: row of a data frame, that was loaded from a fimo.tsv file

        :return: position of the motif in reference to the deletion site
                 'start', 'end' or 'overlap'
    '''
    if row["stop"] <= min(W_SIZE, row["del_start_pos"]):
        return "start"
    elif row["start"] <= min(W_SIZE, row["del_start_pos"]) and \
         row["stop"] > min(W_SIZE, row["del_start_pos"]):
        return "split"
    elif row["start"] > min(W_SIZE, row["del_start_pos"]):
        return "end"


def calc_motif_start(row: object)-> int:
    '''
        Calculates the start position of the motif in reference to the full
        length sequence.
        :param row: row of a data frame, that was loaded from a fimo.tsv file

        :return: Start position of the motif on the full sequence
    '''
    if row["motif_overlap"] == "start":
        # if W_SIZE is bigger than deletion start position take 0 as start
        s = max(row["del_start_pos"] - W_SIZE, 0) + row["start"]
    elif row["motif_overlap"] == "end":
        s = row["del_end_pos"] - W_SIZE + row["start"]
    elif row["motif_overlap"] == "split":
        s = row["del_start_pos"] - W_SIZE + row["start"]
    return s


def calc_motif_end(row: object)-> int:
    '''
        Calculates the end position of the motif in reference to the full
        length sequence.
        :param row: row of a data frame, that was loaded from a fimo.tsv file

        :return: End position of the motif on the full sequence
    '''
    if row["motif_overlap"] == "start":
        # if W_SIZE is bigger than deletion start position take 0 as start
        e = max(row["del_start_pos"] - W_SIZE, 0) + row["stop"]
    elif row["motif_overlap"] == "end":
        e = row["del_end_pos"] - W_SIZE + row["stop"]
    elif row["motif_overlap"] == "split":
        e = row["del_end_pos"] - W_SIZE + row["stop"]
    return e


def check_strains_and_segments(df: object, motif: str)-> None:
    '''
        Groups data frame by strain and segment to see the distribution of a 
        given motif on these. Does this once including duplicates and once
        without. Prints the results to stdout.
        :param df: data frame including the fimo.tsv files
        :param motif: Name of a motif

        :return: None
    '''
    motif_df = df[df["motif_alt_id"] == motif]
    grouped = motif_df.groupby(["strain", "segment"])
    print(grouped.size().reset_index())
    # check for duplicate matches and clean them
    grouped_cleaned = motif_df.groupby(["motif_start", "motif_end", "strain", "segment"])
    grouped_cleaned_df = grouped_cleaned.size().reset_index()
    print(grouped_cleaned_df.groupby(["strain", "segment"]).size().reset_index())


def check_positional_distribution(df: object, motif: str)-> None:
    '''
        Groups data frame by motif position in reference to 'start', 'end' and
        'overlap'. Does this once including duplicates and once whtout. Prints
        the results to stdout.
        :param df: data frame including the fimo.tsv files
        :param motif: Name of a motif

        :return: None
    '''
    motif_df = df[df["motif_alt_id"] == motif]
    grouped = motif_df.groupby(["motif_overlap"])
    print(grouped.size().reset_index())
    # check for duplicate matches and clean them
    grouped_cleaned = motif_df.groupby(["motif_start", "motif_end", "motif_overlap"])
    grouped_cleaned_df = grouped_cleaned.size().reset_index()
    print(grouped_cleaned_df.groupby(["motif_overlap"]).size().reset_index())


def compare_DIPs_of_motifs(df, m1, m2)-> None:
    '''
        Gets the name of two motifs and creates a venn diagramm for the DI
        candidates that belong to those two. Shows the overlap the two given
        motifs have.
        :param df: data frame including the motifs and the names of the
                   DI candidates
        :param m1: Name of the first motif to compare
        :param m2: Name of the second motif to compare
    '''
    fig, axs = plt.subplots(1, 1, figsize=(5, 3), tight_layout=True)

    m1_set = set(df[df["motif_alt_id"] == m1]["sequence_name"])
    m2_set = set(df[df["motif_alt_id"] == m2]["sequence_name"])

    venn2([m1_set, m2_set], set_labels=(m1, m2))
    
    fig.suptitle(f"overlap of DI candidates for motif {m1} & {m2}")

    save_path = os.path.join(RESULTSPATH, "motif_discovery", f"DI_candidates_venn_{m1}_{m2}.png")
    plt.savefig(save_path)
    plt.close()


if __name__ == "__main__":
    # adjust parameters here
    combined = False
    only_remain = True
    motif = "MEME-3"

    # naming of source by parameters and loading of data
    folder = f"window_{W_SIZE}_sequences"
    if combined:
        folder = f"{folder}_combined"
    elif only_remain:
        folder = f"{folder}_onlyremain"
    path = os.path.join(DATAPATH, "meme_suite", "alnaji2019", folder ,"all_xstreme")
    fimo_folders = get_fimo_folders(path)
    fimo_df = load_fimo_files(fimo_folders)

    # drop unused columns
    fimo_df.drop(columns=fimo_df.columns[-1], axis=1, inplace=True)
    fimo_df.drop(columns=["motif_id", "strand"], axis=1, inplace=True)

    # full deletion site merged together
    if combined:
        # split info of sequence_name to columns to work with it
        fimo_df[["strain", "segment", "del_start_pos", "s", "del_end_pos", "e"]] = fimo_df["sequence_name"].str.split("_", expand=True)
        fimo_df.drop(columns=["s", "e"], axis=1, inplace=True)
        fimo_df.del_start_pos = fimo_df.del_start_pos.astype(int)
        fimo_df.del_end_pos = fimo_df.del_end_pos.astype(int)

        # get motif composition and position in reference to full sequence
        fimo_df["motif_overlap"] = fimo_df.apply(get_motif_overlap, axis=1)
        fimo_df["motif_start"] = fimo_df.apply(calc_motif_start, axis=1)
        fimo_df["motif_end"] = fimo_df.apply(calc_motif_end, axis=1)

        check_positional_distribution(fimo_df, motif)
        check_strains_and_segments(fimo_df, motif)

    # only part 1 and part 4 of deletion site is merged
    elif only_remain:
        fimo_df[["strain", "segment", "del_start_pos", "s", "del_end_pos", "e"]] = fimo_df["sequence_name"].str.split("_", expand=True)
        fimo_df.drop(columns=["s", "e"], axis=1, inplace=True)
        fimo_df.del_start_pos = fimo_df.del_start_pos.astype(int)
        fimo_df.del_end_pos = fimo_df.del_end_pos.astype(int)

        # get motif composition and position in reference to full sequence
        fimo_df["motif_overlap"] = fimo_df.apply(get_motif_overlap_remains, axis=1)
        fimo_df["motif_start"] = fimo_df.apply(calc_motif_start, axis=1)
        fimo_df["motif_end"] = fimo_df.apply(calc_motif_end, axis=1)

        check_positional_distribution(fimo_df, motif)
        check_strains_and_segments(fimo_df, motif)

        motif1 = "MEME-5"
        motif2 = "MEME-8"
        motif3 = "MEME-17"
        motif4 = "MEME-32"
        compare_DIPs_of_motifs(fimo_df, motif1, motif2)
        compare_DIPs_of_motifs(fimo_df, motif1, motif3)
        compare_DIPs_of_motifs(fimo_df, motif2, motif3)
        compare_DIPs_of_motifs(fimo_df, motif3, motif4)

    # start and end of sequence used as seperate input sequences
    else:
        fimo_df[["strain", "segment", "del_pos", "class"]] = fimo_df["sequence_name"].str.split("_", expand=True)
        fimo_df.del_pos = fimo_df.del_pos.astype(int)
        # get start and end of motif in correspondence to full sequence
        fimo_df["motif_start"] = fimo_df["del_pos"] - W_SIZE + fimo_df["start"]
        fimo_df["motif_end"] = fimo_df["del_pos"] - W_SIZE + fimo_df["stop"]

        grouped = fimo_df.groupby(["motif_alt_id", "strain", "segment", "motif_start", "matched_sequence"])
        grouped_df = grouped.size().reset_index()
        print(grouped_df)
