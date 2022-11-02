'''
    Analyzes the datasets of PR8 from Alnaji, Pelz and Kupke. Compares if which
    DI candidates can be found in every of the datasets.
'''
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

from matplotlib_venn import venn3

sys.path.insert(0, "..")
sys.path.insert(0, "../relative_occurrence_nucleotides")
sys.path.insert(0, "../direct_repeats")
sys.path.insert(0, "../regression_length_vs_occurrence")

from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import load_alnaji_2021, load_pelz_dataset, load_kupke, create_sequence_library
from composition_junction_site import nucleotide_occurrence_analysis
from search_direct_repeats import direct_repeats_analysis
from regression_length_occurrence import linear_regression_analysis


def venn_different_datasets(df1: object, df2: object, df3: object, labels: list)-> None:
    '''
        
        :return: None
    '''
    set1 = set(df1["DI"])
    set2 = set(df2["DI"])
    set3 = set(df3["DI"])

    fig, axs = plt.subplots(1, 1, figsize=(10, 10), tight_layout=True)

    venn3([set1, set2, set3], set_labels=labels, ax=axs)
    axs.set_title(f"Venn diagramm")

    fig.suptitle(f"overlap of DI candidates for PR8 datasets from Alnaji, Pelz and Kupke")
        
    save_path = os.path.join(RESULTSPATH, "di_rna_conservation", f"venn_alnaji_pelz_kupke.png")
    plt.savefig(save_path)
    plt.close()


if __name__ == "__main__":
    alnaji_df = load_alnaji_2021()["PR8"]
    kupke_df = load_kupke(corrected=True)["PR8"]
    pelz_df = load_pelz_dataset()["PR8"]

    pelz_df["DI"] = pelz_df["Segment"] + "_" + pelz_df["Start"].astype(str) + "_" + pelz_df["End"].astype(str)

    # generate venn diagramm to compare the three datasets
    labels = ["Alnaji", "Kupke", "Pelz"]
    venn_different_datasets(alnaji_df, kupke_df, pelz_df, labels)
 
    PR8_dict = dict({"Alnaji": alnaji_df, "Kupke": kupke_df, "Pelz": pelz_df})
    for l in labels:
        df = PR8_dict[l]
        seq_dict = create_sequence_library({"PR8": df})
        # linear regression
        linear_regression_analysis("PR8", df, author=l)

        # nuc occurrence
        for s in SEGMENTS:
            nucleotide_occurrence_analysis(seq_dict, s, author=l)

        # direct repeats
        direct_repeats_analysis(seq_dict, 1, author=l)
        direct_repeats_analysis(seq_dict, 2, author=l)
        direct_repeats_analysis(seq_dict, 1, correction=True, author=l)

