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
from utils import DATAPATH, RESULTSPATH
from utils import load_alnaji_2021, load_pelz_dataset, load_kupke


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


def get_intersection(df: object, col_name: str, choices: list)-> dict:
    '''
        Gets only the rows of the data frame that are observed in every set of
        the three choices in a given column.
        :param df: data frame of the DI candidates
        :param col_name: name of the column that splits the data
        :param choices: list with three entries, indicating the parameter that
                        defines the three sets

        :return: data frame only containing the selected rows
    '''
    set1 = set(df[df[col_name] == choices[0]]["DI"])
    set2 = set(df[df[col_name] == choices[1]]["DI"])
    set3 = set(df[df[col_name] == choices[2]]["DI"])

    inter_set = set1.intersection(set2, set3)
    inter_df = df[df["DI"].isin(inter_set)]
    
    return inter_df


if __name__ == "__main__":
    alnaji_df = load_alnaji_2021()["PR8"]
    kupke_df = load_kupke(corrected=True)["PR8"]
    pelz_df = load_pelz_dataset()["PR"]

    pelz_df["DI"] = pelz_df["Segment"] + "_" + pelz_df["Start"].astype(str) + "_" + pelz_df["End"].astype(str)

    # generate venn diagramm to compare the three datasets
    labels = ["Alnaji", "Kupke", "Pelz"]
    venn_different_datasets(alnaji_df, kupke_df, pelz_df, labels)
    
    # set label for each dataset and merge them together
    for df, l in zip([alnaji_df, kupke_df, pelz_df], labels):
        df["dataset"] = l

    full_df = pd.concat([alnaji_df, kupke_df, pelz_df], join="inner", ignore_index=True)

    #for kupke and alnaji merge all categories/timepoints
    #get intersecting DI candidates and do analysis with them
    inter_df = get_intersection(full_df, "dataset", labels)
    print(inter_df)

