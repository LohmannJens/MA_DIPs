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


def venn_different_timepoints(data: dict)-> None:
    '''
        Draws a venn diagramm for a given dataset with the groups built by 
        different timepoints (3hpi, 6hpi, 24hpi)
        :param data: dict of the dataset with the strain name as key and a
                     data frame as value
        
        :return: None
    '''
    for k, v in data.items():
        fig, axs = plt.subplots(4, 1, figsize=(8, 10), tight_layout=True)
        for i, t in enumerate(["3hpi", "6hpi", "24hpi", "all"]):
            if t == "all":
                set1 = set(v[v["Timepoint"] == "3hpi"]["DI"])
                set2 = set(v[v["Timepoint"] == "6hpi"]["DI"])
                set3 = set(v[v["Timepoint"] == "24hpi"]["DI"])
                labels = ("3hpi", "6hpi", "24hpi")
            else:
                v_t = v[v["Timepoint"] == t].copy()
                set1 = set(v_t[v_t["Replicate"] == "Rep1"]["DI"])
                set2 = set(v_t[v_t["Replicate"] == "Rep2"]["DI"])
                set3 = set(v_t[v_t["Replicate"] == "Rep3"]["DI"])
                labels = ("Rep1", "Rep2", "Rep3")

            venn3([set1, set2, set3], set_labels=labels, ax=axs[i])
            axs[i].set_title(f"{t}")

        fig.suptitle(f"overlap of replicates at different timepoints for {k}")
        
        save_path = os.path.join(RESULTSPATH, "di_rna_conservation", f"venn_diagramm_{k}.png")
        plt.savefig(save_path)
        plt.close()


def slice_by_occurrence(data: dict, thresh: int, filepath: str)-> None:
    '''
        Allows to slice a dataset by the number of occurrences for each DI
        candidate. Counts the number of occurrences for each DI Name and then
        prints the candidates out, that are above a given threshhold.
        :param data: data frame containing the DI candidates
        :param thresh: occurrences larger or equal to this parameter are
                       included in the writen file
        :param filepath: indicating a path where to save the selected DI
                         candidates as .csv file
        :return: None
    '''
    df = data["PR8"]
    occur_df = pd.DataFrame(df.groupby(["DI"]).size())
    occur_df = occur_df.rename(columns={0: "Occurrences"})

    #print(occur_df.groupby(["occurrences"]).size())

    select_list = occur_df[occur_df["Occurrences"] >= thresh].index.values.tolist()
    selected_df = df[df["DI"].isin(select_list)].copy()
    selected_df = selected_df.groupby("DI").sum() # combine NGS_read_count
    selected_df.reset_index(inplace=True)

    write_df = selected_df[["DI", "NGS_read_count"]].copy()
    write_df[["Segment", "Start", "End"]] = write_df["DI"].str.split("_", expand=True)
    write_df.drop(["DI"], axis=1, inplace=True)
    write_df = write_df[["Segment", "Start", "End", "NGS_read_count"]] # reorder cols

    path = os.path.join(filepath, f"occurrence_{thresh}.csv")
    write_df.to_csv(path, index=False)


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
### compare different timepoints in alnaji 2021 ###
    data_dict = load_alnaji_2021()
    venn_different_timepoints(data_dict)

    # give them out as .csv and use DIP-DSA to investigate them
    # maybe also use other thressholds 6-9
#    slice_by_occurrence(data_dict, 9, ".")

    # also give out all DI candidates that are at least once in each timepoint:
#    get_intersection(data_dict["PR8"], "Timepoint", ["3hpi", "6hpi", "24hpi"])


### compare the three PR8 datasets with each other ###
    alnaji_df = data_dict["PR8"]
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

