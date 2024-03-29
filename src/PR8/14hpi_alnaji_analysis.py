'''
    Analyzes the data set from Alnaji et al. 2021. Main aspect is the
    difference of the intra and extra group at 14 hpi.
'''
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

from matplotlib_venn import venn2

sys.path.insert(0, "..")
sys.path.insert(0, "../di_rna_conservation")
sys.path.insert(0, "../regression_length_vs_occurrence")

from utils import DATAPATH, RESULTSPATH
from utils import load_hpi14_alnaji
from regression_length_occurrence import linear_regression_analysis


def venn_analysis(data: dict)-> None:
    '''
        Draws a venn diagramm for a given dataset with the groups 'internal'
        and 'external'.
        :param data: dict of the dataset with the strain name as key and a
                     data frame as value
        
        :return: None
    '''
    for k, v in data.items():
        fig, axs = plt.subplots(1, 1, figsize=(10, 10), tight_layout=True)
        internal = set(v[v["Class"] == "internal"]["DI"])
        external = set(v[v["Class"] == "external"]["DI"])
    
        venn2([internal, external], set_labels=("internal", "external"))

        fig.suptitle(f"overlap of internal and external DIs for {k}")
 
        save_path = os.path.join(RESULTSPATH, "PR8", "14hpi_venn_diagramm.png")
        plt.savefig(save_path)
        plt.close()


if __name__ == "__main__":
    data_dict = load_hpi14_alnaji()

    venn_analysis(data_dict)

    # Linear regression analysis
    strain = "PR8"
    df = data_dict["PR8"]
    inter_df = df.loc[df["Class"] == "internal"]
    extra_df = df.loc[df["Class"] == "external"]
    del_indices = [4]
    
    linear_regression_analysis(strain, df, del_indices=del_indices, author="Alnaji14hpi")
    linear_regression_analysis(strain, inter_df, del_indices=del_indices, author="Alnaji14hpi-intra")
    linear_regression_analysis(strain, extra_df, del_indices=del_indices, author="Alnaji14hpi-extra")

