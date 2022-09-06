'''
    This script calculates the number of duplicates between the two lineages of
    Alnaji et al. 2019. It is done for each of the strains.
'''
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

from matplotlib_venn import venn2

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import load_alnaji_excel, load_short_reads


def lineage_overlap(data: dict)-> None:
    '''
        Loops over the four strains and calculates the number of duplicate DI
        RNAs in the two lineages. Creates a venn diagramm as a result.
        :param data: dictionary with dataframes as values and strain name as
                     keys

        :return: None
    '''
    fig, axs = plt.subplots(4, 1, figsize=(8, 10), tight_layout=True)
    for i, name in enumerate(["Cal07", "NC", "Perth", "B_LEE"]):
        l1 = data[name + "_l1"]
        l2 = data[name + "_l2"]

        l1_df = set(l1[l1.columns[0:3]].apply(lambda x: "_".join(x.astype(str)), axis=1))
        l2_df = set(l2[l2.columns[0:3]].apply(lambda x: "_".join(x.astype(str)), axis=1))

        venn2([l1_df, l2_df], set_labels=("l1", "l2"), ax=axs[i])
        axs[i].set_title(name)

    fig.suptitle(f"overlap of lineage 1 and 2 for alnaji 2019 dataset")
    save_path = os.path.join(RESULTSPATH, "deletion_length_and_position", f"venn_diagramm_alnaji.png")
    plt.savefig(save_path)
    plt.close()


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    lineage_overlap(cleaned_data_dict)

