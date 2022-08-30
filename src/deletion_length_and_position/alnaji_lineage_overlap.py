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
    for name in ["Cal07", "NC", "Perth", "B_LEE"]:
        l1 = data[name + "_l1"]
        l2 = data[name + "_l2"]

        l1_df = l1[l1.columns[0:3]].apply(lambda x: "_".join(x.astype(str)), axis=1)
        l2_df = l2[l2.columns[0:3]].apply(lambda x: "_".join(x.astype(str)), axis=1)
        l1_l2_df = pd.concat([l1_df, l2_df], ignore_index=True)

        n_l1 = l1_df.nunique()
        n_l2 = l2_df.nunique()
        n_l1l2 = n_l1 + n_l2
        l1l2_all = l1_l2_df.nunique() # containing the duplicates once
        l1l2_duplicates = n_l1l2 - l1l2_all
        l1l2_unique = l1l2_all - l1l2_duplicates # only containing unique values

        venn2(subsets=(n_l1, n_l2, l1l2_duplicates),
              set_labels=("l1", "l2"))
        plt.title(f"{name} (quantity of duplicates: {l1l2_duplicates/l1l2_all*100:.3} %)")
        plt.show()



if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    lineage_overlap(cleaned_data_dict)

