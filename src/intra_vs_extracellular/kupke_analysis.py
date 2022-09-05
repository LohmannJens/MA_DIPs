'''
    Analyzes the data set from Kupke et al. 2020. Main aspect is the difference
    of the pre and post group. Which indicates the difference between intra-
    and extracellular DI RNAs.
'''
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

from matplotlib_venn import venn2

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH
from utils import load_kupke


def venn_analysis(data: dict)-> None:
    '''
        Draws a venn diagramm for a given dataset with the groups 'pre' and
        'post'. Calculates the sizes of the groups and the duplicates that can
        be found in both.
        :param data: dict of the dataset with the strain name as key and a
                     data frame as value
        
        :return: None
    '''
    for k, v in data.items():
        fig, axs = plt.subplots(3, 1, figsize=(8, 10), tight_layout=True)
        for i, s in enumerate(["PB2", "PB1", "PA"]):
            v_s = v[v["Segment"] == s].copy()
            pre = v_s[v_s["infection"] == "pre"]
            post = v_s[v_s["infection"] == "post"]
            pre_post = pd.concat([pre["DI"], post["DI"]], ignore_index=True)

            pre_count = pre["DI"].nunique()
            post_count = post["DI"].nunique()
            pre_post_all = pre_post.nunique()
            duplicates = pre_count + post_count - pre_post_all
            pre_post_unique = pre_post_all - duplicates

            venn2(subsets=(pre_count-duplicates, post_count-duplicates, duplicates),
                  set_labels=("pre", "post"), ax=axs[i])
            axs[i].set_title(f"{s}")
        
        fig.suptitle(f"overlap of pre and post groups for {k}")
        
        save_path = os.path.join(RESULTSPATH, "intra_vs_extracellular", f"venn_diagramm_{k}.png")
        plt.savefig(save_path)
        plt.close()


if __name__ == "__main__":
    data_dict = load_kupke()
    venn_analysis(data_dict)
