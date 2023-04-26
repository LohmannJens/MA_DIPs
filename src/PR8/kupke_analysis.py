'''
    Analyzes the data set from Kupke et al. 2020. Main aspect is the difference
    of the pre and post group. Which indicates the difference between intra-
    and extracellular DI RNAs.
'''
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

from matplotlib_venn import venn3

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
        fig, axs = plt.subplots(4, 1, figsize=(8, 10), tight_layout=True)
        for i, s in enumerate(["PB2", "PB1", "PA", "all"]):
            if s == "all":
                v_s = v.copy()
            else:
                v_s = v[v["Segment"] == s].copy()
            pre = set(v_s[v_s["infection"] == "pre"]["DI"])
            post = set(v_s[v_s["infection"] == "post"]["DI"])
            post_sc = set(v_s[v_s["infection"] == "post_sc"]["DI"])

            venn3([pre, post, post_sc],
                  set_labels=("pre", "post", "post_sc"), ax=axs[i])
            axs[i].set_title(f"{s}")

        fig.suptitle(f"overlap of pre and post groups for {k}")
        
        save_path = os.path.join(RESULTSPATH, "PR8", f"venn_diagramms_kupke.png")
        plt.savefig(save_path)
        plt.close()


if __name__ == "__main__":
    data_dict = load_kupke(corrected=False)
    venn_analysis(data_dict)

