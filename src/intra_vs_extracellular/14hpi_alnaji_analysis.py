'''
    Analyzes the data set from Alnaji et al. 2021. Main aspect is the
    difference of the intra and extra group.
'''
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

from matplotlib_venn import venn2

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH
from utils import load_hpi14_alnaji


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
 
        save_path = os.path.join(RESULTSPATH, "intra_vs_extracellular", "14hpi_venn_diagramm.png")
        plt.savefig(save_path)
        plt.close()


if __name__ == "__main__":
    data_dict = load_hpi14_alnaji()

    venn_analysis(data_dict)
