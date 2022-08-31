'''
    Influenza virus RNA polymerase simulation script that searches for t-loops
    and checks up/down stream for intermol bp.
    run_tloop_analysis() inspired by AJ te Velthuis, Sept 2020 from French
    et al. 2021
    used their main function and adapted it to automize it for alnajis dataset
'''
import os
import sys
import RNA

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, "..")
from utils import  DATAPATH, SEGMENTS
from utils import get_sequence
from French_sliding_window import run_tloop_analysis



if __name__ == "__main__":
    rna = "AGUAGAAACAAGGGUAUUUUUCUUUACUAGUCCGGUUGUUUUGGUUGCCACUAGUCUACCCUGCUUUUGCU"
    filename = "sanity_1_1.csv"
    run_tloop_analysis(rna, filename)

    path = os.path.join(DATAPATH, "energy_calculation", "sliding_window", filename)
    energy_df = pd.read_csv(path)

    plt.plot(energy_df["position"], energy_df["delta_G"])
    plt.show()
