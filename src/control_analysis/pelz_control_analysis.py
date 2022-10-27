'''
    Repeats the analyses on the dataset by Pelz et al. 2021. This is done to
    check if the findings are generalizeable.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats
from decimal import Decimal, ROUND_HALF_UP

sys.path.insert(0, "..")
sys.path.insert(0, "../relative_occurrence_nucleotides")
sys.path.insert(0, "../direct_repeats")
sys.path.insert(0, "../regression_length_vs_occurrence")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, QUANT, S_ROUNDS
from utils import load_pelz_dataset, get_stat_symbol, get_sequence, generate_sampling_data, create_sequence_library
from composition_junction_site import nucleotide_occurrence_analysis
from search_direct_repeats import count_direct_repeats_overall, include_correction
from regression_length_occurrence import format_dataset_for_plotting, fit_models_and_plot_data


def load_top_gain_de_novo()-> dict:
    '''
        Loads the top gain DI RNAs and top de novo DI RNAs
        
        :return: dictionary with strain as key and data frame as value
    '''
    path = os.path.join(DATAPATH, "Pelz2021", "Top_DI_RNA.xlsx")
    data = pd.read_excel(path)
    data_dict = dict({"PR8": data})
    return data_dict


def direct_repeats_analysis(seq_dict: dict, mode: int, top: bool=False, correction: bool=False, savepath: str="")-> None:
    '''
        Calculates the direct repeats of all sequences of the Pelz dataset.
        :param seq_dict: dictionary with the sequences
        :param mode: states which calculation mode is used. Check
                     composition_junction_site.py for more info.
        :param top: states if the whole dataset or just the top DI RNA are used
        :param correction: if True a correction calculation is made
        :param savepath: allows the user to give a path for saving

        :return: None
    '''
    fig, axs = plt.subplots(4, 2, figsize=(5, 10), tight_layout=True)
    for k, v in seq_dict.items():
        j = 0
        for i, s in enumerate(SEGMENTS):
            v_s = v.loc[(v["Segment"] == s)]
            seq = get_sequence(k, s)
            nuc_overlap_dict, _ = count_direct_repeats_overall(v_s, seq, mode)  
            n = len(v_s.index)
            if n <= 1:
                continue

            if correction:
                nuc_overlap_dict = include_correction(nuc_overlap_dict)

            start = (int(v_s.Start.quantile(QUANT)), int(v_s.Start.quantile(1-QUANT)))
            end = (int(v_s.End.quantile(QUANT)), int(v_s.End.quantile(1-QUANT)))
            sampling_data = generate_sampling_data(seq, start, end, n*S_ROUNDS)
            exp, _ = count_direct_repeats_overall(sampling_data, seq, mode)

            x = list(nuc_overlap_dict.keys())
            h = np.array(list(nuc_overlap_dict.values()))
            h_exp = np.array(list(exp.values()))

            # test statistical significance
            f_obs = list()
            f_exp = list()
            for a in x:
                if correction:
                    f_obs.extend([a]*int(Decimal(nuc_overlap_dict[a]).to_integral_value(rounding=ROUND_HALF_UP)))
                else:
                    f_obs.extend([a]*nuc_overlap_dict[a])
                f_exp.extend([a]*exp[a])
            f_obs = np.array(f_obs)
            f_exp = np.array(f_exp)

            res = stats.mannwhitneyu(f_obs, f_exp)
            symbol = get_stat_symbol(res.pvalue)

            axs[i%4, j].bar(x=x, height=h/h.sum(), width=-0.4, align="edge", label="observed")
            axs[i%4, j].bar(x=x, height=h_exp/h_exp.sum(), width=0.4, align="edge", label="expected")
            axs[i%4, j].set_xlabel("number of overlapping nucleotides")
            axs[i%4, j].set_ylabel("relative occurrence")
            axs[i%4, j].set_title(f"{s} (n={n}) {symbol}")
            axs[i%4, j].legend(loc="upper right")
            axs[i%4, j].set_ylim(bottom=0.0, top=1.0)

            if i == 3:
                j = 1
    if savepath == "":
        if top:
            fname = f"direct_repeats_{k}_top_mode{mode}.pdf"
        else:
            fname = f"direct_repeats_{k}_mode{mode}.pdf"
        if correction:
            fname = f"direct_repeats_{k}_mode{mode}_corr.pdf"
        savepath = os.path.join(RESULTSPATH, "control_analysis", fname)
    plt.savefig(savepath)
    plt.close()

def linear_regression_analysis(strain: str, df: object, dst: str) -> None:
    '''
        Runs the linear regression analysis
        :param strain: name of the influenza strain
        :param df: data frame with the necessary data
        :param dst: path to the save destination

        :return: None
    '''
    x, y, err = format_dataset_for_plotting(df, strain)
    y_exp = y
    fit_models_and_plot_data(x, y, y_exp, err, strain)
    src = os.path.join(RESULTSPATH, "regression_length_vs_occurrence", f"{strain}_regression_analysis.png")
    os.rename(src, dst)


if __name__ == "__main__":
    # analysis for top gain and top de novo DI RNAs
    
    top_di_rna = load_top_gain_de_novo()  
    top_seq_dict = create_sequence_library(top_di_rna)
    direct_repeats_analysis(top_seq_dict, 1, True)
    direct_repeats_analysis(top_seq_dict, 2, True)
    
    # analysis for the whole Pelz dataset
    data_dict = load_pelz_dataset()
    seq_list_dict = create_sequence_library(data_dict)
    for s in SEGMENTS:
        nucleotide_occurrence_analysis(seq_list_dict, s)
        src = os.path.join(RESULTSPATH, "relative_occurrence_nucleotides", f"PR_{s}.png")
        dst = os.path.join(RESULTSPATH, "control_analysis", f"PR_{s}_nucleotide_occurrence.png")
        if os.path.exists(src):
            os.rename(src, dst)

    direct_repeats_analysis(seq_list_dict, 1)
    direct_repeats_analysis(seq_list_dict, 2)
    direct_repeats_analysis(seq_list_dict, 1, correction=True)
    
    # do the linear regression analysis (length vs NGS count)
    strain = "PR8"
    df = data_dict[strain]
    dst = os.path.join(RESULTSPATH, "control_analysis", f"PR8_regression_analysis.png")
    linear_regression_analysis(strain, df, dst)
