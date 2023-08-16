
import sys

import pandas as pd
import numpy as np

#plotting
import matplotlib.pyplot as plt
import seaborn as sns

from collections import Counter


sys.path.insert(0, "..")
from utils import get_sequence, load_alnaji_excel, load_short_reads, load_full_alnaji2021, load_pelz_dataset, load_kupke, generate_sampling_data
from utils import SEGMENTS, QUANT, N_SAMPLES

sys.path.insert(0, "../direct_repeats")
from search_direct_repeats import calculate_direct_repeat


def get_p_value_symbol(p):
    if p < 0.00001:
        return "*****"
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

def preprocess(strain, df):
    '''
    
    '''
    df["DI"] = df["Segment"] + "_" + df["Start"].map(str) + "_" + df["End"].map(str)
    return df

def plot_heatmap(y,x,vals,ax, format='.2f', cmap='coolwarm', vmin=0, vmax=1, cbar=False,cbar_ax=None, cbar_kws=None):
    '''Plot heatmap for values in vals, with x (columns) and y (rows) as labels.'''
    df = pd.DataFrame({'x':x,'y':y,'vals':vals})
    df = pd.pivot_table(df, index='x', columns='y', values='vals', sort=False)
    ax = sns.heatmap(df, fmt=format, annot=True, vmin=vmin, vmax=vmax, ax=ax, cbar=cbar, cmap=cmap, cbar_ax=cbar_ax, cbar_kws=cbar_kws)
    return ax


def generate_overlap_matrix_plot(dfs: list, dfnames: list):
    '''
    
    '''
    # initialize an empty matrix
    matrix_size = len(dfs)
    matrix = [[0] * matrix_size for _ in range(matrix_size)]

    # calculate the differences and populate the matrix
    for i in range(matrix_size):
        for j in range(matrix_size):
            set1 = set(dfs[i]["DI"])
            set2 = set(dfs[j]["DI"])

            matrix[i][j] = len(set1 & set2) / (max(len(set1), len(set2)))

    plt.imshow(matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    plt.xticks(np.arange(len(dfnames)), dfnames)
    plt.yticks(np.arange(len(dfnames)), dfnames)
    plt.xticks(rotation=45)
    plt.show()


def generate_maximum_overlap_candidates(dfs: list, dfnames: list):
    '''
    
    '''
    all_candidates = list()
    for df in dfs:
        all_candidates.extend(df["DI"].tolist())

    candidate_counts = Counter(all_candidates)

    counts = list()
    dir_reps = list()
    for cand, count in candidate_counts.items():
        if count > 10:
            seg, s, e = cand.split("_")
            seq = get_sequence("PR8", seg)
            dir_rep, _ = calculate_direct_repeat(seq, int(s), int(e), w_len=10, m=1)
            print(f"{cand}:\t{count}\t{dir_rep}")

            counts.append(count)
            dir_reps.append(dir_rep)

    d = dict({"counts": counts, "dir_reps": dir_reps})
    df = pd.DataFrame(d)
    
    grouped = df.groupby(["counts"])
    c = grouped.mean()
    
    print(c)

    plt.scatter(counts, dir_reps, alpha=0.05)
    plt.show()


if __name__ == "__main__":
    # load all data frames and preprocess with sequence_df(df)
    dfs = list()
    dfnames = list()

    '''
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    for k, v in all_reads_dict.items():
        dfs.append(preprocess(k, v))
        dfnames.append(k)
    
    df_alnaji = load_full_alnaji2021()
    df_alnaji = df_alnaji.groupby(["DI"]).sum(["NGS_read_count"]).reset_index()
    dfs.append(preprocess("PR8", df_alnaji))
    dfnames.append("Alnaji2021")
    
    '''
    alnaji_dict = load_full_alnaji2021()

    for k, v in alnaji_dict.items():
        for t in ["3hpi", "6hpi", "14hpi_internal", "14hpi_external", "24hpi"]:
            for r in ["Rep1", "Rep2", "Rep3"]:
                df = v[(v["Timepoint"] == t) & (v["Replicate"] == r)].copy()
                df = df.groupby(["DI", "Segment", "Start", "End"]).sum(["NGS_read_count"]).reset_index()
                dfs.append(preprocess(k, df))
                dfnames.append(f"Alnaji2021_{t}_{r}")   


    df_pelz = load_pelz_dataset(long_dirna=True)
    for k, v in df_pelz.items():
        dfs.append(preprocess(k, v))
        dfnames.append("Pelz_long")

    generate_overlap_matrix_plot(dfs, dfnames)
    generate_maximum_overlap_candidates(dfs, dfnames)