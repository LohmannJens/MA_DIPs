
import sys

import pandas as pd
import numpy as np

#plotting
import matplotlib.pyplot as plt
import seaborn as sns

from collections import Counter

sys.path.insert(0, "..")
from utils import get_sequence, load_full_alnaji2021, load_pelz_dataset, load_pelz_di244_dataset, load_pelz_PB2PB1PA_di244_dataset

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


def generate_maximum_overlap_candidates(dfs: list, thresh: int):
    '''
    
    '''
    all_candidates = list()
    for df in dfs:
        all_candidates.extend(df["DI"].tolist())

    candidate_counts = Counter(all_candidates)

    counts = list()
    dir_reps = list()
    for cand, count in candidate_counts.items():
        if count > thresh:
            seg, s, e = cand.split("_")
            seq = get_sequence("PR8", seg)
            dir_rep, _ = calculate_direct_repeat(seq, int(s), int(e), w_len=10, m=1)
            print(f"{cand}:\t{count}\t{dir_rep}")

            counts.append(count)
            dir_reps.append(dir_rep)

    '''
    d = dict({"counts": counts, "dir_reps": dir_reps})
    df = pd.DataFrame(d)
    
    grouped = df.groupby(["counts"])
    c = grouped.mean()
    
    print(c)

    plt.scatter(counts, dir_reps, alpha=0.05)
    plt.show()
    '''


if __name__ == "__main__":
    # load all data frames and preprocess with sequence_df(df)
    dfs = list()
    dfnames = list()

    '''
    alnaji_dict = load_full_alnaji2021()
    for k, v in alnaji_dict.items():
        for t in ["3hpi", "6hpi", "14hpi_internal", "14hpi_external", "24hpi"]:
            for r in ["Rep1", "Rep2", "Rep3"]:
                df = v[(v["Timepoint"] == t) & (v["Replicate"] == r)].copy()
                df = df.groupby(["DI", "Segment", "Start", "End"]).sum(["NGS_read_count"]).reset_index()
                dfs.append(preprocess(k, df))
                dfnames.append(f"Alnaji2021_{t}_{r}")   

    generate_overlap_matrix_plot(dfs, dfnames)
    generate_maximum_overlap_candidates(dfs, thresh=10)
    '''

    df_pelz = load_pelz_dataset(long_dirna=True)
    for k, v in df_pelz.items():
        dfs.append(preprocess(k, v))
        dfnames.append("Pelz_long")

    df_pelz = load_pelz_dataset(follow_up=True)
    for k, v in df_pelz.items():
        dfs.append(preprocess(k, v))
        dfnames.append("Pelz_followup")

    df_pelz = load_pelz_di244_dataset(by_time=False)
    for k, v in df_pelz.items():
        dfs.append(preprocess(k, v))
        dfnames.append("Pelz_Di244")
    
    df_pelz = load_pelz_PB2PB1PA_di244_dataset()
    for k, v in df_pelz.items():
        for exp in ["B", "C"]:
            for r in ["Rep1", "Rep2", "Rep3"]:
                df = v[(v["Experiment"] == exp) & (v["Replicate"] == r)].copy()
                dfs.append(preprocess(k, df))
                dfnames.append(f"Pelz_Di244_{exp}_{r}")

    generate_overlap_matrix_plot(dfs, dfnames)
    generate_maximum_overlap_candidates(dfs, thresh=6)

