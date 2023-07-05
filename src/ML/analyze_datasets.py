'''
    Analyzing the datasets that are used in the ML
    - checking the distribution
    - check for duplicates
    - dimensionality reduction
    - comparing start and end positions
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH
from ml_utils import load_all_sets, ngs_set_labels, generate_features, select_datasets


def check_distributions(df: pd.DataFrame)-> None:
    '''
        Checks the distribution of the NGS counts for all given datasets.
        Creates a histogram for that. A second one is using the log() function.
        :param df: data frame including all data sets

        :return: None
    '''
    dataset_names = df["dataset_name"].unique().tolist()
    dataset_names.append("all")
    for name in dataset_names:
        if name == "all":
            n_df = df
        else:
            n_df = df[df["dataset_name"] == name]

        fig, axs = plt.subplots(1, 2, figsize=(9, 4), tight_layout=True)

        axs[0].hist(n_df["NGS_read_count"], bins=100)
        axs[0].set_xlabel("NGS count")
        axs[0].set_ylabel("# occurrences")
        axs[0].set_title(f"NGS distribution for {name}")

        axs[1].hist(n_df["NGS_log_norm"], bins=100)
        axs[1].set_xlabel("log(NGS count)")
        axs[1].set_ylabel("# occurrences")
        axs[1].set_title(f"log(NGS) distribution for {name}")
    
        save_path = os.path.join(RESULTSPATH, "ML", f"distribution_{name}.png")
        plt.rc("font", size=14)
        plt.savefig(save_path)
        plt.close()


def check_stat_parameters(df: pd.DataFrame)-> None:
    '''
        Checks statistics for all data sets (like mean, max, min, etc). Creates
        a LaTeX table as a result and saves it.
        :param df: data frame including all data sets

        :return: None
    '''
    d = dict()
    g = df.groupby("dataset_name")
    d["Datasetname"] = g.size().index
    d["N_samples"] = g.size().tolist()
    d["Max"] = g.max("NGS_read_count")["NGS_read_count"].tolist()
    d["Min"] = g.min("NGS_read_count")["NGS_read_count"].tolist()
    d["Mean"] = g.mean("NGS_read_count")["NGS_read_count"].tolist()
    d["Median"] = g.median("NGS_read_count")["NGS_read_count"].tolist()
    r_df = pd.DataFrame(d)

    all_l = ["all", len(df.index),
             df["NGS_read_count"].max(), df["NGS_read_count"].min(),
             df["NGS_read_count"].mean(), df["NGS_read_count"].median()]
    r_df.loc[len(r_df.index)] = all_l
    r_df = r_df.set_index("Datasetname")
    
    path = os.path.join(RESULTSPATH, "ML", "stat_param_datasets.tex")
    r_df.to_latex(path, float_format="%.2f")


def check_duplicates(df: pd.DataFrame)-> None:
    '''
        For the Segment, Start, End combinations that occur more than once this
        function checks how many of them have the same label. Distinguish these
        groups more by checking if the matches are both from a PR8 data set or
        not.
        :param df: data frame including all data points

        :return: None
    '''
    n_bins = 3
    label_style = "median"

    df["DI"] = df["Segment"] + "_" + df["Start"].map(str) + "_" + df["End"].map(str)
    entries = df["DI"].tolist()
    dupl = [e for e in entries if entries.count(e) > 1]

    df["class"] = ngs_set_labels(df, n_bins, label_style)

    PR8_l = ["Pelz", "Kupke", "Alnaji2021"]

    all = 0
    equal = 0
    unequal = 0
    equal_PR8 = 0
    equal_not_PR8 = 0
    unequal_PR8 = 0
    unequal_not_PR8 = 0
    diff_list = list()
    scatter1_list = list()
    scatter2_list = list()
    for d in dupl:
        l_list = df.loc[df["DI"] == d]["class"].tolist()
        d_list = df.loc[df["DI"] == d]["dataset_name"].tolist()
        ngs_list =  df.loc[df["DI"] == d]["NGS_log_norm"].tolist()
        if len(l_list) == 2:
            diff_list.append(abs(ngs_list[0] - ngs_list[1]))
            scatter1_list.append(ngs_list[0])
            scatter2_list.append(ngs_list[1])
            all = all + 1
            if l_list[0] == l_list[1]:
                equal = equal + 1
                if d_list[0] in PR8_l and d_list[1] in PR8_l:
                    equal_PR8 = equal_PR8 + 1
                else:
                    equal_not_PR8 = equal_not_PR8 + 1
            else:
                unequal = unequal + 1
                if d_list[0] in PR8_l and d_list[1] in PR8_l:
                    unequal_PR8 = unequal_PR8 + 1
                else:
                    unequal_not_PR8 = unequal_not_PR8 + 1
                
    print(f"{all=}")
    print(f"++{unequal=}")
    print(f"++++{unequal_PR8=}")
    print(f"++++{unequal_not_PR8=}")
    print(f"++{equal=}")
    print(f"++++{equal_PR8=}")
    print(f"++++{equal_not_PR8=}")

    print(f"{equal/all}")

    print(np.mean(diff_list))

    fig, axs = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)
    axs.scatter(scatter1_list, scatter2_list)
    axs.plot([0, 1], [0, 1], c="r", linewidth=0.5, linestyle="--")
    if n_bins == 2:
        m = df["NGS_log_norm"].median()
        axs.hlines(m, 0, 1)
        axs.vlines(m, 0, 1)
    elif n_bins == 3:
        perc1 = df["NGS_log_norm"].quantile(q=0.33)
        axs.hlines(perc1, 0, 1)
        axs.vlines(perc1, 0, 1)
        perc2 = df["NGS_log_norm"].quantile(q=0.66)
        axs.hlines(perc2, 0, 1)
        axs.vlines(perc2, 0, 1)
        
    axs.set_xlabel("normalized NGS count of d1")
    axs.set_ylabel("normalized NGS count of d2")
    axs.set_title("comparing normalized NGS count for duplicates")

    save_path = os.path.join(RESULTSPATH, "ML", f"compare_duplicates_scatter.png")
    plt.rc("font", size=14)
    plt.savefig(save_path)
    plt.close()


def check_label_split(df: pd.DataFrame)-> None:
    '''
        Checks postion of the labeling thresholds for the different approaches
        and compares them to the NGS counts for all given datasets.
        :param df: data frame including all data sets

        :return: None
    '''
    dataset_names = df["dataset_name"].unique().tolist()
    dataset_names.append("all")
    for name in dataset_names:
        if name == "all":
            n_df = df
        else:
            n_df = df[df["dataset_name"] == name]

        fig, axs = plt.subplots(2, 2, figsize=(10, 10), tight_layout=True)
        for i in [0, 1]:
            if i == 0:
                n_bins = 2
                l = ["high", "low"]
            elif i == 1:
                n_bins = 3
                l = ["high", "mid", "low"]
            for j in [0, 1]:
                if j == 0:
                    label_style = "median"
                elif j == 1:
                    label_style = "pd.cut"

                axs[i,j].hist(n_df["NGS_log_norm"], bins=100)
                axs[i,j].set_xlabel("log(NGS count)")
                axs[i,j].set_ylabel("# occurrences")

                x_max = axs[i,j].get_ylim()[1]
                if n_bins == 2:
                    if j == 0:
                        m = n_df["NGS_log_norm"].median()
                    elif j == 1:
                        _, bins = pd.cut(n_df["NGS_log_norm"],
                                         bins=n_bins,
                                         labels=l,
                                         retbins=True,
                                         ordered=False)
                        m = bins[1]
                    axs[i,j].vlines(m, 0, x_max, color="red")
                elif n_bins == 3:
                    if j == 0:
                        perc1 = n_df["NGS_log_norm"].quantile(q=0.33)
                        perc2 = n_df["NGS_log_norm"].quantile(q=0.66)
                    elif j == 1:
                        _, bins = pd.cut(n_df["NGS_log_norm"],
                                         bins=n_bins,
                                         labels=l,
                                         retbins=True,
                                         ordered=False)
                        perc1 = bins[1]
                        perc2 = bins[2]

                    axs[i,j].vlines(perc1, 0, x_max, color="orange")
                    axs[i,j].vlines(perc2, 0, x_max, color="green")

                axs[i,j].set_title(f"{label_style} (bins={n_bins})")
    
        save_path = os.path.join(RESULTSPATH, "ML", f"check_labels_{name}.png")
        plt.rc("font", size=20)
        plt.savefig(save_path)
        plt.close()


def principal_component_analysis(X: pd.DataFrame,
                                 y: pd.Series,
                                 name: str
                                 )-> None:
    '''
        Runs PCA for two dimensions and plots the results.
        :param X: input features as data frame
        :param y: output vector as series
        :param name: string of all used datasets

        :return: None
    '''
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(X)
    pca_df = pd.DataFrame(data=principal_components, columns=["pc1", "pc2"])

    for l in y.unique():
        indices = y == l
        pc1 = pca_df.loc[indices, "pc1"]
        pc2 = pca_df.loc[indices, "pc2"]
        plt.scatter(pc1, pc2, alpha=0.2, label=l)

    plt.legend()
    plt.title("PCA")
    save_path = os.path.join(RESULTSPATH, "ML", f"pca_{name}.png")
    plt.savefig(save_path)
    plt.close()


def tsne(X: pd.DataFrame,
         y: pd.Series,
         name: str
         )-> None:
    '''
        Runs t-SNE for two dimensions and plots the results
        :param X: input features as data frame
        :param y: output vector as series
        :param name: string of all used datasets

        :return: None
    '''
    tsne = TSNE(n_components=2)
    X_embedded = tsne.fit_transform(X)
    tsne_df = pd.DataFrame(data=X_embedded, columns=["f1", "f2"])

    for l in y.unique():
        indices = y == l
        f1 = tsne_df.loc[indices, "f1"]
        f2 = tsne_df.loc[indices, "f2"]
        plt.scatter(f1, f2, alpha=0.2, label=l)

    plt.legend()
    plt.title("t-SNE")
    save_path = os.path.join(RESULTSPATH, "ML", f"tsne_{name}.png")
    plt.savefig(save_path)
    plt.close()


def run_dim_reduction(df: pd.DataFrame)-> None:
    '''
        runs two dimensionality reduction analyses to check resulting grouping
        induced by the feature values.
        :param df: data frame including all datasets

        :return: None
    '''
    n_bins = 2
    label_style = "median"
    y_column = "NGS_log_norm"

    features = ["Segment", "DI_Length", "Direct_repeat", "3_5_ratio", "length_proportion", "delta_G", "Peptide_Length"]
    df, feature_cols = generate_features(df, features, load_precalc=True)

    t_datasets = ["Alnaji2021"]
    v_datasets = list()

    X, y, X_val, y_val = select_datasets(df, t_datasets, v_datasets, feature_cols, n_bins, label_style, y_column)
    X = pd.concat([X, X_val])
    y = pd.concat([y, y_val])

    plt.rc("font", size=20)
    name = "_".join(t_datasets)
    principal_component_analysis(X, y, name)
    tsne(X, y, name)


def compare_start_vs_end(df: pd.DataFrame)-> None:
    '''
        Runs the plotting rountine for different slices of the dataset.
        :param df: data frame including all datasets

        :return: None
    '''
    slices = list(["all", "Alnaji2021", "Alnaji2021_cutoff", "Pelz_long"])
    
    for slc in slices:
        plot_start_vs_end(df, slc)


def plot_start_vs_end(df: pd.DataFrame, slc: str)-> None:
    '''
        Filters given dataframe by dataset names and plots start and end
        positions in a scatter.
        :param df: data frame including all datasets
        :param slc: Name of the datasets indicating how to slice the df

        :return: None
    '''
    n_bins = 2
    label_style = "median"
    y_column = "NGS_log_norm"

    if slc == "all":
        t_datasets = ["Alnaji2019_Cal07", "Alnaji2019_NC", "Alnaji2019_Perth", "Alnaji2019_BLEE", "Pelz", "Kupke"]
    elif slc == "Alnaji2021":
        t_datasets = ["Alnaji2021"]
    elif slc == "Alnaji2021_cutoff":
        t_datasets = ["Alnaji2021"]
        df = df[df["NGS_read_count"] > 5]
    elif slc == "Pelz_long":
        t_datasets = ["Pelz_long"]

    v_datasets = list()

    feature_cols = ["Start", "End", "Segment"]

    X, y, X_val, y_val = select_datasets(df, t_datasets, v_datasets, feature_cols, n_bins, label_style, y_column)
    X = pd.concat([X, X_val])
    y = pd.concat([y, y_val])
    X["y"] = y

    plt.rc("font", size=20)
    fig, ax = plt.subplots(1, 1, figsize=(7, 7), tight_layout=True)

    for l in y.unique():
        X_p = X[X["y"] == l]
        ax.scatter(X_p["Start"], X_p["End"], s=7, alpha=0.2, label=l)

    ax.set_xlabel("Start position")
    ax.set_ylabel("End position")
    ax.legend()
    ax.set_title(slc)

    save_path = os.path.join(RESULTSPATH, "ML", f"start_end_{slc}.png")
    plt.savefig(save_path)
    plt.close()


if __name__ == "__main__":
    all_df = load_all_sets()
    check_distributions(all_df)
    check_stat_parameters(all_df)
    check_duplicates(all_df)
    check_label_split(all_df)
    run_dim_reduction(all_df)
    compare_start_vs_end(all_df)

