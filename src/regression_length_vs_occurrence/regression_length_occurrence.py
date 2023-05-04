'''
    Does a linear and exponential regression for data from Schwartz 2016 and 
    Alnaji 2019. Data is normalized by sum of y values for all data sets.
    Expected value is calculated by dividing length of each segment with sum of
    the length of all segements.

    Also creates a model for all three IAV strains together.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, STRAINS
from utils import get_seq_len, load_alnaji_excel, load_short_reads


def format_dataset_for_plotting(df: pd.DataFrame,
                                dataset_name: str,
                                del_indices: list=[]
                                )-> (list, list, list):
    '''
        Formats the dataset to have it ready for plotting and doing the linear
        regression.
        :param df: data frame including the data to prepare
        :param dataset_name: indicates which data set is loaded and will be
                             formatted
        :param del_indices:

        :return: tupel with three entries:
                    x values to plot
                    y values to plot
                    error values (if available)
                    segments used (without excluded ones)
    '''
    if dataset_name == "schwartz":
        x = df["Length"]
        y = df["average"]
        err = np.array(df["standard_deviation"]) / np.array(y).sum()
        segments = SEGMENTS
    else:
        x = list()
        y = list()
        err = list()
        all_sum = df["NGS_read_count"].sum()
        segments = list()
        for i, s in enumerate(SEGMENTS):
            if i in del_indices:
                continue
            else:
                df_s = df.loc[df["Segment"] == s]
                y.append(df_s["NGS_read_count"].sum())
                x.append(get_seq_len(dataset_name, s))
                if df_s.size == 0:
                    err.append(0)
                else:
                    err.append(np.std(df_s["NGS_read_count"]) / all_sum)
                segments.append(s)

    return np.array(x), np.array(y) / np.array(y).sum(), err, segments


def fit_models_and_plot_data(x: list,
                             y: list,
                             y_exp: list,
                             err: list,
                             k: str,
                             author: str="",
                             segments: list=[]
                             )-> None:
    '''
        Creates the linear and exponential model for the given data and plots
        the results.
        :param x: data for x axis (segment length)
        :param y: data for y axis (DI occurrence) as relative values 
        :param y_exp: cleaned y values for the exponential model (no 0 values)
        :param err: data for the error bar (only for schwartz dataset)
        :param k: name of the strain/data set
        :param author: authors name, used to distinguish PR8 datasets
        :param segments: list of the segments to include

        :return: None
    '''
    def label_scatter(x, y, k, segments, a=""):
        if a != "":
            k = f"{a}_{k}"

        # This is for adjusting the labels by hand
        if k == "Cal07":
            y[4] = y[4] + 0.020 # NP
        elif k == "NC":
            y[3] = y[3] + 0.020 # HA
            y[4] = y[4] + 0.020 # NP
            x[4] = x[4] - 80 # NP
            x[5] = x[5] - 140 # NA
        elif k == "Perth":
            y[4] = y[4] + 0.020 # NP
        elif k == "BLEE":
            x[2] = x[2] - 50 # PA
            y[3] = y[3] + 0.005 # HA
            x[4] = x[4] - 130 # NP
        elif k == "three IAV strains":
            y[8] = y[8] - 0.015 # PB2
            x[8] = x[8] - 210 # PB2
            x[2] = x[2] - 130 # PA
            y[19] = y[19] - 0.015 # HA
            y[4] = y[4] - 0.015 # NP
            y[12] = y[12] - 0.015 # NP
            x[20] = x[20] + 30 # NP
            y[13] = y[13] - 0.015 # NA
            y[6] = y[6] - 0.015 # M
            y[14] = y[14] - 0.015 # M
            y[7] = y[7] - 0.015 # NS
            y[15] = y[15] - 0.015 # NS
            x[23] = x[23] - 70 # NS
        elif k == "Alnaji_PR8":
            x[2] = x[2] - 50 # PA
            y[4] = y[4] - 0.015 # NP
            y[6] = y[6] + 0.017 # M
        elif k == "Pelz_PR8":
            y[4] = y[4] + 0.020 # NP
            y[5] = y[5] + 0.020 # NA
        elif k == "Kupke_PR8":
            y[5] = y[5] + 0.023 # NA
            x[7] = x[7] - 50 # PA

        for i, s in enumerate(segments):
            ax.annotate(s, (x[i], y[i]))
            if k == "three IAV strains":
                ax.annotate(s, (x[i+8], y[i+8]))
                ax.annotate(s, (x[i+16], y[i+16]))

    x_plot = np.linspace(0, 2341, num=100).reshape((-1, 1))
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)
    
    ax.scatter(x, y, label="observed", marker=".")
    label_scatter(x.copy(), y.copy(), k, segments, author)
    ax.errorbar(x, y, yerr=err, fmt=".", capsize=5.0)
 
    # include expected values (perfect correlation DI count and length)
    if k == "three IAV strains":
        y_expected = x.copy().astype("float64")
        y_expected[0:8] = y_expected[0:8] / y_expected[0:8].sum()
        y_expected[8:16] = y_expected[8:16] / y_expected[8:16].sum()
        y_expected[16:24] = y_expected[16:24] / y_expected[16:24].sum()
    else:
        y_expected = x / x.sum()

    ax.scatter(x, y_expected, label="expected", color="grey", marker=".")

    # fit the models
    model = LinearRegression().fit(x.reshape((-1, 1)), y)
    exp_model = np.polyfit(x, np.log(y_exp), 1)

    # predict values using the fitted models
    y_pred = model.predict(x.reshape((-1, 1)))
    y_exp_pred = np.exp(exp_model[1]) * np.exp(exp_model[0] * x_plot)

    inter = model.intercept_
    coef = model.coef_[0]
    inter_p = -inter/coef

    # plotting the results
    score = model.score(x[:-2].reshape((-1, 1)), y[:-2])
    label = f"f(x) = {coef:.6f}*x {inter:.3f} (R²: {score:.2f})"
    ax.plot(np.insert(x, 0, inter_p), np.insert(y_pred, 0, 0), label=label, color="green")
    ax.plot(np.insert(x, 0, 0), np.insert(y_expected, 0, 0), color="grey")

    # mark x axis intersection
    ax.plot(inter_p, 0, 'ro')
    ax.annotate(f"{inter_p:.2f}", (inter_p, 0))

    # set labels and title
    ax.legend(loc="upper left")
    if k in STRAINS.keys():
        ax.set_title(f"{STRAINS[k]}")
    else:
        ax.set_title(k)
    ax.set_xlim(left=0, right=2600)
    ax.set_ylim(bottom=0, top=0.65)
    ax.set_xlabel("sequence length")
    ax.set_ylabel("relative DI RNA occurrence")

    # save final figure
    fname = f"{k}_regression_analysis.png"
    if author != "":
        fname = f"{author}_{fname}"
    save_path = os.path.join(RESULTSPATH, "regression_length_vs_occurrence", fname)
    plt.rc("font", size=12)
    plt.savefig(save_path)
    plt.close()


def clean_for_exp_analysis(y_exp: list)-> list:
    '''
        cleans the data to make an exponential regression possible. THis means
        to replace 0 with values very close to 0.
        :param y_exp: y values for the exponential regression

        :return: list with values equal to zero replaced by 0.000001
    '''
    i = 0
    while (i < len(y_exp)):
        if y_exp[i] == 0:
            y_exp[i] = 0.000001
        i += 1

    return y_exp


def linear_regression_analysis(strain: str,
                               df: pd.DataFrame,
                               del_indices: list=[],
                               author: str=""
                               )-> None:
    '''
        Runs the linear regression analysis
        :param strain: name of the influenza strain
        :param df: data frame with the necessary data
        :param del_indices: indices for datapoints to exclude in analysis
        :param author: authors name, used to distinguish PR8 datasets

        :return: None
    '''
    x, y, err, segments = format_dataset_for_plotting(df, strain, del_indices)
    y_exp = clean_for_exp_analysis(y)
    fit_models_and_plot_data(x, y, y_exp, err, strain, author, segments)


def perform_alnaji2019_regression_analysis(data: dict)-> None:
    '''
        Loops over all datasets. Creates a linear and an exponential model for
        all and plots these. Also creates one for all three Influenza strains
        together.
        :param data: dictionary with name as key and data frame as value

        :return: None
    '''
    for k, v in data.items():
        if k == "Perth":
            del_indices = [7]
        if k in ["Cal07", "schwartz", "B_Lee"]:
            del_indices = [6, 7]

        x, y, err, segments = format_dataset_for_plotting(v, k, del_indices)
        # clean data for exponential model (y == 0 is invalid)
        y_exp = clean_for_exp_analysis(y)
        fit_models_and_plot_data(x, y, y_exp, err, k, segments=segments)

    # do regression for all three IAV strains together
    for k in ["Cal07", "NC", "Perth"]:
        v = data[k]
        x, y, err, segments = format_dataset_for_plotting(v, k)
        y_exp = clean_for_exp_analysis(y)
        if k == "Cal07":
            x_IAV = x
            y_IAV = y
            y_IAV_exp = y_exp
            y_IAV_expected = x / x.sum()
            IAV_err = err
        elif k in ["NC", "Perth"]:
            x_IAV = np.concatenate((x_IAV, x))
            y_IAV = np.concatenate((y_IAV, y))
            y_IAV_exp = np.concatenate((y_IAV_exp, y_exp))
            y_IAV_expected = np.concatenate((y_IAV_expected, x / x.sum()))
            IAV_err = np.concatenate((IAV_err, err))

    fit_models_and_plot_data(x_IAV, y_IAV, y_IAV_exp, IAV_err, "three IAV strains", segments=segments)


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    
    xlsx_path = os.path.join(DATAPATH, "schwartz2016", "SchwartzLowen_Fig3a_MdckP3.xlsx")
    schwartz_data = pd.read_excel(xlsx_path, skiprows=[9,10,11,12])
    all_reads_dict["schwartz"] = schwartz_data
    
    perform_alnaji2019_regression_analysis(all_reads_dict)

