'''
    Does a linear and exponential regression for data from Schwartz 2016 and 
    Alnaji 2019. Data is normalized by sum of y values for all data sets.
    Expected value is calculated by dividing length of each segment with sum of
    the length of all segements.

    Also creates a model for all three IV A strains together.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import get_seq_len, load_alnaji_excel, load_short_reads


def format_dataset_for_plotting(df: pd.DataFrame,
                                dataset_name: str,
                                del_indices: list=[]
                                )-> (list, list, list):
    '''
        formats the dataset to have it ready for plotting and doing the linear
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
        :param author:

        :return: None
    '''
    def label_scatter(x, y, k, segments):
        i = 0
        for s in segments:
            if k == "exp" and s == "PB1":
                y[i] = y[i] - 0.007
            ax.annotate(s, (x[i], y[i]))
            if k == "IVA":
                ax.annotate(s, (x[i+8], y[i+8]))
                ax.annotate(s, (x[i+16], y[i+16]))
            i += 1

    x_plot = np.linspace(0, 2341, num=100).reshape((-1, 1))
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), tight_layout=True)
    
    ax.scatter(x, y, label="observation")
    label_scatter(x, y, k, segments)
    ax.errorbar(x, y, yerr=err, fmt="o", capsize=5.0)
 
    # include expected values (perfect correlation DI count and length)
    if k == "IVA":
        y_expected = x.copy().astype("float64")
        y_expected[0:8] = y_expected[0:8] / y_expected[0:8].sum()
        y_expected[8:16] = y_expected[8:16] / y_expected[8:16].sum()
        y_expected[16:24] = y_expected[16:24] / y_expected[16:24].sum()
    else:
        y_expected = x / x.sum()

    ax.scatter(x, y_expected, label="expected", color="grey")
    label_scatter(x, y_expected, "exp", segments)

    # fit the models
    model = LinearRegression().fit(x.reshape((-1, 1)), y)
    exp_model = np.polyfit(x, np.log(y_exp), 1)

    # predict values using the fitted models
    y_pred = model.predict(x.reshape((-1, 1)))
    y_exp_pred = np.exp(exp_model[1]) * np.exp(exp_model[0] * x_plot)

    # plotting the results
    score = model.score(x[:-2].reshape((-1, 1)), y[:-2])
 #   formula = f"f(x) = x * {model.coef_} + {model.intercept_}"
    formula = ""
    ax.plot(x, y_pred, label=f"linear model {formula} (R²: {score:.2f})", color="green")
    ax.plot(x_plot, y_exp_pred, label="exponential model", color="orange")
    ax.plot(np.insert(x, 0, 0), np.insert(y_expected, 0, 0), color="grey")

    # mark x axis intersection
    ax.plot((-model.intercept_/model.coef_), 0, 'ro')
    ax.annotate(f"{(-model.intercept_/model.coef_)[0]:.2f}", ((-model.intercept_/model.coef_), 0))

    # set labels and title
    ax.legend(loc="upper left")
    ax.set_title(f"{k}")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0, top=0.65)
    ax.set_xlabel("sequence position")
    ax.set_ylabel("relative DI occurrence")

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

    # do regression for all three IV A strains together
    for k in ["Cal07", "NC", "Perth"]:
        v = data[k]
        x, y, err, segments = format_dataset_for_plotting(v, k)
        y_exp = clean_for_exp_analysis(y)
        if k == "Cal07":
            x_IVA = x
            y_IVA = y
            y_IVA_exp = y_exp
            y_IVA_expected = x / x.sum()
            IVA_err = err
        elif k in ["NC", "Perth"]:
            x_IVA = np.concatenate((x_IVA, x))
            y_IVA = np.concatenate((y_IVA, y))
            y_IVA_exp = np.concatenate((y_IVA_exp, y_exp))
            y_IVA_expected = np.concatenate((y_IVA_expected, x / x.sum()))
            IVA_err = np.concatenate((IVA_err, err))

    fit_models_and_plot_data(x_IVA, y_IVA, y_IVA_exp, IVA_err, "IVA", segments=segments)


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    
    xlsx_path = os.path.join(DATAPATH, "schwartz2016", "SchwartzLowen_Fig3a_MdckP3.xlsx")
    schwartz_data = pd.read_excel(xlsx_path, skiprows=[9,10,11,12])
    all_reads_dict["schwartz"] = schwartz_data
    
    perform_alnaji2019_regression_analysis(all_reads_dict)

