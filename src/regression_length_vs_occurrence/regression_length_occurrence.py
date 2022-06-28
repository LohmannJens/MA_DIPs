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

from Bio import SeqIO
from sklearn.linear_model import LinearRegression

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, get_seq_len, load_excel, load_short_reads


def format_dataset_for_plotting(df, dataset_name: str)-> (list, list, list):
    '''
        formats the dataset to have it ready for plotting and doing the linear
        regression.
        :param df: data frame including the data to prepare
        :param dataset_name: indicates which data set is loaded and will be
                             formatted

        :return: tupel with three entries:
                    x values to plot
                    y values to plot
                    error values (if available)
    '''
    if dataset_name == "schwartz":
        x = df["Length"]
        y = df["average"]
        err = np.array(df["standard_deviation"]) / np.array(y).sum()
    else:
        x = list()
        y = list()
        err = np.empty(8)
        for s in SEGMENTS:
            y.append(df.loc[df["Segment"] == s]["NGS_read_count"].sum())
            x.append(get_seq_len(dataset_name, s))
    return np.array(x), np.array(y) / np.array(y).sum(), err


def fit_models_and_plot_data(x: list, y: list, y_exp: list, err: list, k: str)-> None:
    '''
        Creates the linear and exponential model for the given data and plots
        the results.
        :param x: data for x axis (segment length)
        :param y: data for y axis (DI occurrence) as relative values 
        :param y_exp: cleaned y values for the exponential model (no 0 values)
        :param err: data for the error bar (only for schwartz dataset)
        :param k: name of the strain/data set

        :return: None
    '''
    def label_scatter(x, y, k):
        for i, s in enumerate(SEGMENTS):
            if k == "exp" and s == "PB1":
                y[i] = y[i] - 0.007
            ax.annotate(s, (x[i], y[i]))
            if k == "IVA":
                ax.annotate(s, (x[i+8], y[i+8]))
                ax.annotate(s, (x[i+16], y[i+16]))

    x_plot = np.linspace(0, 2341, num=100).reshape((-1, 1))
    fig, ax = plt.subplots(1, 1, figsize=(10, 10), tight_layout=True)
 
    ax.scatter(x, y, label="observation")
    label_scatter(x, y, k)
    if k == "schwartz":
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
    label_scatter(x, y_expected, "exp")

    # excluding outliners to maximize R²
    x_exp = x
    if k == "Perth":
        x = x[:-1]
        y = y[:-1]
    if k in ["Cal07", "schwartz", "B_Lee"]:
        x = x[:-2]
        y = y[:-2]
    elif k == "IVA":
        x = np.delete(x, [5, 6, 7, 14, 15, 22, 23])
        y = np.delete(y, [5, 6, 7, 14, 15, 22, 23])

    # fit the models
    model = LinearRegression().fit(x.reshape((-1, 1)), y)
    exp_model = np.polyfit(x_exp, np.log(y_exp), 1)

    # predict values using the fitted models
    y_pred = model.predict(x.reshape((-1, 1)))
    y_exp_pred = np.exp(exp_model[1]) * np.exp(exp_model[0] * x_plot)

    # plotting the results
    score = model.score(x[:-2].reshape((-1, 1)), y[:-2])
    print(f"{k}: {score}")
    ax.plot(x, y_pred, label=f"linear model (R²: {score:.2f})", color="green")
    ax.plot(x_plot, y_exp_pred, label="exponential model", color="orange")

    # mark x axis intersection
    ax.plot((-model.intercept_/model.coef_), 0, 'ro')
    ax.annotate(f"{(-model.intercept_/model.coef_)[0]:.2f}", ((-model.intercept_/model.coef_), 0))

    # set labels and title
    ax.legend(loc="upper left")
    ax.set_title(f"{k}")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("sequence position")
    ax.set_ylabel("relative DI occurrence")

    # save final figure
    save_path = os.path.join(RESULTSPATH, "regression_length_count", f"{k}_regression_analysis.pdf")
    plt.savefig(save_path)
    plt.close()


def perform_regression_analysis(data: dict)-> None:
    '''
        Loops over all datasets. Creates a linear and an exponential model for
        all and plots these. Also creates one for all three Influenza strains
        together.
        :param data: dictionary with name as key and data frame as value

        :return: None
    '''
    for k, v in data.items():
        x, y, err = format_dataset_for_plotting(v, k)

        # clean data for exponential model (y == 0 is invalid)
        y_exp = y
        idx = 0
        while (idx < len(y_exp)):
            if y_exp[idx] == 0:
                y_exp[idx] = 0.000001
            idx += 1

        # collect data for the model over three IV A strains
        if k in ["Cal07", "NC", "Perth"]:
            if "x_IVA" in locals():
                x_IVA = np.concatenate((x_IVA, x))
                y_IVA = np.concatenate((y_IVA, y))
                y_IVA_exp = np.concatenate((y_IVA_exp, y_exp))
                y_IVA_expected = np.concatenate((y_IVA_expected, x / x.sum()))
            else:
                x_IVA = x
                y_IVA = y
                y_IVA_exp = y_exp
                y_IVA_expected = x / x.sum()

        fit_models_and_plot_data(x, y, y_exp, err, k)

    # do regression for all three IV A strains together
    fit_models_and_plot_data(x_IVA, y_IVA, y_IVA_exp, err, "IVA")


if __name__ == "__main__":
    filepath = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    cleaned_data_dict = load_excel(filepath)
    short_reads_filepath = os.path.join(DATAPATH, "alnaji2019", "Small_deletionSize_FA.xlsx")
    all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)

    xlsx_path = os.path.join(DATAPATH, "schwartz2016", "SchwartzLowen_Fig3a_MdckP3.xlsx")
    schwartz_data = pd.read_excel(xlsx_path, skiprows=[9,10,11,12])
    all_reads_dict["schwartz"] = schwartz_data
    
    perform_regression_analysis(all_reads_dict)

