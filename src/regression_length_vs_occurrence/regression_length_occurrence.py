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
        err = np.array(df["standard_deviation"])
    else:
        x = list()
        y = list()
        err = np.empty(8)
        for s in SEGMENTS:
            y.append(df.loc[df["Segment"] == s]["NGS_read_count"].sum())
            x.append(get_seq_len(dataset_name, s))
    return np.array(x), np.array(y) / np.array(y).sum(), err / np.array(y).sum()


def perform_regression_analysis(data: dict)-> None:
    '''
        loops over all datasets. Creates a linear and an exponential model for
        all and plots these. Also creates one for all three Influenza strains
        together.
        :param data: dictionary with name as key and data frame as value

        :return: None
    '''
    for k, v in data.items():
        fig, ax = plt.subplots(1, 1, figsize=(10, 10), tight_layout=True)
        x, y, err = format_dataset_for_plotting(v, k)

        ax.scatter(x, y, label="observation")
        if err.size != 0:
            ax.errorbar(x, y, yerr=err, fmt="o", capsize=5.0)
        for i, s in enumerate(SEGMENTS):
            ax.annotate(s, (x[i], y[i]))

        # clean data for exponential model (y == 0 is invalid)
        y_exp = y
        idx = 0
        while (idx < len(y_exp)):
            if y_exp[idx] == 0:
                y_exp[idx] = 0.000001
            idx += 1
        
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

        x_plot = np.linspace(0, 2341, num=100).reshape((-1, 1))

        model = LinearRegression().fit(x.reshape((-1, 1)), y)
        exp_model = np.polyfit(x, np.log(y_exp), 1)

        y_pred = model.predict(x.reshape((-1, 1)))
        y_exp_pred = np.exp(exp_model[1]) * np.exp(exp_model[0] * x_plot)

        ax.plot(x, y_pred, label=f"linear model (R²: {model.score(x.reshape((-1, 1)), y):.2f})")
        ax.plot(x_plot, y_exp_pred, label="exponential model")

        ax.scatter(x, x / x.sum(), label="expected", color="grey")

        ax.legend(loc="upper left")
        ax.set_title(f"{k}")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("sequence position")
        ax.plot((-model.intercept_/model.coef_), 0, 'ro')
        ax.annotate(f"{(-model.intercept_/model.coef_)[0]:.2f}", ((-model.intercept_/model.coef_), 0))

        save_path = os.path.join(RESULTSPATH, "regression_length_count", f"{k}_regression_analysis.pdf")
        plt.savefig(save_path)
        plt.close()

    # do regression for all three IV A strains together
    fig, ax = plt.subplots(1, 1, figsize=(10, 10), tight_layout=True)
    ax.scatter(x_IVA, y_IVA, label="observation")

    model = LinearRegression().fit(x_IVA.reshape((-1, 1)), y_IVA)
    exp_model = np.polyfit(x_IVA, np.log(y_IVA_exp), 1)

    y_pred = model.predict(x_IVA.reshape((-1, 1)))
    y_exp_pred = np.exp(exp_model[1]) * np.exp(exp_model[0] * x_plot)

    ax.plot(x_IVA, y_pred, label=f"linear model (R²: {model.score(x_IVA.reshape((-1, 1)), y_IVA):.2f})")
    ax.plot(x_plot,  y_exp_pred, label="exponential model")

    ax.scatter(x_IVA, y_IVA_expected, label="expected", color="grey")

    for i, s in enumerate(SEGMENTS):
        ax.annotate(s, (x_IVA[i], y_IVA[i]))
        ax.annotate(s, (x_IVA[i+8], y_IVA[i+8]))
        ax.annotate(s, (x_IVA[i+16], y_IVA[i+16]))
    ax.legend(loc="upper left")
    ax.set_title(f"All three IV A strains")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("sequence position")
    ax.plot((-model.intercept_/model.coef_), 0, 'ro')
    ax.annotate(f"{(-model.intercept_/model.coef_)[0]:.2f}", ((-model.intercept_/model.coef_), 0))

    save_path = os.path.join(RESULTSPATH, "regression_length_count", f"IVA_regression_analysis.pdf")
    plt.savefig(save_path)
    plt.close()


if __name__ == "__main__":
    filepath = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    cleaned_data_dict = load_excel(filepath)
    short_reads_filepath = os.path.join(DATAPATH, "alnaji2019", "Small_deletionSize_FA.xlsx")
    all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)

    xlsx_path = os.path.join(DATAPATH, "schwartz2016", "SchwartzLowen_Fig3a_MdckP3.xlsx")
    schwartz_data = pd.read_excel(xlsx_path, skiprows=[9,10,11,12])
    all_reads_dict["schwartz"] = schwartz_data
    
    perform_regression_analysis(all_reads_dict)

