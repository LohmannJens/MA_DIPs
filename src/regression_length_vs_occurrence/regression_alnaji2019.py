'''
Loads the Start and End points of the deletion sides from Alnaji 2019 and
does a linear and exponential regression.
'''

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from sklearn.linear_model import LinearRegression

sys.path.insert(0, "..")
from utils import SEGMENTS, get_seq_len, load_excel


if __name__ == "__main__":
    filepath = os.path.join("..", "..", "data", "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    cleaned_data_dict = load_excel(filepath)

    # plotting number of deletions per strain and segment against length of segment
    fig, axs = plt.subplots(8, 1, figsize=(5, 40), tight_layout=True)
    i = 0
    for key, value in cleaned_data_dict.items():
        x = list()
        y = list()
        for s in SEGMENTS:
            y.append(value.loc[value["Segment"] == s]["NGS_read_count"].sum())
            x.append(get_seq_len(key, s))
        axs[i].scatter(x, y, label=f"{key}")
        axs[i].set_title(f"{key}")
        axs[i].set_xlim(left=0)
        axs[i].set_xlabel("sequence position")

        x = np.array(x)
        y = np.array(y)
        x_exp = x
        y_exp = y
        idx = 0
        while (idx < len(y_exp)):
            if y_exp[idx] == 0:
         #       x_exp = np.delete(x_exp, idx)
          #      y_exp = np.delete(y_exp, idx)
                y_exp[idx] = 1
                idx += 1
            else:
                idx += 1
        
        x_plot = np.linspace(0, 2341, num=100).reshape((-1, 1))

        model = LinearRegression().fit(x.reshape((-1, 1)), y)
        axs[i].plot(x, model.predict(x.reshape((-1, 1))), label=f"linear model 1 (R²: {model.score(x.reshape((-1, 1)), y):.2f})")
        exp_model = np.polyfit(x_exp, np.log(y_exp), 1)
        axs[i].plot(x_plot, np.exp(exp_model[1]) * np.exp(exp_model[0] * x_plot), label="exponential model")
        axs[i].legend()

        for ix, s in enumerate(SEGMENTS):
            axs[i].annotate(s, (x[ix], y[ix]))
        i += 1

    save_path = os.path.join("results", "regression_alnaji2019.pdf")
    plt.savefig(save_path)

