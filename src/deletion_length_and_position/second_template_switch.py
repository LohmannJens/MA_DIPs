'''
    Loads the Start and End points of the deletion sides from Alnaji 2019 and
    gives insights about the data distribution.

    1.
    Creates a histogram for each line in each strain containing the length of
    the deletion sides multiplied by their occurence.
    
    2.
    Creates a plot where it shows the location of the start and the end points
    of the deletion site in reference to the full length segments

    3.
    Plots length of Start and End part of DI RNA as a scatter plot. Shows if
    they are equally distributed.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib.lines import Line2D

sys.path.insert(0, "..")
from utils import load_pelz_dataset, load_short_reads, load_alnaji_excel, load_WSN_data, load_full_alnaji2021, get_seq_len
from utils import SEGMENTS


def long_di_rna_comparision(d: dict):
    '''
    
    '''
    for k, v in d.items():
        df = v[v["NGS_read_count"] >= 30].copy()
        print(k)
        for s in ["PB2", "PB1", "PA"]:
            df_s = df[df["Segment"] == s].copy()
            short_DIs = df_s[(df_s["Start"] <1000) & (df_s["End"] > 1000)]

            long_DIs = df_s[~((df_s["Start"] <1000) & (df_s["End"] > 1000))].copy()
            long_DIs["deletion_len"] = long_DIs["End"] - long_DIs["Start"]



            short_DIs_n = len(short_DIs)
            long_DIs_n = len(df_s) - short_DIs_n

            print(f"{s}\t{long_DIs_n/len(df_s)}")
            print(long_DIs["deletion_len"].mean())
            print(long_DIs["deletion_len"].median())


def run_comparision_all(ds: list, dnames: list):
    '''
    
    '''
    for d, name in zip(ds, dnames):
        print(name)
        long_di_rna_comparision(d)


def create_start_end_connection_plot(df: pd.DataFrame,
                                     strain: str,
                                     segment: str,
                                     cutoff: int=1):
    '''
    
    '''
    df = df[(df["Segment"] == segment) & (df["NGS_read_count"] >= cutoff)].copy()
    max_val = get_seq_len(strain, segment)
    cm = plt.get_cmap("tab10")
    colors = [cm(1.*i/2) for i in range(2)]

    fig, ax = plt.subplots(figsize=(10, 5))

    ''' #plot with circles
    for i, row in df.iterrows():
        center = row["Start"] + (row["End"] - row["Start"]) / 2
        radius = (row["End"] - row["Start"]) / 2
        start_angle = 0
        end_angle = 180
        color = colors[0]
        if radius < 200:
            start_angle = 180
            end_angle = 0
            color = colors[1]
        half_cirlce = patches.Arc((center, 0), radius*2, radius*2, angle=0, theta1=start_angle, theta2=end_angle, color=color)
        ax.add_patch(half_cirlce)

    # add boxes for start and end of DI RNA sequence
    ax.add_patch(plt.Rectangle((0, -0.1), max_val, 0.1, alpha=0.7, color="black"))
    ax.add_patch(plt.Rectangle((0, 1), max_val, 0.1, alpha=0.7, color="black"))

    # change some values to improve figure
    ax.set_xlim(0, max_val)
    ax.set_ylim(-200, max_val / 2)
    ax.set_xticks(np.arange(0, max_val, 200))
    ax.set_yticks([])
    ax.set_xlabel("Nucleotide position")
    '''
    y_start = -max_val/2
    y_end = max_val/2

    for i, row in df.iterrows():
        if (row["Start"] < max_val/2) and (row["End"] > max_val/2):
            x_start = row["Start"]
            x_end = max_val - row["End"]
            line = Line2D([x_start, x_end], [y_start+10, y_end], linewidth=1)
            ax.add_line(line)
        
        else: # draw cicrle
            radius = (row["End"] - row["Start"]) / 2
            if row["Start"] < max_val/2: # is on bottom
                center = row["Start"] + radius
                start_angle = 180
                end_angle = 0
                y = y_start
            else: # is on top
                center = max_val - row["End"] + radius
                start_angle = 0
                end_angle = 180
                y = y_end

            half_cirlce = patches.Arc((center, y), radius*2, radius*2, angle=0, theta1=start_angle, theta2=end_angle)
            ax.add_patch(half_cirlce)
    
    # add boxes for start and end of DI RNA sequence
    ax.add_patch(plt.Rectangle((0, y_start), max_val/2, 10, alpha=0.7, color="black"))
    ax.add_patch(plt.Rectangle((max_val/2, y_start), 10, max_val+10, alpha=0.7, facecolor="white", edgecolor="black", hatch=r"//"))
    ax.add_patch(plt.Rectangle((0, y_end), max_val/2, 10, alpha=0.7, color="black"))

    ax.set_xlim(0, y_end*1.3)
    ax.set_ylim(y_start*1.3, y_end*1.3)



    # save figure
 #   plt.tight_layout()
  #  save_path = os.path.join(RESULTSPATH, "figure2", f"{strain}_{segment}_{cutoff}.png")
   # plt.savefig(save_path)
    plt.show()
    exit()

if __name__ == "__main__":
    plt.style.use('seaborn')
    ds = list()
    dnames = list()

    ds.append(load_pelz_dataset(long_dirna=True))
    dnames.append("Pelz")

    ds.append(load_short_reads(load_alnaji_excel()))
    dnames.append("Alnaji2019")

    ds.append(load_full_alnaji2021())
    dnames.append("Alnaji2021")

    ds.append(load_WSN_data("Mendes"))
    dnames.append("Mendes")

    ds.append(load_WSN_data("Boussier"))
    dnames.append("Boussier")

    #run_comparision_all(ds, dnames)

    for s in SEGMENTS:
        df = ds[0]["PR8"]
        create_start_end_connection_plot(df, "PR8", s, cutoff=1)

