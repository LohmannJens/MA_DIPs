'''
Loads the Start and End points of the deletion sides from Alnaji 2019 and gives
insights about the data distribution.

1. Creates a histogram for each line in each strain containing the length of
   the deletion sides multiplied by their occurence.

2. Creates a plot where the start and end point of the deletion sides are
   plotted onto the sequence together with the NP density data.

3. Correlates occurrence of the position of the start and end points to the NP
   density.
'''

import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, load_excel, load_short_reads


def load_density_data(path: str)-> dict:
    '''
        Loads the density data for all eight segments.
        :param path: path to the location where the csv files are stored

        :return: dictionary with segment name as key and data frame as value
    '''
    density_dict = dict()
    for s in SEGMENTS:
        df = pd.read_csv(os.path.join(path, f"{s}.csv"), names=["x", "y"])
        density_dict[s] = df
    return density_dict


def load_WSN_data(dir: str)-> dict:
    '''
        Loads junction sites data for WSN strain from Mendes 2021.
        Formats it in the same way as the dataset from Alnaji 2019.
        :param dir: directory to the different .tsv files of the experiments

        :return: dictionary with one key (strain), value (data frame) pair
    '''
    dfs = list()
    for f in os.scandir(dir):
        df = pd.read_csv(f.path, sep="\t", na_values=["", "None"], keep_default_na=False)
        dfs.append(df)

    data = dict({"WSN": pd.concat(dfs)})
    return data


def plot_deletion_lengths(data: dict)-> None:
    '''
        creates a histogram for each strain, indicating the length of the
        deletions.
        :param data: dictionary with a data frame for each strain

        :return: None
    '''
    for key, value in data.items():
        # create a dict for each segment including the NGS read count
        count_dict = dict()
        for s in SEGMENTS:
            count_dict[s] = dict()

        for i, r in value.iterrows():
            if r["Length"] in count_dict[r["Segment"]]:
                count_dict[r["Segment"]][r["Length"]] += r["NGS_read_count"]
            else:
                count_dict[r["Segment"]][r["Length"]] = r["NGS_read_count"]

        # create a subplot for each key, value pair in count_dict
        fig, axs = plt.subplots(8, 1, figsize=(10, 20), tight_layout=True)
        fig.suptitle(f"absolute occurrences of deletions for {key}", x=0.3)
        for i, s in enumerate(SEGMENTS):
            axs[i].hist(count_dict[s].keys(), weights=count_dict[s].values(), bins=100)
            axs[i].set_title(f"{s}")
            axs[i].set_xlim(left=0)
            axs[i].set_xlabel("deletion length")

        save_path = os.path.join(RESULTSPATH, "deletion_length_and_position", f"{key}_length_del_hist.pdf")
        plt.savefig(save_path)
        plt.close()


def map_positions_to_density(data: dict, density_data: dict)-> dict:
    '''
        Maps the NP density to the start and end position of the deletion
        sites.
        :param data: dict with information about start and end position
        :param density_data: dict with density data (key is segment name)

        :return: counts positions found in NGS data
    '''
    NGS_dict = dict()
    for key, value in data.items():
        # create a dict for each segment using Start and End
        count_dict = dict()
        for s in SEGMENTS:
            count_dict[s] = dict()
        for i, r in value.iterrows():
            if r["Start"] in count_dict[r["Segment"]]:
                count_dict[r["Segment"]][r["Start"]] += r["NGS_read_count"]
            else:
                count_dict[r["Segment"]][r["Start"]] = r["NGS_read_count"]
            if r["End"] in count_dict[r["Segment"]]:
                count_dict[r["Segment"]][r["End"]] += r["NGS_read_count"]
            else:
                count_dict[r["Segment"]][r["End"]] = r["NGS_read_count"]

        # create a subplot for each key, value pair in count_dict
        fig, axs = plt.subplots(8, 1, figsize=(7, 14), tight_layout=True)
        fig.suptitle(f"position of deletions on sequence for {key}", x=0.3)
        for i, s in enumerate(SEGMENTS):
            l1 = axs[i].bar(count_dict[s].keys(), height=count_dict[s].values(), label="count")
            axs[i].set_ylabel("number of occurrences")
            l2, = axs[i].twinx().plot(density_data[s]["x"], density_data[s]["y"], label="NP density", alpha=0.5, color="red", fillstyle="full")
            axs[i].set_title(f"{s}")
            axs[i].set_xlim(left=0)
            axs[i].set_xlabel("sequence position")

        fig.legend([l1, l2], ["count", "NP density"])        
        save_path = os.path.join(RESULTSPATH, "deletion_length_and_position", f"{key}_del_position.pdf")
        plt.savefig(save_path)
        plt.close()

        NGS_dict[key] = count_dict

    return NGS_dict


def correlate_position_with_density(data: dict, density_data: dict)-> None:
    '''
        Does a correlation analysis for the NP density data and the NGS count
        of the DI RNA. normalizes the count of start and end points before
        doing it.
        :param data: NGS count data
        :param density_data: NP density data

        :return: None
    '''
    def get_expected_density(p, v):
        i = 0
        while v["x"][i] < p:
            i += 1
        return (v["y"][i] + v["y"][i-1]) / 2

    for k, v in data.items():
        fig, axs = plt.subplots(4, 2, figsize=(4, 8), tight_layout=True)
        j = 0
        for i, s in enumerate(SEGMENTS):
            # normalize NGS_read_count
            count_dict = v[s]
            if len(count_dict) != 0:
                x = np.array(list(count_dict.values())) / max(count_dict.values()) * 100
                y = list()
                for p in count_dict.keys():
                    try:
                        y.append(get_expected_density(p, density_data[s]))
                    except:
                        y.append(-1)
            axs[i%4][j].scatter(x, y, label=s)
            axs[i%4][j].set_title(s)
            axs[i%4][j].set_xlabel("normalized NGS count")
            axs[i%4][j].set_ylabel("NP density")
            axs[i%4][j].set_xlim(0, 100)
            axs[i%4][j].set_ylim(0, 100)

            if i == 3:
                j = 1

        fig.suptitle(k)

        savepath = os.path.join(RESULTSPATH, "deletion_length_and_position", f"{k}_count_density_correlation.pdf")
        fig.savefig(savepath)


if __name__ == "__main__":
    filepath = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    short_reads_filepath = os.path.join(DATAPATH, "alnaji2019", "Small_deletionSize_FA.xlsx")
    density_path = os.path.join(DATAPATH, "Lee2017", "csv_NPdensity")

    cleaned_data_dict = load_excel(filepath)
    all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)
    plot_deletion_lengths(all_reads_dict)

    # Plotting NP density against junction sites
    #    Cal07 data from Alnaji 2019
    Cal07_dens_path = os.path.join(density_path, "Cal07")
    Cal07_density_data = load_density_data(Cal07_dens_path)
    NGS_count_dict = map_positions_to_density(all_reads_dict, Cal07_density_data)
    correlate_position_with_density(NGS_count_dict, density_data)
    
    #    WSN data from Mendes 2021
    WSN_count_path = os.path.join(DATAPATH, "Mendes2021", "beta_sorting_BLASTresults")
    WSN_reads_dict = load_WSN_data(WSN_count_path)
    WSN_dens_path = os.path.join(density_path, "WSN")
    WSN_dens_data = load_density_data(WSN_dens_path)
    _ = map_positions_to_density(WSN_reads_dict, WSN_dens_data)

