'''
    Loads the Start and End points of the deletion sides from Alnaji 2019 and
    compares their position to NP density data.

    1.
    Creates a plot where the start and end point of the deletion sides are
    plotted onto the sequence together with the NP density data.

    2.
    Correlates occurrence of the position of the start and end points to the NP
    density.
'''
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, QUANT, S_ROUNDS
from utils import load_alnaji_excel, load_short_reads, get_sequence, get_seq_len, get_stat_symbol, generate_sampling_data


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


def load_my_density_data(path: str)-> dict:
    '''
        Loads the density data for all eight segments.
        :param path: path to the location where the tsv files are stored

        :return: dictionary with segment name as key and data frame as value
    '''
    seg_mapper = dict({"PB2": "CY121687.1", "PB1": "CY121686.1",
                       "PA": "CY121685.1", "HA": "CY121680.1",
                       "M": "CY121681.1", "NA": "CY121682.1",
                       "NP":"CY121683.1", "NS": "CY121684.1"})

    density_dict = dict()
    filepath = os.path.join(path, "peaks_Sage2019.tsv")
    df = pd.read_csv(filepath, sep="\t", names=["Segment", "Start", "End", "Peak", "Height", "Sign"])
    for s in SEGMENTS:
        s_df = df.loc[df["Segment"] == seg_mapper[s]]

        start_df = pd.concat([s_df[["Start", "Height"]],
                              (s_df["Start"] - 1).to_frame()],
                             ignore_index=True)
        end_df = pd.concat([s_df[["End", "Height"]],
                            (s_df["End"] + 1).to_frame()],
                           ignore_index=True)
        start_df.rename(columns={"Start": "x", "Height": "y"}, inplace=True)
        end_df.rename(columns={"End": "x", "Height": "y"}, inplace=True)

        final_df = pd.concat([start_df, end_df], ignore_index=True, sort=True)
        extra_points_df = pd.DataFrame({"x": [0, get_seq_len("Cal07", s)], "y": [0, 0]})
        final_df = pd.concat([final_df, extra_points_df], ignore_index=True)
        final_df = final_df.fillna(0)
        final_df.sort_values(by=["x"], inplace=True)

        density_dict[s] = final_df
    return density_dict


def load_WSN_data(dir: str)-> dict:
    '''
        Loads junction sites data for WSN strain from Boussier 2020.
        Formats it in the same way as the dataset from Alnaji 2019.
        :param dir: directory to excel file with the data set

        :return: dictionary with one key (strain), value (data frame) pair
    '''
    df = pd.read_excel(dir, sheet_name=3, na_values=["", "None"], keep_default_na=False)
    df = df[df["Virus"] == "WT"].reset_index(drop=True)
    data = dict({"WSN": df})
    return data


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
        for i, s in enumerate(SEGMENTS):
            l1 = axs[i].twinx().bar(count_dict[s].keys(), height=count_dict[s].values(), label="count")
            l2, = axs[i].plot(density_data[s]["x"], density_data[s]["y"], label="NP density", alpha=0.5, color="green", fillstyle="full")
            axs[i].set_title(f"{s}")
            axs[i].set_xlim(left=0, right=max(density_data[s]["x"]))
#            axs[i].set_ylim(bottom=0, top=100)
            axs[i].set_xlabel("sequence position")
            axs[i].set_ylabel("high/low NP area")
 #           axs[i].axhline(y=5.0, color="red", linestyle="--")

        fig.suptitle(f"deletion position against NP areas for {key}", x=0.3)
        fig.legend([l1, l2], ["count", "NP density"])

        save_path = os.path.join(RESULTSPATH, "NP_density", f"{key}_del_position_NP_density.pdf")
        plt.savefig(save_path)
        plt.close()

        NGS_dict[key] = count_dict

    return NGS_dict


def compare_position_with_density(data: dict, density_data: dict, all_reads: dict)-> None:
    '''
        Checks how many of the junction sites are at a position of low NP
        density.
        :param data: NGS count data
        :param density_data: NP density data
        :param all_reads: dictionary with all start and end points

        :return: None
    '''
    def get_dens_at_pos(p, v):
        for _, (x, y) in v.iterrows():
            if x == p:
                return y
            elif x > p:
                return y
        return 0

    for k, v in data.items():
        fig, ax = plt.subplots(1, 1, figsize=(10, 5), tight_layout=True)

        for i, s in enumerate(SEGMENTS):
            count_dict = v[s]
            n = len(count_dict)
            if len(count_dict) == 0:
                continue

            # get expected values by sampling approach
            seq = get_sequence(k, s)
            seg_reads = all_reads[k].loc[all_reads[k]["Segment"] == s]
            start = (int(seg_reads.Start.quantile(QUANT)), int(seg_reads.Start.quantile(1-QUANT)))
            end = (int(seg_reads.End.quantile(QUANT)), int(seg_reads.End.quantile(1-QUANT)))
            n_sampling = n * S_ROUNDS
            sampling_data = generate_sampling_data(seq, start, end, n_sampling)
            samp_pos = sampling_data["Start"].tolist() + sampling_data["End"].tolist()

            # grouping by NP high/low
            y = [get_dens_at_pos(p, density_data[s]) for p in count_dict.keys()]
            y_exp = [get_dens_at_pos(p, density_data[s]) for p in samp_pos]
            below = y.count(0)
            below_exp = y_exp.count(0)
            obs_ratio = below/len(y)
            exp_ratio = below_exp/len(y_exp)

            # statistical testing
            result = stats.binomtest(below, n, exp_ratio)
            symbol = get_stat_symbol(result.pvalue)

            # plotting of the results
            ax.bar([f"{s} obs", f"{s} exp"], [obs_ratio, exp_ratio])
            ax.annotate(f"(n={n}) {symbol}", (i*2+0.5,
                        max(obs_ratio, exp_ratio)),
                        horizontalalignment="center")
            ax.set_xlabel("Segments")
            ax.set_ylabel("Ratio: sites below threshold/all sites")
            ax.set_xticks(ticks=np.arange(0,16), labels=["obs", "exp"]*8)

        plt.legend(SEGMENTS)
        fig.suptitle(k)

        savepath = os.path.join(RESULTSPATH, "NP_density", f"{k}_high_low_NP_areas.pdf")
        fig.savefig(savepath)
        plt.close()


if __name__ == "__main__":
    density_path = os.path.join(DATAPATH, "Lee2017", "csv_NPdensity")
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    del all_reads_dict["NC"]
    del all_reads_dict["Perth"]
    del all_reads_dict["B_LEE"]

    # Plotting NP data (from Sage 2019) against junction sites
    #    Cal07 data from Alnaji 2019
    Cal07_dens_path = os.path.join(density_path, "Cal07")
    Cal07_dens_data = load_my_density_data(Cal07_dens_path)
    NGS_count_dict = map_positions_to_density(all_reads_dict, Cal07_dens_data)
    compare_position_with_density(NGS_count_dict, Cal07_dens_data, all_reads_dict)
        
    #    WSN data from Mendes 2021
    WSN_count_file = os.path.join(DATAPATH, "Boussier2020", "Supplemental_Table_S2.xlsx")
    WSN_reads_dict = load_WSN_data(WSN_count_file)
    WSN_dens_path = os.path.join(density_path, "WSN")
    WSN_dens_data = load_my_density_data(WSN_dens_path)
    WSN_NGS_count_dict = map_positions_to_density(WSN_reads_dict, WSN_dens_data)
    compare_position_with_density(WSN_NGS_count_dict, WSN_dens_data, WSN_reads_dict)
    
