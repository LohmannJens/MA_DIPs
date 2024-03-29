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
from matplotlib.patches import Rectangle

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, STRAINS, QUANT, N_SAMPLES
from utils import load_alnaji_excel, load_short_reads, get_sequence, get_seq_len, get_stat_symbol, generate_sampling_data, load_full_alnaji2021, load_pelz_dataset, load_WSN_data


def load_Lee_density_data(path: str)-> dict:
    '''
        Loads the density data for all eight segments. Uses the Lee et al. 2017
        paper as source.
        :param path: path to the location where the csv files are stored

        :return: dictionary with segment name as key and data frame as value
    '''
    density_dict = dict()
    for s in SEGMENTS:
        df = pd.read_csv(os.path.join(path, f"{s}.csv"), names=["x", "y"])
        density_dict[s] = df
    return density_dict


def load_Sage_density_data(path: str)-> dict:
    '''
        Loads the density data for all eight segments. Uses the Le Sage et al.
        2019 paper as source.
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


def load_Williams_density_data(path: str, mode: str)-> dict:
    '''
        Loads the density data for all eight segments. Uses the Williams et al.
        2018 paper as source.
        :param path: path to the location where the tsv files are stored

        :return: dictionary with segment name as key and data frame as value
    '''
    density_dict = dict()
    filepath = os.path.join(path, "Williams_2018.csv")
    df = pd.read_csv(filepath,
                     na_values=["", "None"],
                     keep_default_na=False,
                     converters={"Segment": str, "Start": int,"End": int, "Peak": int})

    if mode == "high":
        n = 100
    elif mode == "low":
        n = -100

    df = df[df["Peak"] == n]

    for s in SEGMENTS:
        s_df = df.loc[df["Segment"] == s]

        start_df = pd.concat([s_df[["Start", "Peak"]],
                              (s_df["Start"] - 1).to_frame()],
                             ignore_index=True)
        end_df = pd.concat([s_df[["End", "Peak"]],
                            (s_df["End"] + 1).to_frame()],
                           ignore_index=True)
        start_df.rename(columns={"Start": "x", "Peak": "y"}, inplace=True)
        end_df.rename(columns={"End": "x", "Peak": "y"}, inplace=True)

        final_df = pd.concat([start_df, end_df], ignore_index=True, sort=True)

        extra_points_df = pd.DataFrame({"x": [0, get_seq_len("PR8", s)], "y": [0, 0]})
        final_df = pd.concat([final_df, extra_points_df], ignore_index=True)
        final_df = final_df.fillna(0)
        final_df.sort_values(by=["x"], inplace=True)

        density_dict[s] = final_df
    return density_dict


def map_positions_to_density(data: dict,
                             density_data: dict,
                             density_data_2: dict=dict(),
                             author: str=""
                             )-> dict:
    '''
        Maps the NP density to the start and end position of the deletion
        sites.
        :param data: dict with information about start and end position
        :param density_data: dict with density data (key is segment name)

        :return: counts positions found in NGS data
    '''
    plt.rc("font", size=14)
    NGS_dict = dict()
    for k, v in data.items():
        # create a dict for each segment using Start and End
        count_dict = dict()
        for s in SEGMENTS:
            count_dict[s] = dict()
        for i, r in v.iterrows():
            if r["Start"] in count_dict[r["Segment"]]:
                count_dict[r["Segment"]][r["Start"]] += r["NGS_read_count"]
            else:
                count_dict[r["Segment"]][r["Start"]] = r["NGS_read_count"]
            if r["End"] in count_dict[r["Segment"]]:
                count_dict[r["Segment"]][r["End"]] += r["NGS_read_count"]
            else:
                count_dict[r["Segment"]][r["End"]] = r["NGS_read_count"]
        
        # create a subplot for each key, value pair in count_dict
        fig, axs = plt.subplots(8, 1, figsize=(10, 14), tight_layout=True)
        for i, s in enumerate(SEGMENTS):
            counts = count_dict[s].values()
            axs[i].bar(count_dict[s].keys(), height=counts, label="count")            
            l = density_data[s]["x"].tolist()
            for j in range(2, len(l)-1, 4):
                x = l[j]
                y = 0
                width = l[j+1] - x
                height = max(counts) if len(counts) != 0 else 1
                axs[i].add_patch(Rectangle((x, y), width, height, alpha=0.2, color="r"))
            
            if density_data_2:
                l = density_data_2[s]["x"].tolist()
                for j in range(2, len(l)-1, 4):
                    x = l[j]
                    y = 0
                    width = l[j+1] - x
                    height = max(counts) if len(counts) != 0 else 1
                    axs[i].add_patch(Rectangle((x, y), width, height, alpha=0.2, color="g"))
                
            axs[i].set_title(f"{s}")
            axs[i].set_xlim(left=0, right=max(density_data[s]["x"]))
            axs[i].set_ylim(bottom=0)
            axs[i].set_xlabel("sequence position")
            axs[i].set_ylabel("NGS count")

        save_path = os.path.join(RESULTSPATH, "NP_density", f"{k}_{author}_del_position_NP_density.pdf") # leave as .pdf, because saving as .png loses some bars
        plt.savefig(save_path)
        plt.close()

        NGS_dict[k] = count_dict

    return NGS_dict


def map_dens_to_dens(strain: str,
                     dens_1: dict,
                     dens_2: dict
                     )-> None:
    '''
        Plots the NP density of two different papers as source together.
        :param data: dict with information about start and end position
        :param dens_1: dict with dens data (key is segment name) from Lee
        :param dens_2: dict with dens data (key is segment name) from Le Sage

        :return: None
    '''
    for k, v in dens_1.items():
        fig, axs = plt.subplots(8, 1, figsize=(7, 14), tight_layout=True)
        for i, s in enumerate(SEGMENTS):
            l1, = axs[i].plot(dens_1[s]["x"], dens_1[s]["y"], label="Lee", alpha=0.5, color="green", fillstyle="full")
            l2, = axs[i].plot(dens_2[s]["x"], dens_2[s]["y"], label="Le Sage", alpha=0.5, color="blue", fillstyle="full")

            axs[i].set_title(f"{s}")
            axs[i].set_xlim(left=0, right=max(dens_1[s]["x"]))
            axs[i].set_ylim(bottom=0, top=105)
            axs[i].set_xlabel("sequence position")
            axs[i].set_ylabel("NP density")

        fig.suptitle("Mapping of the data source for NP density data")
        fig.legend([l1, l2], ["Lee", "Le Sage"])

        save_path = os.path.join(RESULTSPATH, "NP_density", f"{strain}_map_density_data_sources.png")
        plt.savefig(save_path)
        plt.close()


def compare_position_with_density(data: dict,
                                  density_data: dict,
                                  all_reads: dict,
                                  author: str=""
                                  )-> None:
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

    thresh_dict = dict({
        "PB2": 800, "PB1": 600, "PA": 700, "HA": 700, "NP": 600, "NA": 600, "M": 900, "NS": 500
    })

    obs_ratios = list()
    exp_ratios = list()
    symbols = list()

    SEGMENTS = ["PB2", "PB1", "PA", "HA"]

    plt.rc("font", size=18)
    for k, v in data.items():
        fig, ax = plt.subplots(1, 1, figsize=(10, 5), tight_layout=True)

        for i, s in enumerate(SEGMENTS):
            count_dict = v[s]
            if len(count_dict) == 0:
                obs_ratios.append(0.0)
                exp_ratios.append(0.0)
                symbols.append("")
                continue

            # get expected values by sampling approach
            seq = get_sequence(k, s)
            seg_reads = all_reads[k].loc[all_reads[k]["Segment"] == s]
            start = (int(seg_reads.Start.quantile(QUANT)), int(seg_reads.Start.quantile(1-QUANT)))
            end = (int(seg_reads.End.quantile(QUANT)), int(seg_reads.End.quantile(1-QUANT)))
            sampling_data = generate_sampling_data(seq, start, end, N_SAMPLES)
            samp_pos = sampling_data["Start"].tolist()
                
            rel_pos = count_dict.keys()
            rel_pos = [element for element in rel_pos if element < thresh_dict[s]]
            rel_pos_samp = [element for element in samp_pos if element < thresh_dict[s]]
            n = len(rel_pos)

            # grouping by NP high/low
            y = [get_dens_at_pos(p, density_data[s]) for p in rel_pos]
            y_exp = [get_dens_at_pos(p, density_data[s]) for p in rel_pos_samp]
            below = y.count(0)
            below_exp = y_exp.count(0)
            obs_ratio = below/len(y)
            exp_ratio = below_exp/len(y_exp)
          
            obs_ratios.append(obs_ratio)
            exp_ratios.append(exp_ratio)

            # statistical testing
            result = stats.binomtest(below, n, exp_ratio)
            symbol = get_stat_symbol(result.pvalue)
            ax.annotate(symbol, (i, max(obs_ratio, exp_ratio)), horizontalalignment="center")
            symbols.append(symbol)

        bar_width = 0.35
        x = np.arange(len(SEGMENTS))
        ax.bar(x - bar_width/2, obs_ratios, bar_width, label="expected")
        ax.bar(x + bar_width/2, exp_ratios, bar_width, label="observed")
        ax.set_xticks(x, SEGMENTS)
        ax.set_xlabel("Segment")
        ax.set_ylabel("$r_{NP}$")
        plt.legend(bbox_to_anchor=(1.0, 1.0))
        fig.suptitle(f"{STRAINS[k]}")

        savepath = os.path.join(RESULTSPATH, "NP_density", f"{k}_{author}_high_low_NP_areas.png")
        fig.savefig(savepath)
        plt.close()

    plot_data = dict({"obs": obs_ratios, "exp": exp_ratios, "symbol": symbols})

    return plot_data


def plot_ratios_together(data: dict)-> None:
    '''

    '''
    plt.rc("font", size=18)
    fig, ax = plt.subplots(1, 1, figsize=(10, 4), tight_layout=True)
    cm = plt.get_cmap('tab10')
    ax.set_prop_cycle('color', [cm(1.*i/10) for i in range(10)])
    for i, (k, v) in enumerate(data.items()):
        obs = v["obs"]
        exp = v["exp"]
        symbol = v["symbol"]

        bar_width = 0.05
        xs = np.arange(len(obs))
        move = i * 0.25 - 0.3 - (0.1*i)
        for x in xs:
            ax.annotate(symbol[x], (x - bar_width/2 + move + 0.025, max(obs[x], exp[x])), horizontalalignment="center", fontsize=8)

        ax.bar(xs - bar_width/2 + move, obs, bar_width, label=f"{k} exp.")
        ax.bar(xs + bar_width/2 + move, exp, bar_width, label=f"{k} obs.")

    SEGMENTS = ["PB2", "PB1", "PA", "HA"]

    ax.set_xticks(xs, SEGMENTS)
    ax.set_xlabel("Segment")
    ax.set_ylabel("$r_{NP}$")
    plt.legend(bbox_to_anchor=(1.0, 1.0))
    fig.suptitle("fraction of deletions in low NP areas")

    savepath = os.path.join(RESULTSPATH, "NP_density", f"all_high_low_NP_areas.png")
    fig.savefig(savepath)
    plt.close()


if __name__ == "__main__":
    plt.style.use('seaborn')
    ratio_data = dict()

    density_path = os.path.join(DATAPATH, "Lee2017", "csv_NPdensity")
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    del all_reads_dict["NC"]
    del all_reads_dict["Perth"]
    del all_reads_dict["BLEE"]
    
    # Plotting NP data (from Sage 2019) against junction sites
    #    Cal07 data from Alnaji 2019
    Cal07_dens_path = os.path.join(density_path, "Cal07")
    Cal07_dens_data = load_Sage_density_data(Cal07_dens_path)
    map_dens_to_dens("Cal07", load_Lee_density_data(Cal07_dens_path), Cal07_dens_data)
    NGS_count_dict = map_positions_to_density(all_reads_dict, Cal07_dens_data)
    ratio_data["Cal07"] = compare_position_with_density(NGS_count_dict, Cal07_dens_data, all_reads_dict)
    
    #    WSN data from Mendes 2021 or Boussier 2020
    source = "Mendes"
    WSN_reads_dict = load_WSN_data(source)
    WSN_dens_path = os.path.join(density_path, "WSN")
    WSN_dens_data = load_Sage_density_data(WSN_dens_path)
    map_dens_to_dens("WSN", load_Lee_density_data(WSN_dens_path), WSN_dens_data)
    WSN_NGS_count_dict = map_positions_to_density(WSN_reads_dict, WSN_dens_data)
    ratio_data["WSN Mendes"] = compare_position_with_density(WSN_NGS_count_dict, WSN_dens_data, WSN_reads_dict)
    
    source = "Boussier"
    WSN_reads_dict = load_WSN_data(source)   
    WSN_NGS_count_dict = map_positions_to_density(WSN_reads_dict, WSN_dens_data, author=source)
    ratio_data["WSN Boussier"] = compare_position_with_density(WSN_NGS_count_dict, WSN_dens_data, WSN_reads_dict, author=source)

    # NP data from Williams 2021 they define high and low areas
    #   PR8 data from Pelz and Alnaji2021
    PR8_dens_path = os.path.join(density_path, "PR8")
    PR8_high_dens_data = load_Williams_density_data(PR8_dens_path, "high")
    PR8_low_dens_data = load_Williams_density_data(PR8_dens_path, "low")
    PR8_reads_dict = load_full_alnaji2021()
    PR8_NGS_count_dict = map_positions_to_density(PR8_reads_dict, PR8_high_dens_data, PR8_low_dens_data)
    ratio_data["PR8 Alnaji"] = compare_position_with_density(PR8_NGS_count_dict, PR8_high_dens_data, PR8_reads_dict)

    pelz_data = load_pelz_dataset()
    PR8_NGS_count_dict = map_positions_to_density(pelz_data, PR8_high_dens_data, PR8_low_dens_data, author="Pelz")
    ratio_data["PR8 Pelz"] = compare_position_with_density(PR8_NGS_count_dict, PR8_high_dens_data, PR8_reads_dict, author="Pelz")

    plot_ratios_together(ratio_data)