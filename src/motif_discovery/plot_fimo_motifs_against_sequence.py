'''
    Takes the fimo results from xstreme run and plots the found motifs together
    with the start and end positions of the junction sites.
'''
import os
import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib import cm

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import load_alnaji_excel, load_short_reads, get_sequence, join_lineages, create_sequence_library


def get_sequence_id(strain: str,
                    seg: str
                    )-> str:
    '''
        Returns the sequence ID/name of a given segment of a strain.
        :param strain: name of the strain
        :param seg: name of the segment

        :return: ID of the segement of the strain
    '''
    record = get_sequence(strain, seg, full=True)
    return record.id


def get_fimo_folders(folder: str)-> list:
    '''
        Loads all folders that include fimo results of a specific run of
        xstreme.
        :param folder: path to the folder with the xstreme results

        :return: list all paths with fimo results
    '''
    subfolders = [f.path for f in os.scandir(folder) if f.is_dir() and f.name[:4] == "fimo"]
    return subfolders


def load_fimo_files(folders: list):
    '''
        Gets a list of paths with fimo .tsv files. Loads each of them and 
        merges them together in one pandas data frame.
        :param folders: list of paths to fimo folders

        :return: pandas data frame with the fimo data (motif name, sequence ID,
                 start and end of motif on sequence)
    '''
    li = list()
    for p in folders:
        f = os.path.join(p, "fimo.tsv")
        df = pd.read_csv(f, index_col=None, header=0, skipfooter=3, engine="python", sep="\t")
        li.append(df)

    final_df = pd.concat(li, axis=0, ignore_index=True)
    return final_df


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="plot position of xstreme motifs agains start and end points of junction site")
    p.add_argument("--cropped", "-c", action="store_true")
    p.add_argument("--data", "-d", type=str, help="define which data to use; should be 'all', 'IVA', or 'seg1-3'")
    p.add_argument("--weighted", "-w", action="store_true")
    args = p.parse_args()

    if args.data == "seg1-3":
        SEGMENTS = SEGMENTS[:3]
        fsize = (10, 7)
        plt.rc("font", size=14)
    else:
        fsize = (10, 10)

    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = join_lineages(cleaned_data_dict)
    
    # get all fimo files
    seq_folder = "cropped_sequences" if args.cropped else "full_sequences"
    xstreme_folder = f"{args.data}_xstreme"
    fimo_path = os.path.join(DATAPATH, "meme_suite", "alnaji2019", seq_folder, xstreme_folder)
    fimo_folders = get_fimo_folders(fimo_path)
    fimo_df = load_fimo_files(fimo_folders)

    if args.data == "IVA":
        del all_reads_dict["BLEE"]

    # have 4 figures (for each strain) including 8 subplots (for each segment)    
    for k, v in all_reads_dict.items():
        fig, axs = plt.subplots(len(SEGMENTS), 1, figsize=fsize, tight_layout=True)

        # create color labels for motifs
        color_labels = fimo_df["motif_alt_id"].unique()
        viridis = cm.get_cmap("gist_ncar", len(color_labels))
        color_values = viridis(np.linspace(0,1,len(color_labels)))
        color_map = dict(zip(color_labels, color_values))

        for i, seg in enumerate(SEGMENTS):
            df = v.loc[v["Segment"] == seg]
            id = get_sequence_id(k, seg)
            if not df.empty:
                s = df["Start"]
                e = df["End"]
                if args.weighted:
                    rect_height = max(df["NGS_read_count"])/10
                    axs[i].bar(df["Start"], df["NGS_read_count"])
                    axs[i].bar(df["End"], df["NGS_read_count"])
                else:    
                    # plot start and end point of junctions as scatter
                    rect_height = 0.5
                    axs[i].scatter(x=s, y=np.zeros(len(s)), marker="|")
                    axs[i].scatter(x=e, y=np.zeros(len(e)), marker="|")

            # plot the found motifs as rectangles
            sliced_fimo_df = fimo_df.loc[fimo_df["sequence_name"] == id]
            for r in sliced_fimo_df.iterrows():
                row = r[1]
                s = row["start"]
                e = row["stop"]
                label = row["motif_alt_id"]
                axs[i].add_patch(patches.Rectangle((s, 0), e-s, rect_height, linewidth=0.5, label=label, color=color_map[label]))

            axs[i].set_title(f"{seg}")
            axs[i].set_xlabel("Sequence position")
            axs[i].set_ylabel("NGS count")

        by_label = dict()
        for ax in axs:
            handles, labels = ax.get_legend_handles_labels()
            by_label.update(dict(zip(labels, handles)))

        fig.legend(by_label.values(), by_label.keys(), ncol=6, mode="expand")
        fig.suptitle(f"\n\n\n\n\n{k}")
        fig.subplots_adjust(top=0.2)
        if args.weighted:
            filename = f"{k}_{args.data}_weighted_motif_on_sequence.pdf"
        else:
            filename = f"{k}_{args.data}_motif_on_sequence.png"
        savepath = os.path.join(RESULTSPATH, "motif_discovery", filename)
        plt.savefig(savepath)

