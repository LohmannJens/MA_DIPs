'''
    Using the polymerase footprints as motifs and search if they occur more
    often in the sequences.
'''
import os
import sys
import random

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from scipy import stats
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import get_sequence, load_alnaji_excel, load_short_reads, load_alnaji_2021


def search_footprint(strain: str)-> None:
    '''
        Takes the polymerase footprints and searches on each segment of a given
        strain for them. Creates a LaTeX table and saves the positions where 
        the motifs were found in it.
        :param strain: name of the strain

        :return: None
    '''
    footprints = ["AGCAAAAGCAGG", "AGCGAAAGCAGG", "CCUUGUUUCUACU"]
    d = dict({"segments" : SEGMENTS})
    for f in footprints:
        d[f] = list()
    for s in SEGMENTS:
        seq = get_sequence(strain, s)
        for f in footprints:
            d[f].append([i for i in range(len(seq)) if seq.startswith(f, i)])

    df = pd.DataFrame(d)
    path = os.path.join(RESULTSPATH, "motif_discovery", f"polymerase_footprint_{strain}.tex")
    df.to_latex(path)


def extended_footprint_search(strain: str)-> (float, float, float):
    '''
        Search for the match with the highest score for each footprint on the
        sequences of a given strain. Compares the score to randomly shuffled
        footprints.
        :param strain: name of the influenza strain

        :return: tuple with four entries
                 mean score of the matches
                 mean score of the matches of the randomly shuffled footprints
                 p-value of the t-test comparing both means
                 data frame including start and end positions for all matches
    '''
    # defining polymerase footprints
    f1 = "AGCAAAAGCAGG"
    f2 = "AGCGAAAGCAGG"
    f3 = "CCUUGUUUCUACU"
    footprints = [Seq(f1), Seq(f2), Seq(f3)]

    # generate random footprints by shuffling to have a comparision
    rand_footprints = list()
    for f in [f1, f2, f3]:
        rand_footprints = rand_footprints + [Seq("".join(random.sample(f,len(f)))) for i in range(50)]

    # start PairwiseAligner and set parameters that no gaps are allowed
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.target_internal_open_gap_score = -9999
    aligner.target_left_open_gap_score = -9999
    aligner.target_right_open_gap_score = -9999
    aligner.query_internal_open_gap_score = -9999
    aligner.query_left_open_gap_score = -9999
    aligner.query_right_open_gap_score = -9999

    # run alignment for all segments with given and random footprints
    scores = list()
    rand_scores = list()
    positions = dict({"Start":list(), "End":list(), "Motif":list(), "Segment":list()})
    for s in SEGMENTS:
        seq = get_sequence(strain, s)
        if strain == "PR8": # leave 100 percent matches out for PR8
            seq = seq[13:len(seq)-13]

        for f in footprints:
            for a in aligner.align(seq, f):
                scores.append(a.score)
                positions["Start"].append(a.path[0][0])
                positions["End"].append(a.path[1][0])
                positions["Motif"].append(f)
                positions["Segment"].append(s)

        for r_f in rand_footprints:
            for a in aligner.align(seq, r_f):
                rand_scores.append(a.score)

    # calculate mean and statistical testing with one sided t-test
    mean = sum(scores)/len(scores)
    rand_mean = sum(rand_scores)/len(rand_scores)
    stat, p = stats.ttest_1samp(scores, rand_mean, alternative="greater")

    return mean, rand_mean, p, pd.DataFrame(positions)


def plot_motif_positions_on_sequence(df: object, ngs_df: object, strain: str)-> None:
    '''
        Gets the positions of the found motif and plots them on the full length
        sequence together with the NGS count.
        :param df: data frame including the start and end of the motifs
        :param ngs_df: data frame including the start and end of the DI 
                       candidates and the NGS count
        :param strain: name of the influenza strain

        :return None:
    '''
    f_col = dict({"AGCAAAAGCAGG": "red",
                  "AGCGAAAGCAGG": "blue",
                  "CCUUGUUUCUACU": "green"})

    fig, axs = plt.subplots(8, 1, figsize=(7, 10), tight_layout=True)
    for i, seg in enumerate(SEGMENTS):
        s_df = df[df["Segment"] == seg]
        ngs_s_df = ngs_df[ngs_df["Segment"] == seg]

        if not ngs_s_df.empty:
            rect_h = max(ngs_s_df["NGS_read_count"])/10
            axs[i].bar(ngs_s_df["Start"], ngs_s_df["NGS_read_count"])
            axs[i].bar(ngs_s_df["End"], ngs_s_df["NGS_read_count"])
        else:
            rect_h = 1

        if not s_df.empty:
            for r in s_df.iterrows():
                row = r[1]
                s = row["Start"]
                e = row["End"]
                l = row["Motif"]
                p = patches.Rectangle((s, 0), e-s, rect_h, label=l, color=f_col[l])
                axs[i].add_patch(p)

        axs[i].set_title(f"{seg}")
        axs[i].set_xlabel("Sequence position")
        axs[i].set_ylabel("NGS count")

    plt.legend()
    path = os.path.join(RESULTSPATH, "motif_discovery", f"footprints_on_sequence_{strain}.pdf")
    plt.savefig(path)


if __name__ == "__main__":
    strains = ["Cal07", "NC", "Perth", "BLEE", "PR8"]
    ngs_data = load_alnaji_2021()
    dict2 = load_alnaji_excel()
    dict2 = load_short_reads(dict2)
    ngs_data.update(dict2)

    d = dict({"header": ["mean", "r_mean", "p"]})
    for strain in strains:
        search_footprint(strain)
        m, rand_m, p, positions_df = extended_footprint_search(strain)
        d[strain] = [m, rand_m, p]

        plot_motif_positions_on_sequence(positions_df, ngs_data[strain], strain)

    df = pd.DataFrame(d)
    path = os.path.join(RESULTSPATH, "motif_discovery", f"polymerase_footprint_extended.tex")
    df.to_latex(path)

