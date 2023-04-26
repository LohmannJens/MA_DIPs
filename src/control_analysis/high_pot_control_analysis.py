'''
    Checks my overall findings by comparing to high potential DI RNA. They are
    taken from Pelz et al. 2021 and the DI244 DI RNA is used.
'''
import os
import sys

import pandas as pd

from Bio import SeqIO

sys.path.insert(0, "..")
sys.path.insert(0, "../relative_occurrence_nucleotides")
sys.path.insert(0, "../direct_repeats")
from utils import DATAPATH, RESULTSPATH
from utils import load_alnaji_excel, get_sequence, get_seq_len
from composition_junction_site import count_nucleotide_occurrence
from search_direct_repeats import calculate_direct_repeat


def load_DI244()-> object:
    '''
        Loads the data of the DI244 candidate

        :return: pandas data frame
    '''
    data = dict({"Name": "DI244",
                 "Segment": "PB2",
                 "Start": 244,
                 "End": 2191
                })
    df = pd.DataFrame(data, index=[0])
    return df


def load_top_DI_RNA_pelz()-> object:
    '''
        Loads the top gain DI RNAs and top de novo DI RNAs from Pelz et al 2021
        
        :return: pandas data frame
    '''
    def create_name(r):
        return "_".join((r["Segment"], str(r["Start"]), str(r["End"])))

    path = os.path.join(DATAPATH, "Pelz2021", "Top_DI_RNA.xlsx")
    data = pd.read_excel(path)
    data["Name"] = data.apply(create_name, axis=1)
    return data


def get_sequence_around_site(segments: list,
                             points: list
                             )-> list:
    '''
        Returns the part of the sequences that are directly around the deletion
        site.
        :param segments: list of segments
        :param points: list of integers indicating the deletion site, is of
                       equal length than segments

        :return: list of sequences around the deletion site
    '''
    sequences = list()
    for s, p in zip(segments, points):
        seq = get_sequence("PR8", s)
        sequences.append(seq[p-5:p+4])
    return sequences


def analyse_top_DI_RNA(df: pd.DataFrame)-> None:
    '''
        Gets a data frame of DI RNA and calculates the lengths of the resulting
        fragments. Also does the sequence overlap analysis in the two different
        modes. Saves the results as .csv file and .tex table
        :param df: pandas data frame with data of DI RNA

        :return: None
    '''
    # calculate length of resulting DI RNA
    for s in ["PB2", "PB1", "PA"]:
        length = df["Start"] + (get_seq_len("PR8", s) - df["End"] + 1)
        df.loc[df["Segment"] == s, "Length"] = length
    df["Length"] = df["Length"].astype(int)
    
    # nucleotide occurrence at junction start and end site
    segments = df["Segment"].tolist()
    starts = df["Start"].tolist()
    ends = df["End"].tolist()
    df["Sequence_At_Start"] = get_sequence_around_site(segments, starts)
    df["Sequence_At_End"] = get_sequence_around_site(segments, ends)

    # calculate direct repeats and resulting sequence
    direct_repeat_m1 = list()
    sequence_m1 = list()
    direct_repeat_m2 = list()
    sequence_m2 = list()
    for _, (_, seg, s, e, _, _, _, _) in df.iterrows():
        seq = get_sequence("PR8", seg)
        c, overlap_seq = calculate_direct_repeat(seq, s, e, w_len=10, m=1)
        direct_repeat_m1.append(c)
        sequence_m1.append(overlap_seq)
        c, overlap_seq = calculate_direct_repeat(seq, s, e, w_len=10, m=2)
        direct_repeat_m2.append(c)
        sequence_m2.append(overlap_seq)

    df["Direct_Repeat_Length_Mode1"] = direct_repeat_m1
    df["Direct_Repeat_Sequence_Mode1"] = sequence_m1
    df["Direct_Repeat_Lenght_Mode2"] = direct_repeat_m2
    df["Direct_Repeat_Sequence_Mode2"] = sequence_m2

    df = df.set_index("Name")

    # output table as latex and csv file
    path = os.path.join(RESULTSPATH, "control_analysis")
    df.to_latex(os.path.join(path, "top_DI_RNA.tex"))
    df.to_csv(os.path.join(path, "top_DI_RNA.csv"))


if __name__ == "__main__":
    fasta_file = os.path.join(DATAPATH, "Dimmock2008", "PB2.fasta")

    pelz_top_df = load_top_DI_RNA_pelz()
    di244_df = load_DI244()

    pelz_top_df = pelz_top_df.drop("NGS_read_count", axis=1)
    full_df = pd.concat([pelz_top_df, di244_df], ignore_index=True)
    name_col = full_df.pop("Name")
    full_df.insert(0, "Name", name_col)

    analyse_top_DI_RNA(full_df)

