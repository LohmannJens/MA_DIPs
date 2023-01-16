'''
    generates FASTA files of the high pot candidates from Pelz and DI244.
'''
import os
import sys

import pandas as pd

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH
from utils import create_sequence_library
from high_pot_control_analysis import load_DI244, load_top_DI_RNA_pelz

sys.path.insert(0, "../free_energy_estimations")
from prepare_sequences_for_RNAvienna import write_sequence


def generate_fastas(df: pd.DataFrame)-> None:
    '''
        Generates an independent FASTA file for each entry in the given Data
        Frame
        :param df: data frame including at least a Name and DelSequence

        :return: None
    '''
    root_folder = os.path.join(DATAPATH, "FIRM-AVP")
    for i, r in df.iterrows():
        id = r["Name"]
        seq = Seq(r["DelSequence"])
        seq = seq.reverse_complement().translate()
        record = SeqRecord(seq, id)

        write_sequence(record, id, root_folder)


def clean_output_files()-> None:
    '''
        In each of the generated FASTA files the sequence is merged into one
        column by removing the newlines and the stop codons are removed.

        :return: None
    '''
    # open each file in given folder
    root_folder = os.path.join(DATAPATH, "FIRM-AVP")
    for f in os.listdir(root_folder):
        file = os.path.join(root_folder, f)
        with open(file) as handle:
            lines = handle.readlines()

            new_lines = list()
            newline = ""
            for l in lines:
                if l.startswith(">"):
                    new_lines.append(l)
                else:
                    newline = newline + l.replace("\n", "")

            newline = newline + "\n"
            newline = newline.replace("*", "")
            new_lines.append(newline)

        with open(file, "w") as handle:
            handle.writelines(new_lines)


if __name__ == "__main__":
    fasta_file = os.path.join(DATAPATH, "Dimmock2008", "PB2.fasta")

    pelz_top_df = load_top_DI_RNA_pelz()
    di244_df = load_DI244()

    pelz_top_df = pelz_top_df.drop("NGS_read_count", axis=1)
    full_df = pd.concat([pelz_top_df, di244_df], ignore_index=True)
    name_col = full_df.pop("Name")
    full_df.insert(0, "Name", name_col)
    seq_df = create_sequence_library({"PR8": full_df})

    generate_fastas(seq_df["PR8"])
    clean_output_files()

