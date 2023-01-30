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
from utils import create_sequence_library, get_sequence
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
        seq = seq.translate(to_stop=True)
        if type(r["Class"]) != "str":
            desc = ""
        else:
            desc = r["Class"]
        record = SeqRecord(seq, id, description=desc)
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
     #       newline = newline.replace("*", "")
            new_lines.append(newline)

        with open(file, "w") as handle:
            handle.writelines(new_lines)

    
def cropped_sequence_library(df: pd.DataFrame)-> pd.DataFrame:
    '''
        Uses the specified start and end points for translation and crops the
        sequences corresponding to that.
        :param df: data frame including the candidates

        :return: data frame with the cropped sequences in column 'DelSequence'
    '''
    indices = dict({"PB2": (28, 2307),
                    "PB1": (25, 2298),
                    "PA":  (25, 2175)
                    })

    del_seq_list = list()
    for i, row in df.iterrows():
        full_seq = get_sequence("PR8", row["Segment"])
        s = indices[row["Segment"]][0]-1
        e = indices[row["Segment"]][1]
        del_seq = full_seq[s:row["Start"]] + full_seq[row["End"]-1:e]
        del_seq_list.append(del_seq)

    df["DelSequence"] = del_seq_list

    # add full sequences
    for k, v in indices.items():
        s = v[0]-1
        e = v[1]
        del_seq = full_seq[s:e]
        
        df.loc[len(df)] = list([f"{k}_full", k, pd.NA, pd.NA, "full", del_seq])

    return df


if __name__ == "__main__":
    fasta_file = os.path.join(DATAPATH, "Dimmock2008", "PB2.fasta")

    pelz_top_df = load_top_DI_RNA_pelz()
    di244_df = load_DI244()

    pelz_top_df = pelz_top_df.drop("NGS_read_count", axis=1)
    full_df = pd.concat([pelz_top_df, di244_df], ignore_index=True)
    name_col = full_df.pop("Name")
    full_df.insert(0, "Name", name_col)
#    seq_df = create_sequence_library({"PR8": full_df})["PR8"]
    seq_df = cropped_sequence_library(full_df)

    generate_fastas(seq_df)
    clean_output_files()

