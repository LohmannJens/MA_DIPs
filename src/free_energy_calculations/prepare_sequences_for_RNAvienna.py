'''
    This script creates FASTA files that will be used for free energy
    estimations by RNAvienna.
'''
import os
import sys
import shutil

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randrange

sys.path.insert(0, "..")
sys.path.insert(0, "../density_and_length_analysis")
from utils import DATAPATH, SEGMENTS
from utils import load_alnaji_excel, load_short_reads, get_sequence
from composition_junction_site import create_sequence_library


def write_sequence(seq, name: str, folder: str)-> None:
    '''
        gets RNA sequence as Biopython SeqRecord and writes it into three
        files. One for all sequences, one for the strains and one for the
        different segments.
        :param seq: RNA sequence as Biopython SeqRecord
        :param name: filename of the FASTA file to be created; looks like:
                        "{strain}_{segment}_{start}_{end}"
        :param folder: path to the save location

        :return: None
    '''
    filename = os.path.join(folder, f"{name}.fasta")

    with open(filename, "w") as f:
        SeqIO.write(seq, f, "fasta")


def random_crop_sequence(s: str, n: int)-> str:
    '''
        Creates a random deletion site by cropping a given sequence at a random
        point.
        :param s: Given sequence
        :param n: length of the resulting sequence

        :return: randomly cropped sequence of length n
    '''
    start = randrange(len(s)-n)
    seq = s[:start] + s[start+n:]
    return seq


def create_cropped_seq_files(d: dict, shuffle: bool=False, random: bool=False)-> None:
    '''
        Creates FASTA files for the cropped sequences of the different strains
        and segments. Cropped sequences are the ones that exclude the deletion
        site of the DI RNA.
        :param d: dict containing sequence and deletion site info
        :param shuffle: indicates if the sequence should be shuffled. Is used
                        to generate a comparision
        :param random: indicates if a randomly cropped sequence should be
                       generated

        :return: None
    '''
    folder = "cropped_sequences"
    if shuffle:
        folder = f"{folder}_shuffled"
    elif random:
        folder = f"{folder}_randomcrop"
    root_folder = os.path.join(DATAPATH, "energy_calculation", folder)
    if os.path.exists(root_folder):
        print(f"{root_folder} already exists! Should it be overwritten?")
        if input("[Y/n]: ") == "Y":
            shutil.rmtree(root_folder)
        else:
            return
    os.makedirs(root_folder)

    for k, v in d.items():
        for r in v.iterrows():
            r = r[1]
            seg = r["Segment"]
            s = r["Start"]
            e = r["End"]
            NGS = r["NGS_read_count"]
            if k == "B_LEE":
                k = "BLEE"
            id = f"{k}_{seg}_{s}_{e}_{NGS}"

            if shuffle:
                seq_list = list(r["DelSequence"])
                np.random.shuffle(seq_list)
                seq = Seq("".join(seq_list))
                id = f"{id}_shuffled"
            elif random:
                full_seq = get_sequence(k, seg)
                del_len = len(full_seq) - len(r["DelSequence"])
                seq = random_crop_sequence(full_seq, del_len)
                id = f"{id}_randomcrop"
            else:
                seq = r["DelSequence"]
                
            record = SeqRecord(seq, id=id)
            write_sequence(record, id, root_folder)


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)

    seq_library = create_sequence_library(all_reads_dict)
    create_cropped_seq_files(seq_library)
    create_cropped_seq_files(seq_library, shuffle=True)
    create_cropped_seq_files(seq_library, random=True)

