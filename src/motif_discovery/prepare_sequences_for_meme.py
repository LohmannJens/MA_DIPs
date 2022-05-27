"""
This script creates FASTA files that will be used for motif discovery by MEME
suite
"""
import os
import sys
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, "..")
sys.path.insert(0, "../density_and_length_analysis")
from utils import DATAPATH, SEGMENTS, load_excel, load_short_reads, get_sequence
from composition_junction_site import create_sequence_library


def write_sequence(seq, strain: str, segment: str, folder: str)-> None:
    '''
        gets RNA sequence as Biopython SeqRecord and writes it into three
        files. One for all sequences, one for the strains and one for the
        different segments.
        :param seq: RNA sequence as Biopython SeqRecord
        :param strain: name of the strain
        :param segment: name of the RNA segment
        :param folder: path to the save location

        :return: None
    '''
    all_file = os.path.join(folder, "all.fasta")
    strain_file = os.path.join(folder, f"{strain}.fasta")
    seg_file = os.path.join(folder, f"{segment}.fasta")

    with open(all_file, "a") as f:
        SeqIO.write(seq, f, "fasta")
    with open(strain_file, "a") as f:
        SeqIO.write(seq, f, "fasta")
    with open(seg_file, "a") as f:
        SeqIO.write(seq, f, "fasta")


def create_full_seq_files(strains: list)-> None:
    '''
        Creates FASTA files for the different strains and segments.
        :param strains: list of the different virus strains

        :return: None
    '''
    root_folder = os.path.join(DATAPATH, "meme_suite", "alnaji2019", "full_sequences")
    if os.path.exists(root_folder):
        print(f"{root_folder} already exists! Should it be overwritten?")
        if input("[Y/n]: ") == "Y":
            shutil.rmtree(root_folder)
        else:
            return
    os.makedirs(root_folder)

    for st in strains:
        for s in SEGMENTS:
            seq = get_sequence(st, s)
            write_sequence(seq, st, s, root_folder)


def create_cropped_seq_files(d: dict)-> None:
    '''
        Creates FASTA files for the cropped sequences of the different strains
        and segments. Cropped sequences are the ones that exclude the deletion
        site of the DI RNA.
        :param d: dict containing sequence and deletion site info

        :return: None
    '''
    root_folder = os.path.join(DATAPATH, "meme_suite", "alnaji2019", "cropped_sequences")
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
            seq = r["DelSequence"]
            seg = r["Segment"]
            s = r["Start"]
            e = r["End"]
            record = SeqRecord(seq, id=f"{k}_{seg}_{s}_{e}")
            write_sequence(record, k, seg, root_folder)


def create_windows_del_site_files(d: dict, n: int)-> None:
    '''
        Creates FASTA files for a n wide window around the start and end of the
        sequences of the different strains and segments. 
        :param d: dict containing sequence and deletion site info
        :param n: half window size (only indicating size in one direction)

        :return: None
    '''
    root_folder = os.path.join(DATAPATH, "meme_suite", "alnaji2019", f"window_{n}_sequences")
    if os.path.exists(root_folder):
        print(f"{root_folder} already exists! Should it be overwritten?")
        if input("[Y/n]: ") == "Y":
            shutil.rmtree(root_folder)
        else:
            return
    os.makedirs(root_folder)

    for k, v in d.items():
        for row in v.iterrows():
            r = row[1]
            seq = r["WholeSequence"]
            seg = r["Segment"]
            s = r["Start"]
            e = r["End"]
            # cropping the two windows at the deletion site and concatenating
            # them. Min/Max operator to avoid conflicts when site is close to
            # the beginning/end of the whole sequence.
            window_seq = seq[max(s-n, 0):s+n] + seq[e-n:min(s+n, len(seq))]
            record = SeqRecord(window_seq, id=f"{k}_{seg}_{s}_{e}")
            write_sequence(record, k, seg, root_folder)


if __name__ == "__main__":
    filepath = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    cleaned_data_dict = load_excel(filepath)

    short_reads_filepath = os.path.join(DATAPATH, "alnaji2019", "Small_deletionSize_FA.xlsx")
    all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)

    create_full_seq_files(list(all_reads_dict.keys()))

    seq_library = create_sequence_library(all_reads_dict)
    create_cropped_seq_files(seq_library)

    create_windows_del_site_files(seq_library, 50)
