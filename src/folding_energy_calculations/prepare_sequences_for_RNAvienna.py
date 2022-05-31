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


def create_cropped_seq_files(d: dict)-> None:
    '''
        Creates FASTA files for the cropped sequences of the different strains
        and segments. Cropped sequences are the ones that exclude the deletion
        site of the DI RNA.
        :param d: dict containing sequence and deletion site info

        :return: None
    '''
    root_folder = os.path.join(DATAPATH, "energy_calculation", "cropped_sequences")
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
            id = f"{k}_{seg}_{s}_{e}"
            record = SeqRecord(seq, id=id)
            write_sequence(record, id, root_folder)


if __name__ == "__main__":
    filepath = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    cleaned_data_dict = load_excel(filepath)

    short_reads_filepath = os.path.join(DATAPATH, "alnaji2019", "Small_deletionSize_FA.xlsx")
    all_reads_dict = load_short_reads(cleaned_data_dict, short_reads_filepath)

    seq_library = create_sequence_library(all_reads_dict)
    create_cropped_seq_files(seq_library)

