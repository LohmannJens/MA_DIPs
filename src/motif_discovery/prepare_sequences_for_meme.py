'''
    This script creates FASTA files that will be used for motif discovery by
    MEME suite.
'''
import os
import sys
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, "..")
from utils import DATAPATH, SEGMENTS
from utils import load_alnaji_excel, load_short_reads, get_sequence, create_sequence_library


def delete_folder(folder: str)-> bool:
    '''
        Asks the user if a given folder should be overwritten.
        :param folder: folder which will be overwritten

        :return: True if folder should be deleted; False otherwise
    '''
    print(f"{folder} already exists! Should it be overwritten?")
    if input("[Y/n]: ") == "Y":
        return True
    else:
        return False


def write_sequence(record: SeqRecord,
                   strain: str,
                   segment: str,
                   folder: str,
                   control: bool=False
                   )-> None:
    '''
        Gets RNA sequence as Biopython SeqRecord and writes it into three
        files. One for all sequences, one for the strains and one for the
        different segments.
        :param seq: RNA sequence as Biopython SeqRecord
        :param strain: name of the strain
        :param segment: name of the RNA segment
        :param folder: path to the save location

        :return: None
    '''
    if control:
        all_file = os.path.join(folder, "all_control.fasta")
        strain_file = os.path.join(folder, f"{strain}_control.fasta")
        seg_file = os.path.join(folder, f"{segment}_control.fasta")
    else:
        all_file = os.path.join(folder, "all.fasta")
        strain_file = os.path.join(folder, f"{strain}.fasta")
        seg_file = os.path.join(folder, f"{segment}.fasta")

    with open(all_file, "a") as f:
        SeqIO.write(record, f, "fasta")
    with open(strain_file, "a") as f:
        SeqIO.write(record, f, "fasta")
    with open(seg_file, "a") as f:
        SeqIO.write(record, f, "fasta")


def create_full_seq_files(strains: list)-> None:
    '''
        Creates FASTA files for the different strains and segments.
        :param strains: list of the different virus strains

        :return: None
    '''
    root_folder = os.path.join(DATAPATH, "meme_suite", "alnaji2019", "full_sequences")
    if os.path.exists(root_folder):
        if delete_folder(root_folder):
            shutil.rmtree(root_folder)
        else:
            return
    os.makedirs(root_folder)

    for st in strains:
        for s in SEGMENTS:
            seq = get_sequence(st, s, full=True)
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
        if delete_folder(root_folder):
            shutil.rmtree(root_folder)
        else:
            return
    os.makedirs(root_folder)

    for k, v in d.items():
        for r in v.iterrows():
            r = r[1]
            seq = Seq(r["DIRNASequence"])
            seg = r["Segment"]
            s = r["Start"]
            e = r["End"]
            record = SeqRecord(seq, id=f"{k}_{seg}_{s}_{e}")
            write_sequence(record, k, seg, root_folder)


def create_windows_del_site_files(d: dict,
                                  n: int,
                                  combine: bool,
                                  only_remain: bool
                                  )-> None:
    '''
        Creates FASTA files for a n wide window around the start and end of the
        sequences of the different strains and segments. 
        :param d: dict containing sequence and deletion site info
        :param n: half window size (only indicating size in one direction)
        :param combine: states if the start and end window should be
                        concatenated or not
        :param only_remain: states if only the first and fourth part around 
                            the deletion site should be kept

        :return: None
    '''
    fname = f"window_{n}_sequences"
    if combine:
        fname = f"{fname}_combined"
    elif only_remain:
        fname = f"{fname}_onlyremain"
    root_folder = os.path.join(DATAPATH, "meme_suite", "alnaji2019", fname)
    if os.path.exists(root_folder):
        if delete_folder(root_folder):
            shutil.rmtree(root_folder)
        else:
            return
    os.makedirs(root_folder)

    for k, v in d.items():
        for row in v.iterrows():
            r = row[1]
            seg = r["Segment"]
            seq = get_sequence(k, seg)
            s = r["Start"]
            e = r["End"]
            # cropping the two windows at the deletion site and concatenating
            # them. Min/Max operator to avoid conflicts when site is close to
            # the beginning/end of the whole sequence.
            if combine:
                window_seq = seq[max(s-n, 0):s+n] + seq[e-n:min(e+n, len(seq))]
                record = SeqRecord(Seq(window_seq), id=f"{k}_{seg}_{s}_start_{e}_end")
                write_sequence(record, k, seg, root_folder)
            elif only_remain:
                print(s, n, e)
                window_seq = seq[max(s-n, 0):s] + seq[e:min(e+n, len(seq))]
                print(window_seq)
                record = SeqRecord(Seq(window_seq), id=f"{k}_{seg}_{s}_start_{e}_end")
                write_sequence(record, k, seg, root_folder)
            else:
                window_seq = seq[max(s-n, 0):s+n]
                record = SeqRecord(Seq(window_seq), id=f"{k}_{seg}_{s}_start")
                write_sequence(record, k, seg, root_folder)

                window_seq = seq[e-n:min(e+n, len(seq))]
                record = SeqRecord(Seq(window_seq), id=f"{k}_{seg}_{e}_end")
                write_sequence(record, k, seg, root_folder)


if __name__ == "__main__":
    cleaned_data_dict = load_alnaji_excel()
    all_reads_dict = load_short_reads(cleaned_data_dict)
    create_full_seq_files(list(all_reads_dict.keys()))

    seq_library = create_sequence_library(all_reads_dict)
    create_cropped_seq_files(seq_library)
    create_windows_del_site_files(seq_library, 50, False, True)

