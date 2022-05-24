"""

"""
import os
import sys

from Bio import SeqIO

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, load_excel, load_short_reads, get_sequence
sys.path.insert(0, "../density_and_length_analysis")
from composition_junction_site import count_nucleotide_occurrence, calculate_overlapping_nucleotides

START = 244
END = 2191
LENGTH = 2341


def load_sequence_as_dict(path: str)-> object:
    '''

    '''
    record = SeqIO.read(path, "fasta")
    dna_seq = record.seq
    rna_seq = dna_seq.transcribe()

    data = dict({"ID": record.id,
                 "Length": LENGTH,
                 "Start": START,
                 "End": END,
                 "WholeSequence": rna_seq,
                 "DelSequence": rna_seq[:START] + rna_seq[END-1:]})

    return data


if __name__ == "__main__":
    # load sequence
    # save as dict
    fasta_file = os.path.join(DATAPATH, "Dimmock2008", "PB2.fasta")
    data = load_sequence_as_dict(fasta_file)

    # count nucleotides before and after junction site with
    start = count_nucleotide_occurrence(data["WholeSequence"], data["Start"])
    end = count_nucleotide_occurrence(data["WholeSequence"], data["End"])

    print(start, end)
    print("1234J6789\t\t\tJ is position 244 (last one in sequence that is not deleted)")
    print(data["WholeSequence"][START-5:START+4])
    print("1234J6789\t\t\tJ is position 2191 (first one after deletion)")
    print(data["WholeSequence"][END-5:END+4])

    # calculate number of overlapping nucleotides
    # gives numver of overlapping nt and sequence of overlap
    count, overlap_seq = calculate_overlapping_nucleotides(data["WholeSequence"], data["Start"], data["End"], w_len=10, m=1)

    print(count, overlap_seq)

    count, overlap_seq = calculate_overlapping_nucleotides(data["WholeSequence"], data["Start"], data["End"], w_len=10, m=2)

    print(count, overlap_seq)



