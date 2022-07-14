"""

"""
import os
import sys

from Bio import SeqIO

sys.path.insert(0, "..")
sys.path.insert(0, "../density_and_length_analysis")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import load_alnaji_excel, get_sequence
from composition_junction_site import count_nucleotide_occurrence, calculate_overlapping_nucleotides

### parameters of DI244 (from Dimmock 2008) ###
START = 244
END = 2191
LENGTH = 2341


def load_sequence_as_dict(path: str)-> object:
    '''

    '''
    record = SeqIO.read(path, "fasta")
    rna_seq = str(record.seq.transcribe())

    data = dict({"ID": record.id,
                 "Length": LENGTH,
                 "Start": START,
                 "End": END,
                 "DelSequence": rna_seq[:START] + rna_seq[END-1:]})

    return data


def check_deletion_site(seq, s, e)-> None:
    '''

    '''
    # print out deletion site
    print("\n#####")
    print(f"1234J6789\t\t\tJ is position {s} (last one in sequence that is not deleted)")
    print(seq[s-5:s+4])
    print(f"1234J6789\t\t\tJ is position {e} (first one after deletion)")
    print(seq[e-5:e+4])

    # check the nucleotide overlap at start and end
    count, overlap_seq = calculate_overlapping_nucleotides(seq, s, e, w_len=10, m=1)
    print(count, overlap_seq)
    count, overlap_seq = calculate_overlapping_nucleotides(seq, s, e, w_len=10, m=2)
    print(count, overlap_seq)

    return


if __name__ == "__main__":
    fasta_file = os.path.join(DATAPATH, "Dimmock2008", "PB2.fasta")
    DI244_dict = load_sequence_as_dict(fasta_file)

    alnaji_dict = load_alnaji_excel()

    # do analysis for DI244 (control DI RNA)
    seq = str(SeqIO.read(fasta_file, "fasta").seq.transcribe())
    s = DI244_dict["Start"]
    e = DI244_dict["End"]
    check_deletion_site(seq, s, e)

    # do analysis for all fragments from alnaji that have NGS_count > 1000
    for k, v in alnaji_dict.items():
        print(k)
        for row in v.iterrows():
            r = row[1]
            if r.loc["NGS_read_count"] > 1000:
                RNA_seq = get_sequence(k[:-3], r["Segment"])
                check_deletion_site(str(RNA_seq), r["Start"], r["End"])

