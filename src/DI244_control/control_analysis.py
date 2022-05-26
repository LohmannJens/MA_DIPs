"""

"""
import os
import sys

from Bio import SeqIO

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS, load_excel, get_sequence
sys.path.insert(0, "../density_and_length_analysis")
from composition_junction_site import count_nucleotide_occurrence, calculate_overlapping_nucleotides

### parameters of DI244 (from Dimmock 2008) ###
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
    # load sequence
    # save as dict
    fasta_file = os.path.join(DATAPATH, "Dimmock2008", "PB2.fasta")
    DI244_dict = load_sequence_as_dict(fasta_file)

    excel_file = os.path.join(DATAPATH, "alnaji2019", "DI_Influenza_FA_JVI.xlsx")
    alnaji_dict = load_excel(excel_file)

    # do analysis for DI244 (control DI RNA)
    seq = DI244_dict["WholeSequence"]
    s = DI244_dict["Start"]
    e = DI244_dict["End"]
    check_deletion_site(seq, s, e)

    # do analysis for all fragments from alnaji that have NGS_count > 1000
    for k, v in alnaji_dict.items():
        print(k)
        for row in v.iterrows():
            r = row[1]
            if r.loc["NGS_read_count"] > 1000:
                record = get_sequence(k[:-3], r["Segment"])
                RNA_seq = record.seq.transcribe()
                check_deletion_site(str(RNA_seq), r["Start"], r["End"])





