"""

"""
import os
import sys

from Bio import SeqIO

sys.path.insert(0, "..")
sys.path.insert(0, "../density_and_length_analysis")
from utils import RESULTSPATH, SEGMENTS
from utils import load_pelz_dataset, get_sequence
from composition_junction_site import create_sequence_library


if __name__ == "__main__":
    data_dict = load_pelz_dataset()

    seq_list_dict = create_sequence_library(data_dict)

