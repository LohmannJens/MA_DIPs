'''
    Using the polymerase footprints as motifs and search if they occur more
    often in the sequences.
'''
import os
import sys

import pandas as pd

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import get_sequence


def search_footprint(strain: str)-> None:
    '''
        Takes the polymerase footprints and searches on each segment of a given
        strain for them. Creates a LaTeX table and saves the positions where 
        the motifs were found in it.
        :param strain: name of the strain

        :return: None
    '''
    footprints = ["AGCAAAAGCAGG", "AGCGAAAGCAGG", "CCUUGUUUCUACU"]
    d = dict({"segments" : SEGMENTS})
    for f in footprints:
        d[f] = list()
    for s in SEGMENTS:
        seq = get_sequence(strain, s)
        for f in footprints:
            d[f].append([i for i in range(len(seq)) if seq.startswith(f, i)])

    df = pd.DataFrame(d)

    path = os.path.join(RESULTSPATH, "motif_discovery", f"polymerase_footprint_{strain}.tex")
    df.to_latex(path)


def extended_footprint_search(strain: str)-> None:
    '''

    '''
    footprints = [Seq("AGCAAAAGCAGG"), Seq("AGCGAAAGCAGG"), Seq("CCUUGUUUCUACU")]

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.target_internal_open_gap_score = -99999
    aligner.target_left_open_gap_score = -99999
    aligner.target_right_open_gap_score = -99999
    aligner.query_internal_open_gap_score = -99999
    aligner.query_left_open_gap_score = -99999
    aligner.query_right_open_gap_score = -99999

    for s in SEGMENTS:
        seq = get_sequence(strain, s)
        for f in footprints:
            for a in aligner.align(seq, f):
                print(a.aligned, a.score)

if __name__ == "__main__":
    strains = ["Cal07", "NC", "Perth", "BLEE", "PR8"]
    for strain in strains:
        print(strain)
       # search_footprint(strain)
        extended_footprint_search(strain)


    extended_footprint_search("PR8")
