'''
    Using the polymerase footprints as motifs and search if they occur more
    often in the sequences.
'''
import os
import sys

import pandas as pd

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
            d[f].append(seq.find(f))

    df = pd.DataFrame(d)

    path = os.path.join(RESULTSPATH, "motif_discovery", f"polymerase_footprint_{strain}.tex")
    df.to_latex(path)


if __name__ == "__main__":
#    strains = ["Cal07", "NC", "Perth", "BLEE", "PR8"]
 #   for strain in strains:
  #      search_footprint(strain)

    search_footprint("PR8")
