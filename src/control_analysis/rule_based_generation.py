'''
    generates a .csv file for potential good PR8 candidates by applying rules.
'''
import os
import sys
import json
import argparse

import numpy as np
import pandas as pd

sys.path.insert(0, "..")
from utils import DATAPATH
from utils import load_alnaji_2021, get_sequence, get_seq_len

sys.path.insert(0, "../ML")
from ml_utils import get_direct_repeat_length
    

def generate_potential_candidates(df: pd.DataFrame,
                                  strain: str,
                                  segment: str
                                  )-> None:
    '''
        Generates potentially good candidates with a combinatorial approach.
        Only candidates with A at position 5 and direct repeat length 5 are
        allowed.
        :param df: data frame to calcualte the quantile from
        :param strain: Name of the strain
        :param segment: Name of the segment
    '''
    # select rows by segment
    df = df[df["Segment"] == segment]

    # define range by packaging signal and quantile of a dataset
    with open(os.path.join(DATAPATH, "Pelz2021", "packaging_signal.json"), "r") as f:
        packaging_signals = json.load(f)
        signals = packaging_signals[segment]
    
    seq_len = get_seq_len(strain, segment)
    s1 = signals["incorporation_start"] + signals["bundling_start"]
    e2 = seq_len - (signals["incorporation_end"] + signals["bundling_end"])
    s2 = int(df.Start.quantile(0.7))
    e1 = int(df.End.quantile(0.3))
    starts = np.arange(s1, s2)
    ends = np.arange(e1, e2)

    # remove all points (start and end) that do not have a A at position 5
    seq = get_sequence(strain, segment)
    starts = np.asarray([p for p in starts if seq[p] == "A"]) + 1
    ends = np.asarray([p for p in ends if seq[p] == "A"]) + 1

    # generate all combinations of start and end positions as df
    m = np.array(np.meshgrid(starts, ends)).T.reshape(-1, 2)
    strain_col = np.full(shape=int(m.size/2), fill_value=strain)
    seg_col = np.full(shape=int(m.size/2), fill_value=segment)
    m = np.c_[m, strain_col, seg_col]
    df = pd.DataFrame(m, columns=["Start", "End", "Strain", "Segment"])
    df.Start = df.Start.astype(int)
    df.End = df.End.astype(int)

    # calculate direct repeat length for each row and sort by them
    df["dir_rep"] = df.apply(get_direct_repeat_length, axis=1)
    df = df.sort_values(by=["dir_rep"], ascending=False)

    # write results to .csv file
    path = os.path.join(DATAPATH, "ML", f"rule_generated_{segment}_{strain}.csv")
    df.to_csv(path, index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Creates artificial DI RNA candidates by applying rules.")
    p.add_argument("--strain", "-s", type=str, help="Write name of strain to use")
    p.add_argument("--segment", "-g", type=str, help="Write name of segment to use")
    args = p.parse_args()

    if args.strain == "PR8":
        df = load_alnaji_2021()[args.strain]
    else:
        exit(f"No defined dataset for {args.strain}.")

    generate_potential_candidates(df, args.strain, args.segment)

