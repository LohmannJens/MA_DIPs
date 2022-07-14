'''
    Estimates the delta G of the sequences using a sliding window approach and
    Vienna RNA. Step size and window size can be changed. Results are saved as
    CSV files.
'''
import os
import sys
import RNA

import pandas as pd

sys.path.insert(0, "..")
from utils import DATAPATH, RESULTSPATH, SEGMENTS
from utils import get_sequence


def delta_G_calculation(seq: str, p: int, w_s: int)-> float:
    '''
        Estimates the delta G of a part of a given sequence. The part is given
        by a middle point p and a window size (w_s).
        :param seq: full RNA sequence to calculate the delta G after cropping
        :param p: start point of the window for the sequence
        :param w_s: window size to calculate delta G

        :return: delta G of the cropped sequence
    '''
    c_seq = seq[p-int(w_s/2): p+int(w_s/2)]
    _, mfe = RNA.fold(c_seq)
    return mfe


def sliding_window_approach(w_s: int, step: int, seq: str)-> object:
    '''
        Function for the main part of the sliding window approach. Loops over 
        The whole sequence and calculates the delta G for each resulting
        subsequence.
        :param w_s: window size to calculate delta G
        :param step: step size for sliding the window
        :param seq: full RNA sequence of the segment

        :return: Data Frame with position on sequence and corresponding delta G
    '''
    results_dict = dict({"position": list(), "delta_G": list()})
    start = int(w_s/2)
    end = len(seq)-int(w_s/2)
    for i in range(start, end, step):
        delta_G = delta_G_calculation(seq, i, w_s)
        results_dict["position"].append(i)
        results_dict["delta_G"].append(delta_G)

    # check for last possible position in sequence
    if i < end:
        delta_G = delta_G_calculation(seq, end, w_s)
        results_dict["position"].append(end)
        results_dict["delta_G"].append(delta_G)
    
    df = pd.DataFrame(results_dict)
    return df


if __name__ == "__main__":
    window_size = 100
    step_size = 1
    for strain in ["Cal07", "NC", "Perth", "B_LEE"]:
        for s in SEGMENTS:
            sequence = get_sequence(strain, s)
            df = sliding_window_approach(window_size, step_size, sequence)
            path = os.path.join(DATAPATH, "energy_calculation", "sliding_window")
            file_path = os.path.join(path, f"{strain}_{s}_{window_size}_{step_size}.csv")
            df.to_csv(file_path)

