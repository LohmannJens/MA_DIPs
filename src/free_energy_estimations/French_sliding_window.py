'''
    Influenza virus RNA polymerase simulation script that searches for t-loops
    and checks up/down stream for intermol bp.
    run_tloop_analysis() inspired by AJ te Velthuis, Sept 2020 from French
    et al. 2021

    Used their main function and adapted it to automize it for other datasets
'''
import os
import sys
import RNA

sys.path.insert(0, "..")
from utils import  DATAPATH, SEGMENTS
from utils import get_sequence


def run_tloop_analysis(rna: str,
                       filename: str
                       )-> None:
    '''
        Loops over a given RNA sequence and calculates the delta G for a
        possible t loop at each position. The output can also be the delta
        delta G of the t loop (compared to a sequence, that can also create
        internal loops.
        :param rna: rna sequence
        :param filename: file to save the calculated values in

        :return: None
    '''
    Footprint = 20 # nts bound by polymerase
    TloopDuplex = 48 # Size of a t-loop
    Duplex = int(TloopDuplex / 2) # Length of one side of a t-loop

    NP = 24
    Uloop = "&"

    path = os.path.join(DATAPATH, "energy_calculation", "sliding_window", filename)
    if os.path.exists(path):
        os.remove(path)

    f = open(path, "w+")
    f.write("position,delta_G\n")

    # invert input sequence to start at 3' end
    NegRNA = rna[::-1]
    Length = len(NegRNA)
    Bubble = Footprint + Duplex
    End = int((Length - Footprint + 1))

    for i in range(1, End-1):
        # This is for the Tloop Duplex
        if i <= Duplex:
            Prime3 = NegRNA[0:Duplex]
        else:
            Prime3 = NegRNA[i-Duplex:i]
        Prime3inv = Prime3[::-1]

        if i <= End:
            Prime5 = NegRNA[Footprint+i:Bubble+i]
        else:
            Prime5 = NegRNA[Footprint+i:Length]
        Prime5inv = Prime5[::-1]

        duplex = RNA.duplexfold(Prime5inv, Prime3inv)

        # This is for the t loops that can occure before and after the polymerase
        if i <= NP:
            Upstream = i
            Downstream = 0
        else:
            Upstream = NP
            Downstream = i-NP
     
        Ahead = NegRNA[Footprint+i:Footprint+NP+i]
        Aheadinv = Ahead[::-1]
        Down = NegRNA[Downstream:i]
        Downinv = Down[::-1]
        Other = Aheadinv + Uloop + Downinv
        Other = Aheadinv

        (ss, mfe_dimer) = RNA.cofold(Other)

        result = duplex.energy + mfe_dimer

        f.write(f"{i},{result}\n")

    f.close()
    print (f"Done {filename}")


if __name__ == "__main__":
    strains = ["Cal07", "NC", "Perth", "BLEE"]
    for strain in strains:
        for s in SEGMENTS:
            rna = get_sequence(strain, s)
            filename = f"{strain}_{s}_1_1.csv"
            run_tloop_analysis(rna, filename)

