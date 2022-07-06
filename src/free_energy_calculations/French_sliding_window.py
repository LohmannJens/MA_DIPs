#Influenza virus RNA polymerase simulation script that searches for t-loops and checks up/down streamfor intermol bp
#run_tloop_analysis() inspired by AJ te Velthuis, Sept 2020 from French et al. 2021
#used their main function and adapted it to automize it for alnajis dataset

import os
import sys
import RNA

sys.path.insert(0, "..")
from utils import  DATAPATH, SEGMENTS
from utils import get_sequence


def run_tloop_analysis(rna: str, filename: str)-> None:
    '''
        Loops over a given RNA sequence and calculates the delta G for a
        possible t loop at each position. The output can also be the delta
        delta G of the t loop (compared to a sequence, that can also create
        internal loops.
        :param rna: rna sequence
        :param filename: file to save the calculated values in

        :return: None
    '''
    Footprint = 20 # polymerase properties
    NP = 24 # NP footprint
    TloopDuplex = 48
    Duplex = int(TloopDuplex / 2)
    Uloop = "&" #Use & for co-fold to compute long-range interactions between upstream and downstream sequences.
    Swindow = 1 #size of sliding window

    path = os.path.join(DATAPATH, "energy_calculation", "sliding_window")
    f = open(os.path.join(path, filename), "w+")
    f.write("position,delta_G\n")

    # invert input sequence to start at 3' end
    NegRNA = rna[::-1]
    Length = len(NegRNA)
    Bubble = Footprint + Duplex
    End = int((Length - Footprint + 1) / Swindow)

    for i in range(1, End-1, Swindow):
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

        if i <= Duplex:
            Prime3 = NegRNA[0:Upstream]
        else:
            Prime3 = NegRNA[i-Duplex:i]
        Prime3inv = Prime3[::-1]

        if i <= End:
            Prime5 = NegRNA[Footprint+i:Bubble+i]
        else:
            Prime5 = NegRNA[Footprint+i:Length]
        Prime5inv = Prime5[::-1]

        #use duplex fold to compute t-loop from Vienna package because it ignores intermol bp
        duplex = RNA.duplexfold(Prime5inv, Prime3inv)

        #use cofold from Vienna package to check for bp in sequence upstream and downstream of t-loop.
        #Various options can be checked separately, including just upstream seq, just downstream seq, or both seq
        Other = Aheadinv + Uloop + Downinv
        #Other = Aheadinv
        #Other = Downinv
        (ss, mfe_dimer) = RNA.cofold(Other)
        DDeltaG = duplex.energy - mfe_dimer

     #   f.write(f"{i},{duplex.energy}\n")
     #   f.write(f"{i},{mfe_dimer}\n")
        f.write(f"{i},{DDeltaG}\n")

    # close .txt file that deltaG values were written to
    f.close()
    print (f"Done {filename}")


if __name__ == "__main__":
    strains = ["Cal07", "NC", "Perth", "B_LEE"]
    for strain in strains:
        for s in SEGMENTS:
            rna = str(get_sequence(strain, s).seq)
            filename = f"{strain}_{s}_1_1.csv"
            run_tloop_analysis(rna, filename)

