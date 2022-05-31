#!/usr/bin/env bash

dir=../../data/energy_calculation
dir_full=$dir/full_sequences
dir_cropped=$dir/cropped_sequences


for f in $dir_cropped/*.fasta;
do
    out=${f%.fasta}.fold
    echo $out
    RNAfold --infile $f --outfile=temp.fold
    mv temp.fold $out
done

# remove temp .ps files
rm *.ps

