#!/usr/bin/env bash

dir=../../data/energy_calculation
dir_full=$dir/full_sequences
dir_cropped=$dir/cropped_sequences
dir_shuffled=$dir/cropped_sequences_shuffled
dir_random=$dir/cropped_sequences_randomcrop

for f in $dir_random/*.fasta;
do
    out=${f%.fasta}.fold
    echo $out
    RNAfold --infile $f --outfile=temp.fold
    mv temp.fold $out
done

# remove temp .ps files
rm *.ps

