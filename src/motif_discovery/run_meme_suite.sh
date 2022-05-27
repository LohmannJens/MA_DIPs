#!/usr/bin/env bash

dir="../../data/meme_suite/alnaji2019"
dir_full="$dir/full_sequences"
dir_crop="$dir/cropped_sequences"

#dataset="IVA"

while getopts d:c:p: flag
do
    case "${flag}" in
        d) dataset=${OPTARG};; # name of dataset (without .fasta)
        c) cropped=${OPTARG};; # 'crop', 'full' or 'window'
        p) program=${OPTARG};; # program to run (meme, streme, xstreme, glam)
    esac
done


if [[ $cropped == "full" ]]
then
    path="$dir_full/$dataset"
elif [[ $cropped == "crop" ]]
then
    path="$dir_crop/$dataset"
elif [[ ${cropped:0:6} == "window" ]]
then
    path="$dir/$cropped""_sequences/$dataset"
else
    echo "use 'full', 'crop' or 'window_{n}' as argument for when calling script."
    echo "Input for -c/--cropped was $cropped"
    echo "Exiting script!"
    exit
fi

echo $path

if [[ $program == "meme" ]]
then
    meme "$path.fasta" -o "$path""_meme" -rna -nmotifs 20 -minw 8 -maxw 15
elif [[ $program == "streme" ]]
then
    streme --p "$path.fasta" -o "$path""_streme" -rna -nmotifs 20
elif [[ $program == "xstreme" ]]
then
    xstreme --p "$path.fasta" -o "$path""_xstreme" -rna
elif [[ $program == "glam" ]]
then
    glam2 n "$path.fasta" -o "$path""_glam"

fi

firefox "$path""_$program/$program.html"
