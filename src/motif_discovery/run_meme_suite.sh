#!/usr/bin/env bash

dir="../../data/meme_suite/alnaji2019"
dir_full="$dir/full_sequences"
dir_crop="$dir/cropped_sequences"

# Parsing arguments
# -d dataset -> name of dataset without '.fasta' (e.g. 'all', 'M', 'PB2')
# -c cropped -> defines format of the input sequences
    # full sequences where used: 'full'
    # sequences without deletion site: 'crop' 
    # only a window of size ##  around the deletion site was used: 'window_##'
# -m mode    -> mode that was used to create the sequences
    # only needs to be specified if 'window_##' is used
# -p program -> program to run (meme, streme, xstreme, glam)

while getopts d:c:m:p: flag
do
    case "${flag}" in
        d) dataset=${OPTARG};; # name of dataset (without .fasta)
        c) cropped=${OPTARG};; # 'crop', 'full' or 'window_##'
        m) mode=${OPTARG};; # mode to use; needs to be specified if 'window_##' is used
        p) program=${OPTARG};; # program to run (meme, streme, xstreme, glam)
    esac
done

if ! ([ "$mode" == "" ] || [ "$mode" == "combined" ] || [ "$mode" == "onlyremain" ])
then
    echo "when -m is selected 'combined' or 'onlyremain' has to be written after"
    exit
fi

if [ $cropped == "full" ]
then
    path="$dir_full/$dataset"
elif [ $cropped == "crop" ]
then
    path="$dir_crop/$dataset"
elif [ ${cropped:0:6} == "window" ]
then
    path="$dir/$cropped""_sequences_$mode/$dataset"
elif [ $cropped == "high_ngs" ]
then
    path="$dir/$cropped/$dataset"
else
    echo "use 'full', 'crop' or 'window_{n}' as argument for when calling script."
    echo "Input for -c/--cropped was $cropped"
    echo "Exiting script!"
    exit
fi

echo $path

if [ $program == "meme" ]
then
    meme "$path.fasta" -o "$path""_meme" -rna -nmotifs 20 -minw 8 -maxw 15
elif [ $program == "streme" ]
then
    streme --p "$path.fasta" -o "$path""_streme" -rna -nmotifs 20
elif [ $program == "xstreme" ] && [ $cropped == "high_ngs" ]
then
    xstreme --p "$path.fasta" --n "$path""_control.fasta" -o "$path""_xstreme" -rna
elif [ $program == "xstreme" ]
then
    xstreme --p "$path.fasta" -o "$path""_xstreme" -rna
elif [ $program == "glam" ]
then
    glam2 n "$path.fasta" -o "$path""_glam"
fi

firefox "$path""_$program/$program.html"

