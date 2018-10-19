#!/bin/bash

declare -a seed_arr=("1663" "3250" "4563" "5089" "7210" "7272" "7609" "8757" "8859" "9889")
declare -a missing_arr=("MCAR" "MRR")
declare -a processing_arr=("all_log" "all_un" "microarray" "seq_log" "seq_un")

mkdir imputed/MCAR -p
mkdir imputed/MRR -p

for seed in "${seed_arr[@]}"
do
  for miss in "${missing_arr[@]}"
  do
    for process in "${processing_arr[@]}"
    do
      input_file="data/masked/${miss}/${process}_${seed}_${miss}.pcl"
      echo "Processing ${input_file}"
      output_file="imputed/${miss}/${process}_${seed}_${miss}_knnimpute.pcl"
      KNNImputer -i $input_file -o $output_file -s 0 -l 20
    done
  done
done
