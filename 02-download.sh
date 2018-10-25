#!/bin/bash

mkdir data

# wget all the zebrafish RNA-seq data
mkdir data/rnaseq
wget -i selected_RNAseq_samples.txt -P data/rnaseq

# random selection of zebrafis microarray samples
mkdir data/microarray
wget -i selected_GEO_samples.txt -P data/microarray

# quantile normalization targets for zebrafish
mkdir data/qn_target 
wget https://s3.amazonaws.com/data-refinery-s3-circleci-prod/d8gyazczr661umla1979zxu7_1538153041_target.tsv \
  -P data/qn_target
