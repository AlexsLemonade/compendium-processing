#!/bin/bash

# Make some directories
mkdir data
mkdir results
mkdir data/id_refinery

# Download id_refinery files if they aren't there
FILE=$data/id_refinery/all_1536267482.zip

if [ ! -f $FILE ]; then
wget https://zenodo.org/record/1410647/files/all_1536267482.zip \
-P data/id_refinery
unzip data/id_refinery/all_1536267482.zip -d data/id_refinery/
fi

# Download supported microarray list
wget https://raw.githubusercontent.com/AlexsLemonade/refinebio/33c5275f6cfe5caa124ecff8eb83848d779fcdde/config/supported_microarray_platforms.csv -P data

# Download human microarray list from arrayexpress
wget -O data/array_express_human.txt "https://www.ebi.ac.uk/arrayexpress/ArrayExpress-Experiments.txt?keywords=human+&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=&directsub=on"

# Run the scripts!
Rscript 0-get_human_platform_info.R
Rscript 1a-human_array_platforms_coverage.R
Rscript 1b-human_rnaseq_gene_coverage.R
Rscript 2-combine_data.R
