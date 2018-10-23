# C. Savonen 
# CCDL for ALSF 
# This script is adapted from the example script on ARCHS4
# 
# Objective: get an idea of the coverage of each gene in human
# 
#download.file("https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5", destfile = "archs4.h5")

# R script to download selected samples
# Copy code and run on a local machine to initiate download
# Check for dependencies and install if missing
packages <- c("rhdf5")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  print("Install required packages")
  source("https://bioconductor.org/biocLite.R")
  biocLite("rhdf5")
}
library("rhdf5")

destination_file = "human_matrix.h5"
extracted_expression_file = "example_expression_matrix.tsv"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
  print("Downloading compressed gene expression matrix.")
  url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
  download.file(url, destination_file, quiet = FALSE)
} else{
  print("Local file already exists.")
}

# Selected samples to be extracted
platforms = read.csv("exp_acc_human_only.csv")
platforms$external_accession

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/Sample_geo_accession")
# Identify columns to be extracted
sample_locations = which(samples %in% samp)

tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
genes = h5read(destination_file, "meta/genes")
series = h5read(destination_file, "meta/Sample_series_id")

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]

# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))