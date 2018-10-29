# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Purpose: Obtain list of genes measured reliably on rna-seq platforms
# This script is generally adapted from the example script on ARCHS4
#
# This script is intended to be ran from using the bash script:
# `human_missingness/run_scripts.sh`
# 
# Objective: get an idea of the coverage of each gene in human
library(org.Hs.eg.db)
library(rhdf5)

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if (!file.exists("data/human_matrix.h5")) {
  download.file("https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5",
                "data/human_matrix.h5", quiet = FALSE)
  } else {
  print("Local file already exists.")
}

# Retrieve information from compressed data
metadata <- data.frame(platform = h5read("data/human_matrix.h5", "meta/Sample_platform_id"),
                        tissue = h5read("data/human_matrix.h5", "meta/Sample_source_name_ch1"),
                        sample.id = h5read("data/human_matrix.h5", "meta/Sample_geo_accession"),
                        series = h5read("data/human_matrix.h5", "meta/Sample_series_id"))

# Store how many rna.seq samples there are
n.rna.seq.samples <- nrow(metadata)

# Get a list of the genes
genes <- h5read("data/human_matrix.h5", "meta/genes")

#----Get RNA-Seq gene list for a subset of samples - Due to RAM constraints----#

# Make a random subset of the samples
samples <- runif(200, min = 1, max = nrow(metadata))

# extract gene expression from compressed data
expression <- h5read("data/human_matrix.h5", "data/expression", index = list(1:length(genes),
                                                                      samples))
H5close()
rownames(expression) <- genes
colnames(expression) <- metadata$sample.id[samples]

rna.seq.perc.zeroes <- apply(expression, 1, function (x) {
  (sum(x == 0)/ncol(expression))
})

#------------------ Get RNA-Seq gene list for all samples---------------------#
#
########## Don't run this unless you have to. Takes a lot of ram and time.###### 
if (!file.exists("data/ARCHS4_human_matrix_percent_zeroes.txt")) {
# Jackie this code in order to get the txt file we read in:  
  # Establish file name
  human_file <- "data/human_matrix.h5"
  
  # Read in all the expression data 
  gene_mat <- rhdf5::h5read(human_file, "/data/expression")
  
  # Get the number of samples of each row that are 0's
  gene_zero_counts <- apply(gene_mat, 1, function(x) sum(x == 0) / length(x))
  
  # Retrieve the entrez_ids
  gene_ids <- rhdf5::h5read(human_file, "/meta/gene_entrezid")
  
  # Make a data frame with the percent zeroes and entrez ids
  genes_df <- data.frame(entrez_id = gene_ids, percent_zero_counts = gene_zero_counts)
  
  # Write the percent of zeroes per gene to a text file
  readr::write_tsv(genes_df, "ARCHS4_human_matrix_percent_zeroes.txt")
}

# Read in the percent zeroes per gene table 
rna.seq.perc.zeroes <- read.table("data/ARCHS4_human_matrix_percent_zeroes.txt", sep = "\t",
                          skip = 1, stringsAsFactors = FALSE)

# Get rid of rows without a Entrez ID
rna.seq.perc.zeroes <- rna.seq.perc.zeroes[!is.na(rna.seq.perc.zeroes$V1),]

# Convert to a list
rna.seq.perc.zeroes <- data.frame('ensembl' = mapIds(org.Hs.eg.db, 
                                          keys = as.character(rna.seq.perc.zeroes$V1),
                                          column = "ENSEMBL", keytype = "ENTREZID"),
                        'perc.zeroes' = rna.seq.perc.zeroes$V2, 
                        stringsAsFactors = FALSE)

# Get rid of rows without an ensembl gene ID
rna.seq.perc.zeroes <- rna.seq.perc.zeroes[!is.na(rna.seq.perc.zeroes$ensembl), ]

# Save this info to an RData file
save(list = c("rna.seq.perc.zeroes", "n.rna.seq.samples"), file = "rna.seq.genes.RData")
