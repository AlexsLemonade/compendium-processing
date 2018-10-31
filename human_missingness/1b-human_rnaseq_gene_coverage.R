# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Purpose: Obtain list of genes measured reliably on rna-seq platforms
# This script is generally adapted from the example script on ARCHS4
# 
# This script is intended to be ran from using the bash script:
# `human_missingness/run_scripts.sh`

# Objective: get an idea of the coverage of each gene in human

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Attach libraries
library(org.Hs.eg.db)
library(rhdf5)
#------------------Set Up Functions------------------#

# Directory set up
results.dir <- "results"
plots.dir <- file.path("results", "plots")

# human file directory path
human.file <- file.path("data", "human_matrix.h5")

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if (!file.exists(human.file)) {
  download.file("https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5",
                human.file, quiet = FALSE)
  } else {
  print("Local file already exists.")
}

# Retrieve information from compressed data
metadata <- data.frame(platform = h5read(human.file,
                                         "meta/Sample_platform_id"),
                        tissue = h5read(human.file,
                                        "meta/Sample_source_name_ch1"),
                        sample.id = h5read(human.file,
                                           "meta/Sample_geo_accession"),
                        series = h5read(human.file,
                                        "meta/Sample_series_id"))

# Store how many rna.seq samples there are
n.rna.seq.samples <- nrow(metadata)

# Get a list of the genes
genes <- h5read(human.file, "meta/genes")

#----Get RNA-Seq gene list for a subset of samples - Due to RAM constraints----#

# Make a random subset of the samples
samples <- runif(200, min = 1, max = nrow(metadata))

# extract gene expression from compressed data
expression <- h5read(human.file, "data/expression", index = list(1:length(genes),
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
if (!file.exists(file.path("data", "ARCHS4_human_matrix_percent_zeroes.txt"))) {
# Jackie used this code in order to get the txt file we read in:  
  # Establish file name
  human_file <- file.path("data", "human_matrix.h5")
  
  # Read in all the expression data 
  gene.mat <- rhdf5::h5read(human.file, "/data/expression")
  
  # Get the number of samples of each row that are 0's
  gene.zero.counts <- apply(gene.mat, 1, function(x) sum(x == 0) / length(x))
  
  # Retrieve the entrez.ids
  gene.ids <- rhdf5::h5read(human.file, "/meta/gene_entrezid")
  
  # Make a data frame with the percent zeroes and entrez ids
  genes.df <- data.frame(entrez.id = gene.ids, percent.zero.counts = gene.zero.counts)
  
  # Write the percent of zeroes per gene to a text file
  readr::write_tsv(genes.df, "ARCHS4_human_matrix_percent_zeroes.txt")
}

# Read in the percent zeroes per gene table 
rna.seq.perc.zeroes <- read.table(
                        file.path("data","ARCHS4_human_matrix_percent_zeroes.txt"),
                                  sep = "\t", skip = 1, stringsAsFactors = FALSE)

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
saveRDS(rna.seq.perc.zeroes, file = file.path(results.dir, "rna_seq_genes.RDS"))
saveRDS(n.rna.seq.samples, file = file.path(results.dir, "n_rna_seq_samples.RDS"))
