# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Purpose: Combine gene lists and sample counts from array and rna-seq platform
# data and make graphs and tables
# 
# This script is intended to be ran from using the bash script:
# `human_missingness/run_scripts.sh`
# 
# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Directory set up
results.dir <- "results"
plots.dir <- file.path("results", "plots")

#-----------------------------Get genes per platform---------------------------#
# Read in human platform list
platforms <- read.csv(file.path("data", "exp_acc_human_only.csv"), stringsAsFactors = FALSE)

# Read in the lists from previously
genes.per.illum <- readRDS(file.path(results.dir, "genes_per_illumina_array.RDS"))
genes.per.affy <- readRDS(file.path(results.dir,"genes_per_affy_array.RDS"))
rna.seq.perc.zeroes <- readRDS(file.path(results.dir,"rna_seq_genes.RDS"))
n.rna.seq.samples <- readRDS(file.path(results.dir,"n_rna_seq_samples.RDS"))

# Combine lists
genes.per.platform <- c(genes.per.affy, genes.per.illum)
genes.per.platform[["rnaseq"]] <- unlist(as.character(rna.seq.perc.zeroes$ensembl))

# Get a list of all genes covered by over all the platforms
all.genes <- as.character(unique(unlist(genes.per.platform)))
all.genes <- grep("ENS", all.genes, value = TRUE)

# Get a list of the number of genes covered by each platform
genes.covered <- lapply(genes.per.platform,
                        function(x) length(!is.na(match(x, all.genes))))

# Make the list into a ratio
perc.genes.per.platform <- unlist(genes.covered)/length(all.genes)

# Keep the names of the platforms
names(perc.genes.per.platform) <- names(genes.per.platform)

# Print out the stats as a barplot
jpeg(file.path(plots.dir, "percent_genes_covered_per_platform.jpeg"))
par(mar = c(8,4,1,1))
barplot(perc.genes.per.platform *100, las = 2, ylab = "Percent of all array genes covered")
dev.off()

# Write this to a table in a csv file
write.csv(data.frame("percent_all_genes_detected" = perc.genes.per.platform),
          file = file.path(results.dir, "genes_detected_per_platform.csv"),
          quote = FALSE)

#------------------Get number of samples per platform--------------------------#
# Load in sample info
hs.gpl <- readRDS(file.path(results.dir,"GEO_exp_info.RDS"))

# Make a recode key from the external accessions to internal accession
acc.convert <- as.character(platforms$internal_accession)
names(acc.convert) <- platforms$external_accession

# Only keep samples from the platforms refine.bio supports
platform.samples <- hs.gpl$gsm %>% 
  dplyr::filter(from_acc %in% names(acc.convert))

# Recode into our internal accession codes
platform.samples$from_acc <- dplyr::recode(platform.samples$from_acc, !!!acc.convert)

# Find out how many samples per platform
samples.per.array <- summary(factor(platform.samples$from_acc))

# Tack on RNA-seq sample count at the end of the data
samples.per.platform <- c(samples.per.array, n.rna.seq.samples)
names(samples.per.platform)[length(samples.per.platform)] <- "rnaseq"

# Print out the stats as a barplot
jpeg(file.path(plots.dir,"number_of_samples_per_platform.jpeg"))
par(mar = c(12,6,1,1))
barplot(samples.per.array, las = 2)
title(ylab="Number of samples per platform", mgp = c(5,5,10))
dev.off()

write.csv(data.frame("samples_per_platform" = samples.per.platform),
          file = file.path(results.dir, "samples_per_platform.csv"),
          quote = FALSE)

#------------------------ Make everything into one table-----------------------#
samples.per.platform <- data.frame(samples.per.platform, 
                                   'platforms' = names(samples.per.platform))
perc.genes.per.platform  <- data.frame('num_genes' = unlist(genes.covered),
                                       'perc_all_genes' = perc.genes.per.platform , 
                                       'platforms' = names(perc.genes.per.platform))

platform.table <- samples.per.platform %>% dplyr::full_join(perc.genes.per.platform, 
                                                        by = 'platforms')

# Write this info to a master table
write.csv(platform.table, 
          file = file.path(results.dir,"human_master_platform_info.csv"),
          quote = FALSE)

#------------------- Make gene by percent of samples table---------------------#
tmp <- match(names(genes.per.platform), samples.per.platform$platforms)
samples.per.platform <- samples.per.platform[tmp, 1]

# Make a master list with the gene lists and sample numbers
# Only keep ensembl gene ids "ENS"
platform.list <- lapply(genes.per.platform, function (x) {
                      x.ind <- parent.frame()$i[]
                         list("genes" = grep("ENS", x, value = TRUE) , 
                              "n.samples" = samples.per.platform[x.ind])
                      })
# Add extra thing for rna-seq
platform.list$rnaseq$perc.detected <- (1 - rna.seq.perc.zeroes$perc.zeroes)

# Make a samples per gene matrix
samples.per.gene <- data.frame("gene" = all.genes,
                  "tot.samples" = 0,
                  stringsAsFactors = FALSE)

# Make a loop to calculate the number of samples that have measurements for
# each gene
for (platform_ind in 1:length(platform.list)) {
    # Set the platform and type variables by this index
    platform <- platform.list[[platform_ind]]
    type <- names(platform.list)[platform_ind]
      
    # Match genes of platform to genes of samples.per.gene matrix
    gene.indices <- match(platform$genes, all.genes) 
          
    if (type == "rnaseq"){ # if the type of data is rnaseq,
      add <- platform$n.samples*platform$perc.detected # multiply total samples by percent
    } else {
      add <- platform$n.samples # if NOT rnaseq, use sample total as is.
    }
    # Add number of samples to running total
    samples.per.gene$tot.samples[gene.indices] <- 
    samples.per.gene$tot.samples[gene.indices] + add
}

# Calculate the percentages of samples that can detect each gene
total.samples <- sum(samples.per.platform)
samples.per.gene$perc.samples <- samples.per.gene$tot.samples/total.samples

# Print the distrbution of these percentages
jpeg(file = file.path(plots.dir, "detection_percentage_distribution.jpeg"))
plot(density(samples.per.gene$perc.samples), xlab = "Ratio of samples which have the gene",
     main = "Distribution of Detection Percentages for All Genes")
dev.off()

# Make it a histogram
jpeg(file.path(plots.dir, "detection_percentage_histogram.jpeg"))
hist(samples.per.gene$perc.samples, xlab = "Ratio of samples which have the gene",
     main = "Histogram of Detection Percentages for All Genes")
dev.off()

# Save as an RData object
write.csv(samples.per.gene, file = file.path(results.dir, "perc_samples_per_genes.csv"),
          quote = FALSE, row.names = FALSE)
