# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Purpose: Combine gene lists and sample counts from array and rna-seq platform
# data and make graphs and tables
# magrittr pipe
`%>%` <- dplyr::`%>%`

#-----------------------------Get genes per platform---------------------------#
# Read in human platform list
platforms <- read.csv("data/exp_acc_human_only.csv", stringsAsFactors = FALSE)

# Read in the lists from prior
load("genes.per.array.RData")
load("rna.seq.genes.RData")

# Combine lists
genes.per.platform <- c(genes.per.affy, genes.per.illum, list(as.character(rna.seq.perc.zeroes$ensembl)))
# Keep the names of the platforms
names(genes.per.platform ) <- c(gsub("\\.tsv.gz", "", id.files),
                                    names(genes.per.illum), "rnaseq")

# Get a list of all genes covered by over all the platforms
all.genes <- as.character(unique(unlist(genes.per.platform)))
all.genes <- grep("ENS", all.genes, value = TRUE)

# Get a list of the number of genes covered by each platform
genes.covered <- lapply(genes.per.platform, function(x) {
  length(!is.na(match(x, all.genes)))})

# Make the list into a ratio
perc.genes.per.platform <- unlist(genes.covered)/length(all.genes)

# Keep the names of the platforms
names(perc.genes.per.platform) <- c(gsub("\\.tsv.gz", "", id.files),
                       names(genes.per.illum),"rnaseq")

# Print out the stats as a barplot
jpeg("results/Percent_genes_covered_per_platform.jpeg")
par(mar = c(8,4,1,1))
barplot(perc.genes.per.platform *100, las = 2, ylab = "Percent of all array genes covered")
dev.off()

# Write this to a table in a csv file
write.csv(perc.genes.per.platform , file = "results/genes_detected_per_platform.csv", quote = FALSE)

#------------------Get number of samples per platform--------------------------#
# Load in sample info
load("GEO_exp_info.RData")

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
samples.per.platform <- c(samples.per.array, 133776)
names(samples.per.platform)[length(samples.per.platform)] <- "rnaseq"

# Print out the stats as a barplot
jpeg("results/Number_of_samples_per_platform.jpeg")
par(mar = c(12,6,1,1))
barplot(samples.per.array, las = 2)
title(ylab="Number of samples per platform", mgp = c(5,5,10))
dev.off()

# Write this to a table in a csv file
write.csv(samples.per.platform, file = "results/samples_per_platform.csv", quote = FALSE)

#------------------------ Make everything into one table-----------------------#
samples.per.platform <- data.frame(samples.per.platform, 
                                   'platforms' = names(samples.per.platform))
perc.genes.per.platform  <- data.frame('num_genes' = unlist(genes.covered),
                                       'perc_all_genes' = perc.genes.per.platform , 
                                       'platforms' = names(perc.genes.per.platform))

platform.table <- samples.per.platform %>% dplyr::full_join(perc.genes.per.platform, 
                                                        by = 'platforms')

# Write this info to a master table
write.csv(platform.table, file = "results/human_master_platform_info.csv", quote = FALSE)



#------------------- Make gene by percent of samples table---------------------#
tmp <- match(names(genes.per.platform), samples.per.platform$platforms)
samples.per.platform <- samples.per.platform[tmp, ]

# Make a genes vs how many samples have them matrix
mat <- data.frame("genes" = all.genes,
                  "tot.samples" = 0, 
                  stringsAsFactors = FALSE)

# Run through for each platform
for (x in 1:length(genes.per.platform))  {
  
  # Get number of samples in platform
  num.samples <- samples.per.platform[x,1]
  
  # For RNA-seq we have to calculate it differently
  if (names(genes.per.platform)[x] == "rnaseq") {
    
    # Put the perc.zeroes for rnaseq in the same order as the master matrix
    tmp <- match(rna.seq.perc.zeroes$ensembl, mat$genes)
    tmp <- tmp[!is.na(tmp)]
    
    # Add the number of samples that contain the gene in rna-seq to the running total
    mat$tot.samples[tmp] <- mat$tot.samples[tmp] + num.samples*(1-rna.seq.perc.zeroes$perc.zeroes)
  } else {
    # Match the genes to the master matrix
    tmp <- match(genes.per.platform[[x]], mat$genes)
    tmp <- tmp[!is.na(tmp)]
    
    # Add the samples to the respective genes' running totals. 
    mat$tot.samples[tmp] <- mat$tot.samples[tmp] + num.samples
  }
}

# Calculate the percentages of samples that can detect each gene
total.samples <- sum(samples.per.platform$samples.per.platform)
mat$perc.samples <- mat$tot.samples/total.samples

# Print the distrbution of these percentages
jpeg("results/detection_percentage_distribution.jpeg")
plot(density(mat$perc.samples), xlab = "Ratio of samples which have the gene", main = "Distribution of Detection Percentages for All Genes")
dev.off()

# Make it a histogram
jpeg("results/detection_percentage_histogram.jpeg")
hist(mat$perc.samples, xlab = "Ratio of samples which have the gene", main = "Histogram of Detection Percentages for All Genes")
dev.off()

save(mat, file = "perc_samples_per_genes.RData")
