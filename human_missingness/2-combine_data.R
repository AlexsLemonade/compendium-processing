# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Purpose: Combine gene lists and sample counts from array and rna-seq platform
# data and make graphs and tables
# magrittr pipe
`%>%` <- dplyr::`%>%`
#-----------------------------Get genes per platform---------------------------#
# Read in the lists from prior
load("genes.per.array.RData")
load("rna.seq.genes.RData")

# Combine lists
genes.per.platform <- c(genes.per.affy, genes.per.illum, list(as.character(rna.seq.genes)))

# Get a list of all genes covered by over all the platforms
all.genes <- unique(unlist(genes.per.platform))

# Get a list of the number of genes covered by each platform
genes.covered <- lapply(genes.per.platform, function(x) {
  length(!is.na(match(x, all.genes)))})

# Make the list into a ratio
perc.genes.per.platform <- unlist(genes.covered)/length(all.genes)

# Keep the names of the platforms
names(perc.genes.per.platform) <- c(gsub("\\.tsv.gz", "", id.files),
                       gsub("\\.db", "", packages),
                       "RNA-seq")

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
platform.samples$from_acc <- recode(platform.samples$from_acc, !!!acc.convert)

# Find out how many samples per platform
samples.per.array <- summary(factor(platform.samples$from_acc))

# Tack on RNA-seq sample count at the end of the data
samples.per.platform <- c(samples.per.array, 133776)
names(samples.per.platform)[length(samples.per.platform)] <- "RNA-seq"

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

# Rectify the illlumina platforms mismatch. This is very janky. There is probably 
# a more elegant way to do it. 
key.indices <- grep("illuminaHuman", platform.table$platforms)
platform.table$platforms[key.indices] <- paste0(platform.table$platforms[key.indices], ".0")
platform.table$platforms <- gsub("illuminaHuman", "", platform.table$platforms)

for (x in key.indices) {
       tmp <- grep(platform.table$platforms[x], platform.table$platforms[-x], ignore.case = TRUE)
       platform.table[tmp, 3:4] <- platform.table[x, 3:4]
}
# Get rid of the version rows
platform.table <- platform.table[-key.indices, ]

# Write this info to a master table
write.csv(platform.table, file = "results/human_master_platform_info.csv", quote = FALSE)
