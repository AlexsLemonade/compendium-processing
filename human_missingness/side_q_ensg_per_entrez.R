# C. Savonen
# 
# Purpose: how many genes have multiple ensg's per their entrez ids. 

# Directory set up
results.dir <- "results"
plots.dir <- "results/plots"

# Read in the percent zeroes per gene table 
rna.seq.perc.zeroes <- read.table("data/ARCHS4_human_matrix_percent_zeroes.txt", sep = "\t",
                                  skip = 1, stringsAsFactors = FALSE)

# Get rid of rows without a Entrez ID
rna.seq.perc.zeroes <- rna.seq.perc.zeroes[!is.na(rna.seq.perc.zeroes$V1),]

# 
get.all <- mapIds(org.Hs.eg.db, 
                  keys = as.character(rna.seq.perc.zeroes$V1),
                  column = "ENSEMBL", keytype = "ENTREZID", 
                  multiVals = "list")

# Get the number of ensg IDs per entrez ids
genes.per.mapping <- vapply(get.all, length, FUN.VALUE = integer(1))

jpeg(file.path(plots.dir, "distribution_genes_per_entrez.jpeg"))
# Plot the distribution of this
plot(density(genes.per.mapping), main = "Distribution of Number of ENSG IDs")
dev.off()

# Get numbers on how many genes of each of these
x <- summary(factor(genes.per.mapping))
jpeg(file.path(plots.dir, "genes_per_entrez.jpeg"))
plot(names(x), x, ylab = "How many genes with this number of ensg?",
                  xlab = "How many ensg's assoc with one entrez id?")
dev.off()

# How many genes have more than one ENSG ID associated? 
sum(summary(factor(genes.per.mapping))[-1])