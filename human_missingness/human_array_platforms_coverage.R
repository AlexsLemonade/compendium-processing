# C. Savonen 
# CCDL for ALSF 
#
# Identify which supported microarray platforms are Homo sapiens.
#  
#------------------- Get the supported platforms file from GitHub--------------#
download.file(url = "https://raw.githubusercontent.com/AlexsLemonade/refinebio/33c5275f6cfe5caa124ecff8eb83848d779fcdde/config/supported_microarray_platforms.csv",
              destfile = "supported_platforms.csv")
# Read in the file:
platforms <- read.csv("supported_platforms.csv")

# Need dplyr
if (!("dplyr" %in% installed.packages())) {
  install.packages("dplyr")
}
library(dplyr)
#------------------------------ Get info on GEO data: ------------------------#
# Install and attach this package:
if (!("GEOmetadb" %in% installed.packages())) {
source("https://bioconductor.org/biocLite.R")
biocLite("GEOmetadb")
}
library(GEOmetadb)

# Check if the database is downloaded
if(!file.exists('GEOmetadb.sqlite')){
  getSQLiteFile()
}

# If it has been, connect to it
if(file.exists('GEOmetadb.sqlite')){
  con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
}
# Get homo sapiens gpls
hs.gpl <- dbGetQuery(con, paste0('SELECT gpl FROM gpl WHERE gpl.organism LIKE "%homo%"'))

# Get info about these experiments
hs.gpl <- geoConvert(hs.gpl$gpl, out_type = c("gse", "gpl"),
                            sqlite_db_name = "GEOmetadb.sqlite")

# url for getting the array express experiment info
array.express <- "https://www.ebi.ac.uk/arrayexpress/ArrayExpress-Experiments.txt?keywords=human+&organism=Homo+sapiens&exptype%5B%5D=&exptype%5B%5D=&array=&directsub=on"

# Download human array info array
download.file(array.express, destfile = "array_express_human.txt")

# Read in the array express human experiments meta info
hs.array <- read.table("array_express_human.txt", sep = "\t", skip = 1, stringsAsFactors = FALSE)

#-----------------------Only keep platforms that are human---------------------#
# List of both GEO and array express accessions
hs.acc <- c(unique(hs.gpl$gse$from_acc), hs.array$V1)

# Make dataframe with only human acc
platforms <- platforms %>% dplyr::filter(external_accession %in% hs.acc)

# Write info to a csv file
write.csv(platforms, file = "exp_acc_human_only.csv", quote = FALSE,
          row.names = FALSE)

#---------------Identify genes that are measured per platform------------------#
# From Identifer Refinery repo, get the gene lists
# This is ~1.4GB so it may take a bit
if (!file.exists("id_refinery")){
  # Create a directory
  dir.create("id_refinery")
  
  # Download annotation file from identifier_refinery repository
  download.file("https://zenodo.org/record/1410647/files/all_1536267482.zip", 
                destfile = "id_refinery/id_refinery.zip")
  # Unzip it
  system("unzip id_refinery.zip")
}
# Get list of id_refinery files
id.refinery <- dir("id_refinery")

# Only keep those for the platforms we care about (human in this case)
tmp <- match(unique(platforms$internal_accession), gsub(".tsv.gz", "", id.refinery))
id.files <- id.refinery[tmp[!is.na(tmp)]]

# Read in the genes from every file
setwd("id_refinery")
genes.per.platform <- lapply(id.files, function (x) {
  index_path <- paste('gzcat', x)
  suppressWarnings(index_data <- data.table::fread(
    cmd = index_path,
    stringsAsFactors = FALSE,
    sep = "\t",
    header=TRUE,
    autostart = 10,
    data.table = FALSE,
    check.names = FALSE,
    fill = TRUE,
    na.strings = "",
    showProgress = FALSE)
  )
  unique(index_data$ENSEMBL)
  })
# Get a list of all genes covered by over all the platforms
all.genes <- unique(unlist(genes.per.platform))

# Get a list of the number of genes covered by each platform
genes.covered <- lapply(genes.per.platform, function(x) {
  length(!is.na(match(x, all.genes)))})

# Make the list into a ratio
perc.genes <- unlist(genes.covered)/length(all.genes)
names(perc.genes) <- gsub("\\.tsv.gz", "", id.files)

# Print out the stats as a barplot
jpeg("../Percent_genes_covered_per_array.jpeg")
par(mar = c(8,4,1,1))
barplot(perc.genes*100, las = 2, ylab = "Percent of all array genes covered")
dev.off()
