# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Purpose: Identify genes that are measured per array platform 
#
#----------------------------Get Affy Gene Lists-------------------------------#
# magrittr pipe
`%>%` <- dplyr::`%>%`

# Read in human platform list
platforms <- read.csv("data/exp_acc_human_only.csv", stringsAsFactors = FALSE)

# Get list of id_refinery files
id.refinery <- dir("data/id_refinery")

# Only keep those for the platforms we care about (human in this case)
tmp <- match(unique(platforms$internal_accession), gsub(".tsv.gz", "", id.refinery))
id.files <- id.refinery[tmp[!is.na(tmp)]]

# For each id refinery file that matches our list, print out the ensembl IDs
genes.per.affy <- lapply(id.files, function (x) {
  index_path <- paste0('gzcat data/id_refinery/', x)
  # Read the table
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
  # Print out the unique ensembl IDs from each id refinery file
  unique(index_data$ENSEMBL)
})

#---------------------------Get Illumina Gene lists ---------------------------#
packages <- c("illuminaHumanv1.db", "illuminaHumanv2.db", "illuminaHumanv3.db",
              "illuminaHumanv4.db")

# For each illumina package, get out a gene list:
genes.per.illum <- lapply(packages, function(x){
  if (!(x %in% installed.packages())) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(x, suppressUpdates = TRUE)
  }
  # Load library 
  library(x, character.only = TRUE)
  
  # Get the entrez gene IDs that are mapped to an Ensembl ID
  obj <- eval(parse(text = gsub(".db$", "ENSEMBL", x)))
  mapped.genes <- mappedkeys(obj)
  
  # Convert to a list
  genes <- data.frame(obj[mapped.genes])
  
  # Print out ensembl IDs only
  genes$ensembl_id
})

#---------------------------Save Gene Lists to an RData -----------------------#
save(list = c("genes.per.illum", "genes.per.affy"), file = "genes.per.array.RData")

