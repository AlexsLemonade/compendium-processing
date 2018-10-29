# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Purpose: Identify genes that are measured per array platform 
# 
# This script is intended to be ran from using the bash script:
# `human_missingness/run_scripts.sh`
# 
#----------------------------Get Affy Gene Lists-------------------------------#
# Magrittr pipe
`%>%` <- dplyr::`%>%`

library(org.Hs.eg.db)

# Read in human platform list
platforms <- read.csv(file.path("data", "exp_acc_human_only.csv"), stringsAsFactors = FALSE)

# Get list of id.refinery files
id.refinery <- dir("data/id.refinery")

# Only keep those for the platforms we care about (human in this case)
tmp <- match(unique(platforms$internal_accession), gsub(".tsv.gz", "", id.refinery))
id.files <- id.refinery[tmp[!is.na(tmp)]]

# Make a function to read the gzipped id.file and return the ensembl gene ids 
get.ensg.from.gz <- function (id.file) {
  # Formulate the command to read the input file. 
  file.cmd <- paste0('gzcat data/id.refinery/', id.file)
  # Read the table
  suppressWarnings(file.data <- data.table::fread(
    cmd = file.cmd,
    stringsAsFactors = FALSE,
    sep = "\t",
    header=TRUE,
    autostart = 10,
    data.table = FALSE,
    check.names = FALSE,
    fill = TRUE,
    na.strings = "",
    showProgress = FALSE))
  # Return only the unique ensembl IDs from each id refinery file
  unique(file.data$ENSEMBL)
}

# For each id refinery file that matches our list, read the file and return the ensembl IDs
genes.per.affy <- lapply(id.files, get.ensg.from.gz)

# Bring the names of the platforms
names(genes.per.affy) <- gsub(".tsv.gz", "", id.files)

#-------------------Get Illumina Gene lists from GPL files --------------------#
# Get a list of the GPLs from the platforms list
illum.gpls <- platforms$external_accession[grep("Illumina", platforms$internal_accession)]

# Given a GPL ID, Retrieve NM_accessions that correspond to it. 
get.nm.ids <- function(gpl = gpl) {
  
  # Download the GPL file from GEO
  gpl.df <- GEOquery::getGEO(gpl, destdir = "data")
  
  # Take a look at the GPL table
  gpl.df <- gpl.df@dataTable@table
  
  # Usually NM_ids are stored in a column named "GB_ACC" 
  if (!is.null(gpl.df$GB_ACC)) {
    return(gpl.df$GB_ACC)
  
  # If this column doesn't exist, search for the NM_ containing column
  } else {
    genes <- apply(gpl.df, 2, function(gpl) any(grepl("^NM_[0-9]*", gpl)))
    
    # Return only columnns with RefSeq "NM_" accession ids
    return(gpl.df[, genes])
  }
}


# For each GPL, obtain gene accession ids
genes.per.illum <- lapply(illum.gpls, get.nm.ids)

# Keep the platform names
names(genes.per.illum) <- grep("Illumina", platforms$internal_accession, value = TRUE)

#------------Convert Illumina Gene NM_ ids to ENSG's per platform--------------#

# For each platform id, get a list of the unique gene identifiers 
# and convert them to ENSEMBL Ids

# Conversion function to change NM_ accession to ENSG 
convert.nm.ensg <-  function(nm.ids) {
  
  # Keep unique NM_ids, make them a vector
  nm.ids <- unique(unlist(nm.ids))
  
  # Get rid of the version decimal point
  nm.ids <- gsub("\\.[0-9]*", "", nm.ids)

  # Convert the NM accessions to ensembl ids
  ensembl <- mapIds(org.Hs.eg.db, keys = nm.ids, column = "ENSEMBL",
                  keytype = "REFSEQ")
  # We want a vector
  ensembl <- unlist(ensembl)

  # Return non NA ensembl genes
  return(ensembl[!is.na(ensembl)])
}

# Make one gene list per platform that is ENSG's.
# This is collapsing multiple gene lists into a gene unique list of genes all 
# the gpls in the set 
genes.per.illum <- as.list(tapply(genes.per.illum, names(genes.per.illum),
                                  convert.nm.ensg))

#---------------------------Save Gene Lists to an RData -----------------------#
save("genes.per.illum", file = "genes.per.illumina.array.RData")
save("genes.per.affy", file = "genes.per.affy.array.RData")
