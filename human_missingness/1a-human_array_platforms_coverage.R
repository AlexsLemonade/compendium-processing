# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Purpose: Identify genes that are measured per array platform 
# 
# This script is intended to be ran from using the bash script:
# `human_missingness/run_scripts.sh`
# 
#-----------------------------------Set Up-------------------------------------#
# Magrittr pipe
`%>%` <- dplyr::`%>%`

library(org.Hs.eg.db)

#------------------Set Up Functions------------------#
getEnsgFromGz <- function (id.file) {
# Purpose: Read the gzipped identifier refinery file and return the ensembl gene ids
#          within the file 

  # Args: 'id.file' : Path to id_refinery zip file
  
  # Returns: List of unique ensembl gene IDs
  
  # Formulate the command to read the input file. 
  file.cmd <- paste0('gzcat ', id.file)
  
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

getGPLsRefSeqIds <- function(gpl = gpl) {
  # Purpose: Given a GPL ID, Retrieve RefSeq accession ID that corresponds to it.
   
  # Args: 'gpl' : A experiment platform id from gene expression omnibus
  
  # Returns: 1) A vector of RefSeq gene ids is returned (if it exists in the GPL
  #            file)
  #          2) GPL file is downloaded to "data" directory 
  
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

convertRefSeqtoENSG <-  function(refseq.ids) {
  # Purpose: Convert RefSeq accession ids to ensembl gene ids

  # Args: 'refseq.ids' : a vector of refseq ids (typically look like: NM_12345)

  # Returns: A vector of corresponding ensembl gene ids
  
  # Keep unique NM_ids, make them a vector
  refseq.ids <- unique(unlist(refseq.ids))
  
  # Get rid of the version decimal point
  refseq.ids <- gsub("\\.[0-9]*", "", refseq.ids)
  
  # Convert the NM accessions to ensembl ids
  ensembl <- mapIds(org.Hs.eg.db, keys = refseq.ids, column = "ENSEMBL",
                    keytype = "REFSEQ")
  # We want a vector
  ensembl <- unlist(ensembl)
  
  # Return non NA ensembl genes
  return(ensembl[!is.na(ensembl)])
}

#------------------ Directory set up------------------#
results.dir <- "results"
plots.dir <- file.path("results","plots")

#----------------------------Get Affy Gene Lists-------------------------------#
# Read in human platform list
platforms <- read.csv(file.path("data", "exp_acc_human_only.csv"), 
                      stringsAsFactors = FALSE)

# Get list of id.refinery files
id.refinery <- dir(file.path("data", "id.refinery"))

# Only keep those for the platforms we care about (human in this case)
tmp <- match(unique(platforms$internal_accession), gsub(".tsv.gz", "", id.refinery))

# Make vector of file paths to human id_refinery files
id.files <- file.path("data", "id.refinery", id.refinery[tmp[!is.na(tmp)]])

# For each id refinery file that matches our list, read the file and return the 
# ensembl IDs
genes.per.affy <- lapply(id.files, getEnsgFromGz)

# Bring the names of the platforms
names(genes.per.affy) <- gsub(".tsv.gz", "", id.refinery[tmp[!is.na(tmp)]])

#-------------------Get Illumina Gene lists from GPL files --------------------#
# Get a list of the GPLs from the platforms list
illum.gpls <- platforms$external_accession[grep("Illumina", 
                                                platforms$internal_accession)]

# For each GPL, obtain gene accession ids
genes.per.illum <- lapply(illum.gpls, getGPLsRefSeqIds)

# Keep the platform names
names(genes.per.illum) <- grep("Illumina", platforms$internal_accession,
                               value = TRUE)

#------------Convert Illumina Gene NM_ ids to ENSG's per platform--------------#

# For each platform id, get a list of the unique gene identifiers 
# and convert them to ENSEMBL Ids

# Make one gene list per platform that is ENSG's.
# This is collapsing multiple gene lists into a gene unique list of genes all 
# the gpls in the set 
genes.per.illum <- as.list(tapply(genes.per.illum, names(genes.per.illum),
                                  convertRefSeqtoENSG))

#---------------------------Save Gene Lists to an RData -----------------------#
saveRDS(genes.per.illum, file = file.path(results.dir,
                                          "genes_per_illumina_array.RDS"))
saveRDS(genes.per.affy, file = file.path(results.dir,
                                         "genes_per_affy_array.RDS"))
