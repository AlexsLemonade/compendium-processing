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
# Get a list of the GPLs from the platforms list
illum.gpls <- platforms$external_accession[grep("Illumina", platforms$internal_accession)]

# For each GPL, obtain gene accession ids
genes.per.illum <- lapply(illum.gpls, function(x) {
  
  # Download the GPL file from GEO
  gpl.df <- GEOquery::getGEO(x, destdir = "data")
  
  # Take a look at the GPL table
  gpl.df <- gpl.df@dataTable@table
  
  if (!is.null(gpl.df$GB_ACC)) {
    # If this vector exists, return it 
    return(gpl.df$GB_ACC)
  } else { # Otherwise, search for a column with "NM_" RefSeq ids
    genes <- apply(gpl.df, 2, function(x) any(grepl("^NM_[0-9]*", x)))
  
    # Return only columnns with RefSeq "NM_" accession ids
    return(gpl.df[, genes])
  }
})
# Keep the platform names
names(genes.per.illum) <- grep("Illumina", platforms$internal_accession, value = TRUE)

# For each platform id, get a list of the unique gene identifiers 
# and convert them to ENSEMBL Ids
genes.per.illum <- as.list(tapply(genes.per.illum, names(genes.per.illum), function(x) {
  
    # Need to get rid of the version decimal point
    nm_ids <- gsub("\\.[0-9]*", "", unique(unlist(x)))
    
    # Convert the NM accessions to ensembl ids
    ensembl <- mapIds(org.Hs.eg.db, keys = nm_ids, column = "ENSEMBL",
           keytype = "REFSEQ")
    ensembl <- unlist(ensembl)
    
    # Return non NA ensembl genes
    return(ensembl[!is.na(ensembl)])
  }))

#---------------------------Save Gene Lists to an RData -----------------------#
save(list = c("genes.per.illum", "genes.per.affy", "id.files"), file = "genes.per.array.RData")
