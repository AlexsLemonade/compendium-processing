# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Identify which supported microarray platforms are Homo sapiens.
#  
#------------------- Get the supported platforms file from GitHub--------------#
# magrittr pipe
`%>%` <- dplyr::`%>%`

# Read in the file:
platforms <- read.csv("data/supported_microarray_platforms.csv")

#------------------------ Get GEO list of human arrays-------------------------# 
# Install and attach this package:
if (!("GEOmetadb" %in% installed.packages())) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("GEOmetadb")
}
library(GEOmetadb)

# Check if the database is downloaded
if(!file.exists('data/GEOmetadb.sqlite')){
  getSQLiteFile(destdir = "data")
}

# If it has been, connect to it
if(file.exists('data/GEOmetadb.sqlite')){
  con <- dbConnect(SQLite(),'data/GEOmetadb.sqlite')
}
# Get homo sapiens gpls
hs.gpl <- dbGetQuery(con, paste0('SELECT gpl FROM gpl WHERE gpl.organism LIKE "%homo%"'))

# Get info about these experiments
hs.gpl <- geoConvert(hs.gpl$gpl, out_type = c("gse", "gpl", "gsm"),
                     sqlite_db_name = "data/GEOmetadb.sqlite")

# Save this info for future use
save(hs.gpl, file = "GEO_exp_info.RData")

#------------------Get ArrayExpress list of human arrays-----------------------# 
# Read in the array express human experiments meta info
hs.array <- read.table("data/array_express_human.txt", sep = "\t", skip = 1, stringsAsFactors = FALSE)

#---------------Only keep supported platforms that are human-------------------#
#
# List of both GEO and array express accessions
hs.acc <- c(unique(hs.gpl$gse$from_acc), hs.array$V1)

# Make dataframe with only human acc
platforms <- platforms %>% dplyr::filter(external_accession %in% hs.acc)

# Write info to a csv file
write.csv(platforms, file = "data/exp_acc_human_only.csv", quote = FALSE,
          row.names = FALSE)
