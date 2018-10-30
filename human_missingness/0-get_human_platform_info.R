# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# Purpose: Identify which supported microarray platforms are Homo sapiens.
#
# This script is intended to be ran from using the bash script:
# `human_missingness/run_scripts.sh`
# 
#------------------- Get the supported platforms file from GitHub--------------#
# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Directory set up
results.dir <- "results"
plots.dir <- "results/plots"

# Read in the file:
platforms <- read.csv(file.path("data", "supported_microarray_platforms.csv"))

#------------------------ Get GEO list of human arrays-------------------------# 
# Install and attach this package:
if (!("GEOmetadb" %in% installed.packages())) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("GEOmetadb")
}
library(GEOmetadb)

geo.db <- file.path("data", "GEOmetadb.sqlite")
if (!file.exists(geo.db)) {
# Check if the database is downloaded
  getSQLiteFile(destdir = "data")
}

# If it has been, connect to it
if(file.exists(geo.db)){
  con <- dbConnect(SQLite(),geo.db)
}
# Get homo sapiens gpls
hs.gpl <- dbGetQuery(con, paste0('SELECT gpl FROM gpl WHERE gpl.organism LIKE "%homo%"'))

# Get info about these experiments
hs.gpl <- geoConvert(hs.gpl$gpl, out_type = c("gse", "gpl", "gsm"),
                     sqlite_db_name = geo.db)

# Save this info for future use
save(hs.gpl, file = "results/GEO_exp_info.RData")

#------------------Get ArrayExpress list of human arrays-----------------------# 
# Read in the array express human experiments meta info
hs.array <- read.table(file.path("data", "array_express_human.txt"), sep = "\t",
                       skip = 1, stringsAsFactors = FALSE)

#---------------Only keep supported platforms that are human-------------------#
#
# Get a list of erxperiment accessions by combining the GEO and array acc lists
# obtained from above. 
hs.acc <- c(unique(hs.gpl$gse$from_acc), hs.array$V1)

# Make dataframe with only human acc
platforms <- platforms %>% dplyr::filter(external_accession %in% hs.acc)

# Write info to a csv file
write.csv(platforms, file = file.path("data", "exp_acc_human_only.csv"),
          quote = FALSE, row.names = FALSE)
