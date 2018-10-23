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
geo_tables <- dbListTables(con)
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

# Write info to a csv file
write.csv(platforms %>% filter(external_accession %in% hs.acc),
          file = "exp_acc_human_only.csv", quote = FALSE, row.names = FALSE)
