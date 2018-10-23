# C. Savonen 
# CCDL for ALSF 
# 10/22/18
# 
# What genes are measured using Ensembl gene ids from identifier refinery?
#  
# Find out for Affymetrix platforms and the illumanHumanvx.db for Illumina BeadArrays 
# Note that BeadArrays can be tricky because there's a "whole genome" and a 
# "ref seq" for (almost) all version numbers.
#     
# Read in the human accession platform IDs that are supported by refine.bio
platforms <- read.csv("exp_acc_human_only.csv")

# Obtain a list of the unique platform ids
platforms <- unique(platforms$internal_accession)

# this is ~1.4GB so it may take a bit
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
tmp <- match(platforms, gsub(".tsv.gz", "", id.refinery))
id.files <- id.refinery[tmp[!is.na(tmp)]]
id.files <- paste0("id_refinery/", id.files)

genes.per.platform <- lapply( id.files, function (x) {
                      index_path <- paste('zcat', x)
     suppressWarnings(index_data <- data.table::fread(
                                     cmd = index_path,
                                     stringsAsFactors=FALSE,
                                     sep="\t",
                                     header=TRUE,
                                     autostart=10,
                                     data.table=FALSE,
                                     check.names=FALSE,
                                     fill=TRUE,
                                     na.strings="",
                                     showProgress=FALSE)
)})

# Because this is a huge object that takes a long time to make, let's save it 
# to an RData file so we don't have to make it again hopefully. 
save(genes.per.platform, file = "genes.per.platform.RData")
