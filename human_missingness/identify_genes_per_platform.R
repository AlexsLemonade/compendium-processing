# C. Savonen 
# CCDL for ALSF 
# 
# What genes are measured using Ensembl gene ids from identifier refinery?
#  
# Find out for Affymetrix platforms and the illumanHumanvx.db for Illumina BeadArrays 
# Note that BeadArrays can be tricky because there's a "whole genome" and a 
# "ref seq" for (almost) all version numbers.
#     
download.file("https://zenodo.org/record/1322711/files/all_1532646318.zip?download=1",
              destfile = "identifier_refinery.zip")
unzip("identifier_refinery.zip")