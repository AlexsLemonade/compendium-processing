# C. Savonen 
# CCDL for ALSF 
# 
# Use ARCHS4 and/or recount2 to get an idea of the coverage of each gene in human
# RNA-seq data. (Preference for ARCHS4 because they use kallisto.)
# Assume that we will set 0 values to NA in RNA-seq data.
# Assume also that any of our supported sequencing instruments can measure 
# any gene in our reference transcriptome (probably not 100% true but close enough) 
# when asking "how many samples are there."