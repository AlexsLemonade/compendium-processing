FROM rocker/tidyverse:3.4.4

RUN R -e "BiocInstaller::biocLite(c('limma', 'rtracklayer', 'GenomicRanges', 'Rsamtools'), suppressUpdates = TRUE)"
