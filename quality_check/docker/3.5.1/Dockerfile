FROM rocker/tidyverse:3.5.1

RUN R -e "devtools::install_version('XML', version = '3.98-1.16')"
RUN R -e "devtools::install_version('ggfortify', version = '0.4.5')"
RUN R -e "devtools::install_version('cowplot', version = '0.9.3')"
RUN R -e "BiocManager::install(c('ensembldb', 'rWikiPathways', 'qvalue', 'GSA', 'org.Dr.eg.db'), update = FALSE)"
RUN R -e "devtools::install_github('wgmao/PLIER', ref = 'a2d4a2aa343f9ed4b9b945c04326bebd31533d4d', dependencies = TRUE)"