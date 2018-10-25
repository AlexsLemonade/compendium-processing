FROM rocker/tidyverse:3.5.1

RUN R -e "BiocInstaller::biocLite('preprocessCore', 'rhdf5')"
RUN R -e "devtools::install_version('Metrics', version = '0.1.4')"
