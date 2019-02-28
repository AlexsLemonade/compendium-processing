#!/bin/bash

# use ccdl/dr_salmon:v1.4.9
# like so: docker run -it --volume $(pwd):/media ccdl/dr_salmon:v1.4.9 bash

# Generate Danio_rerio.GRCz10.91_tx2gene.tsv
Rscript scripts/get_tx2gene.R \
	--gtf_file transcriptome_index/Danio_rerio.GRCz10.91.gtf \
	--output_file transcriptome_index/Danio_rerio.GRCz10.91_tx2gene.tsv

# get gene-level estimates with tximport
mkdir tximport_output

# first experiment
Rscript scripts/tximport.R \
	--file_list quants/DRP003275/ \
	--gene2txmap transcriptome_index/Danio_rerio.GRCz10.91_tx2gene.tsv \
	--rds_file tximport_output/DRP003275_tximport.RDS \
	--tpm_file tximport_output/DRP003275_lengthScaledTPM.tsv

# second experiment
Rscript scripts/tximport.R \
	--file_list quants/SRP128941/ \
	--gene2txmap transcriptome_index/Danio_rerio.GRCz10.91_tx2gene.tsv \
	--rds_file tximport_output/SRP128941_tximport.RDS \
	--tpm_file tximport_output/SRP128941_lengthScaledTPM.tsv
