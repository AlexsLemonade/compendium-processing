#!/bin/bash

# use ccdl/dr_salmon:v1.4.9
# like so: docker run -it --volume $(pwd):/media ccdl/dr_salmon:v1.4.9 bash

# here's how Danio_rerio.GRCz10.91_tx2gene.tsv was generated in R
#
# gtf_file <- "transcriptome_index/Danio_rerio.GRCz10.91.gtf"
# output_file <- "transcriptome_index/Danio_rerio.GRCz10.91_tx2gene.tsv"
# db_file <- sub(".gtf", ".sqlite", gtf_file)
# if (!file.exists(db_file)) {
#   ensembldb::ensDbFromGtf(gtf_file, outfile = db_file)
# }
# edb <- ensembldb::EnsDb(db_file)
# tx <- ensembldb::transcriptsBy(edb, by = "gene")
# tx.df <- as.data.frame(tx@unlistData)
# tx2gene <- tx.df[, c("tx_name", "gene_id")]
# readr::write_tsv(tx2gene, path = output_file)

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
