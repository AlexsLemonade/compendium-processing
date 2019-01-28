#!/bin/bash

# get refinebio data
mkdir refinebio-data && cd refinebio-data
wget https://s3.amazonaws.com/data-refinery-test-assets/8d102523-af99-4c45-9693-5e9fcb469198_compendia.zip
unzip 8d102523-af99-4c45-9693-5e9fcb469198_compendia.zip
rm 8d102523-af99-4c45-9693-5e9fcb469198_compendia.zip
cd ..

# get data for building salmon indices -- we'll use the cdna data
# this is a bit different from what we do in refine.bio proper
mkdir transcriptome_index 
wget --directory-prefix=transcriptome_index \
	ftp://ftp.ensembl.org/pub/release-91/fasta/danio_rerio/cdna/Danio_rerio.GRCz10.cdna.all.fa.gz
wget --directory-prefix=transcriptome_index \
	ftp://ftp.ensembl.org/pub/release-91/gtf/danio_rerio/Danio_rerio.GRCz10.91.gtf.gz
gunzip transcriptome_index/Danio_rerio.GRCz10.91.gtf.gz

# get the fastq files we're going to use for the comparison
# this is objectively terrible code!

mkdir -p data/fastq && cd data/fastq

# first experiment
mkdir SRP128941 && cd SRP128941
ascp -QT -l 1000m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR646/005/SRR6466455 .
ascp -QT -l 1000m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR646/006/SRR6466456 .
ascp -QT -l 1000m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR646/007/SRR6466457 .
ascp -QT -l 1000m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR646/008/SRR6466458 .

# second experiment 
cd ..
mkdir DRP003275 && cd DRP003275
for i in {132..143}
do
	ascp -QT -l 1000m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/DRR067/DRR067${i} .
done

# FASTA file for GRCz11
cd ../.. && mkdir GRCz11 && cd GRCz11
wget ftp://ftp.ensembl.org/pub/release-95/gff3/danio_rerio/Danio_rerio.GRCz11.95.gff3.gz
gunzip Danio_rerio.GRCz11.95.gff3.gz
# this should get the transcript lengths
# from: https://github.com/stephenturner/oneliners/tree/3233f3cc6c835cc3f70aaf0462f2d2b31244069d#gff3-annotations
# grep $'\tgene\t' Danio_rerio.GRCz11.95.gff3 | cut -s -f 4,5 | perl -ne '@v = split(/\t/); printf("%d\n", $v[1] - $v[0] + 1)'
