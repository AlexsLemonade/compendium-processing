#!/bin/bash

# use ccdl/dr_salmon:v1.4.9
# like so: docker run -it --volume $(pwd):/media ccdl/dr_salmon:v1.4.9 bash

# we can check the read lengths for each experiment like so:
# zcat data/fastq/DRP003275/DRR067132/DRR067132_1.fastq.gz \
#   | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) \
#     {print l, lengths[l]}}'
#
# zcat data/fastq/SRP128941/SRR6466455/SRR6466455_1.fastq.gz \
#    | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) \
#      {print l, lengths[l]}}'

# we can use the long txome index for both of these

for drp_iter in {132..143}
do
  salmon --no-version-check --threads=8 quant \
  	-l A -i transcriptome_index/Danio_rerio_long_index \
  	-1 data/fastq/DRP003275/DRR067${drp_iter}/DRR067${drp_iter}_1.fastq.gz \
  	-2 data/fastq/DRP003275/DRR067${drp_iter}/DRR067${drp_iter}_2.fastq.gz \
  	-o quants/DRP003275/DRR067${drp_iter} \
  	--biasSpeedSamp 5 --seqBias --gcBias --dumpEq
done

for srp_iter in {455..458}
do
  salmon --no-version-check --threads=8 quant \
  	-l A -i transcriptome_index/Danio_rerio_long_index \
  	-1 data/fastq/SRP128941/SRR6466${srp_iter}/SRR6466${srp_iter}_1.fastq.gz \
  	-2 data/fastq/SRP128941/SRR6466${srp_iter}/SRR6466${srp_iter}_2.fastq.gz \
  	-o quants/SRP128941/SRR6466${srp_iter} \
  	--biasSpeedSamp 5 --seqBias --gcBias --dumpEq
done
