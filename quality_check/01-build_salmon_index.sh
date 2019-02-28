#!/bin/bash

# use ccdl/dr_salmon:v1.4.9
# like so: docker run -it --volume $(pwd):/media ccdl/dr_salmon:v1.4.9 bash

# build indices, long and short
# long
salmon --no-version-check --threads=8 index \
	-t transcriptome_index/Danio_rerio.GRCz10.cdna.all.fa.gz \
	-i transcriptome_index/Danio_rerio_long_index \
	--type quasi -k 31

# short
salmon --no-version-check --threads=8 index \
	-t transcriptome_index/Danio_rerio.GRCz10.cdna.all.fa.gz \
	-i transcriptome_index/Danio_rerio_short_index \
	--type quasi -k 23
