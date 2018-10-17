FROM ubuntu:16.04

# Adapted from Greg Way
# https://github.com/greenelab/hgsc_subtypes/blob/9986ddd51c557bf4e1222998d7871ff635475a47/docker/Dockerfile

# Install Ubuntu dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
		git \
		ed \		
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		bison \
		flex \
		apt-utils \
		gawk \
		libcurl4-openssl-dev \
		libxml2-dev \
		libssl-dev

# Install Sleipnir
RUN apt-get install -y --no-install-recommends \
		mercurial \
		gengetopt \
		libboost-regex-dev \
		libboost-graph-dev \
		liblog4cpp5-dev \
		build-essential \
		libgsl0-dev \
		automake

RUN hg clone https://bitbucket.org/libsleipnir/sleipnir
RUN cd sleipnir  && ./gen_auto && ./gen_tools_am
RUN cd sleipnir && ./configure && make && make install
