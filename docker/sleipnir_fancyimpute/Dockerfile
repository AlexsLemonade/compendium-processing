FROM jupyter/tensorflow-notebook:6c85e4b43a26

RUN pip install --ignore-installed 'numpy==1.15.2' 'cvxpy==1.0.10' 'scikit-learn==0.20.0'
RUN pip install 'fancyimpute==0.4.0'

# Adapted from Greg Way
# https://github.com/greenelab/hgsc_subtypes/blob/9986ddd51c557bf4e1222998d7871ff635475a47/docker/Dockerfile

USER root

# Install Ubuntu dependencies
RUN sudo apt-get update && apt-get install -y --no-install-recommends \
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
RUN sudo apt-get install -y --no-install-recommends \
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

USER $NB_UID
