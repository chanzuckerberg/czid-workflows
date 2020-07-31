FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="IDseq Team idseq-tech@chanzuckerberg.com"
LABEL description = "Image for consensus genome by metagenomic sequencing with spiked primer enrichment or amplicon sequencing"


RUN apt-get -qq update && apt-get -qq -y install curl locales zip \
  && locale-gen en_US.UTF-8

# Install miniconda and configure it
ENV MINICONDA_SHA256 bb2e3cedd2e78a8bb6872ab3ab5b1266a90f8c7004a22d8dc2ea5effeb6a439a

RUN curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh -o /tmp/miniconda.sh \
  && echo "$MINICONDA_SHA256 /tmp/miniconda.sh" | sha256sum -c - \
  && bash /tmp/miniconda.sh -bfp /usr/local \
  && rm -rf /tmp/miniconda.sh \
  && conda update conda \
  && conda install -y python=3.7.4 \
  && conda clean -a

# Setup conda environment
RUN conda config --add channels defaults \
  && conda config --add channels conda-forge \
  && conda config --add channels bioconda \
  && conda clean -a

RUN conda install -y bcftools=1.9 biopython=1.76 freebayes=1.3.2 ivar=1.2 kraken2=2.0.8_beta minimap2=2.17 \
  multiqc=1.8 pysam=0.15.4 quast=5.0.2 samtools=1.9 seaborn=0.10.0 seqtk=1.3 trim-galore=0.6.5 \
  && conda clean -a
