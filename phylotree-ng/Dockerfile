FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="IDseq Team idseq-tech@chanzuckerberg.com"
LABEL description = "Image for IDseq phylotree-ng workflow"

RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list; \
    echo 'APT::Install-Recommends "false";' > /etc/apt/apt.conf.d/98idseq; \
    echo 'APT::Install-Suggests "false";' > /etc/apt/apt.conf.d/99idseq

RUN apt-get -qq update && apt-get -qq -y install \
    jq \
    moreutils \
    curl \
    locales \
    zip \
    unzip \
    httpie \
    zlib1g-dev \
    libhts-dev \
    pkg-config \
    apt-utils \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    python3-yaml \
    python3-dateutil \
    python3-numpy \
    python3-scipy \
    python3-matplotlib \
    python3-pandas \
    python3-seaborn \
    iqtree \
    build-essential \
    && locale-gen en_US.UTF-8

RUN pip3 install toytree==2.0.1 toyplot==0.19.0 taxoniq==0.6.0 && \
    pip3 install --upgrade \
    https://github.com/chanzuckerberg/taxoniq/releases/download/v0.6.0/ncbi_genbank_accession_db-2021.4.10-py3-none-any.whl \
    https://github.com/chanzuckerberg/taxoniq/releases/download/v0.6.0/ncbi_genbank_accession_lengths-2021.4.10-py3-none-any.whl \
    https://github.com/chanzuckerberg/taxoniq/releases/download/v0.6.0/ncbi_genbank_accession_offsets-2021.4.10-py3-none-any.whl

RUN curl -L https://github.com/simonrharris/SKA/archive/refs/tags/v1.0.tar.gz | tar -xvz && \
    make -C SKA-1.0 && \
    make -C SKA-1.0 install
