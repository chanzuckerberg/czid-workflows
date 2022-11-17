# syntax=docker/dockerfile:1.4
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="IDseq Team idseq-tech@chanzuckerberg.com"
LABEL description = "Image for consensus genome by metagenomic sequencing with spiked primer enrichment or amplicon sequencing"

RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list; \
    echo 'APT::Install-Recommends "false";' > /etc/apt/apt.conf.d/98czid; \
    echo 'APT::Install-Suggests "false";' > /etc/apt/apt.conf.d/99czid

RUN apt-get -qq update && apt-get -qq -y install \
    apt-utils \
    build-essential \
    cmake \
    curl \
    g++ \
    gdb \
    git \
    httpie \
    jq \
    libfreetype6-dev \
    libhts-dev \
    libz-dev \
    locales \
    make \
    minimap2 \
    moreutils \
    pkg-config \
    python3-dev \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    samtools \
    seqkit \
    seqtk \
    unzip \
    wget \
    zip \
    zlib1g-dev \
    && locale-gen en_US.UTF-8

RUN pip3 install boto3 marisa-trie==0.7.7 requests biopython pandas pysam

ADD https://raw.githubusercontent.com/chanzuckerberg/miniwdl/v1.1.5/examples/clean_download_cache.sh /usr/local/bin
RUN chmod +x /usr/local/bin/clean_download_cache.sh

RUN curl -Ls https://github.com/chanzuckerberg/s3parcp/releases/download/v0.2.0-alpha/s3parcp_0.2.0-alpha_Linux_x86_64.tar.gz | tar -C /usr/bin -xz s3parcp

RUN wget -O fastp http://opengene.org/fastp/fastp.0.23.1 && \
    chmod a+x ./fastp && \
    cp fastp /bin/fastp

RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.3.0/seqkit_linux_386.tar.gz && \
    tar -xvf seqkit_linux_386.tar.gz && \
    rm seqkit_linux_386.tar.gz && \
    cp seqkit /usr/local/bin
ENV PATH=$PATH:$(pwd)

ARG FLYE_VER="2.9"
RUN wget https://github.com/fenderglass/Flye/archive/${FLYE_VER}.tar.gz && \
    tar -xvf ${FLYE_VER}.tar.gz && \
    rm ${FLYE_VER}.tar.gz && \
    cd Flye-${FLYE_VER} && \
    python3 setup.py build && \
    python3 setup.py install

WORKDIR /tmp
RUN git clone --recursive https://github.com/mlin/minimap2-scatter.git
WORKDIR /tmp/minimap2-scatter
RUN make minimap2
RUN mv /tmp/minimap2-scatter/minimap2/minimap2 /usr/local/bin/minimap2-scatter

WORKDIR /tmp
RUN git clone https://github.com/morsecodist/diamond
WORKDIR /tmp/diamond
RUN git checkout minimal-mods
WORKDIR /tmp/diamond/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j6
RUN mv diamond /usr/local/bin

COPY --from=lib idseq-dag /tmp/idseq-dag
RUN pip3 install /tmp/idseq-dag && rm -rf /tmp/idseq-dag

COPY --from=lib idseq_utils /tmp/idseq_utils
RUN pip3 install /tmp/idseq_utils && rm -rf /tmp/idseq_utils

COPY --from=lib bin /usr/local/bin
