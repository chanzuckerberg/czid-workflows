FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="IDseq Team idseq-tech@chanzuckerberg.com"
LABEL description = "Image for consensus genome by metagenomic sequencing with spiked primer enrichment or amplicon sequencing"

RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list; \
    echo 'APT::Install-Recommends "false";' > /etc/apt/apt.conf.d/98czid; \
    echo 'APT::Install-Suggests "false";' > /etc/apt/apt.conf.d/99czid

RUN apt-get -qq update && apt-get -qq -y install \
    wget \
    libz-dev \
    build-essential \
    make \
    jq \
    moreutils \
    curl \
    locales \
    zip \
    unzip \
    httpie \
    zlib1g-dev \
    libfreetype6-dev \
    libhts-dev \
    pkg-config \
    apt-utils \
    minimap2 \
    samtools \
    seqtk \
    seqkit \
    python3-dev \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    && locale-gen en_US.UTF-8

# Match versions from Ubuntu 20.04 LTS (please do not upgrade except to fix bugs)
RUN pip3 install boto3==1.9.253 marisa-trie==0.7.7 requests biopython pandas

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

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.3-beta-Linux_x86_64.zip && \
    unzip centrifuge-1.0.3-beta-Linux_x86_64.zip
ENV PATH="/centrifuge-1.0.3-beta:$PATH"

COPY idseq-dag /tmp/idseq-dag
RUN pip3 install /tmp/idseq-dag && rm -rf /tmp/idseq-dag

COPY idseq_utils /tmp/idseq_utils
RUN pip3 install /tmp/idseq_utils && rm -rf /tmp/idseq_utils
