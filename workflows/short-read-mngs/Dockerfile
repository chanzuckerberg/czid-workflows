# syntax=docker/dockerfile:1.4
FROM ubuntu:18.04 AS s3quilt

RUN apt-get update && apt-get install -y curl python3-pip

ARG GO_VERSION=1.19.2
RUN curl -Lo go${GO_VERSION}.linux-amd64.tar.gz https://go.dev/dl/go${GO_VERSION}.linux-amd64.tar.gz
RUN tar -C . -xzf go${GO_VERSION}.linux-amd64.tar.gz && rm go${GO_VERSION}.linux-amd64.tar.gz
ENV GOPATH /go
ENV PATH $GOPATH/bin:$PATH
RUN /go/bin/go version
WORKDIR /go/src/s3quilt
COPY --from=lib s3quilt/go.mod s3quilt/go.sum s3quilt/s3quilt.go ./

RUN go get golang.org/x/tools/cmd/goimports
RUN go install golang.org/x/tools/cmd/goimports
RUN go get github.com/go-python/gopy
RUN go install github.com/go-python/gopy

RUN pip3 install pybindgen setuptools wheel

RUN gopy build -vm=python3 --output s3quilt .
COPY --from=lib s3quilt/__init__.py s3quilt/
RUN cp -r s3quilt/ /usr/lib/python3/dist-packages/

FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
ARG MINIWDL_VERSION=1.1.5

LABEL maintainer="CZ ID Team <idseq-tech@chanzuckerberg.com>"

RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list; \
        echo 'APT::Install-Recommends "false";' > /etc/apt/apt.conf.d/98czid; \
        echo 'APT::Install-Suggests "false";' > /etc/apt/apt.conf.d/99czid

RUN apt-get -q update && apt-get -q install -y \
        jq \
        moreutils \
        pigz \
        pixz \
        aria2 \
        httpie \
        curl \
        wget \
        zip \
        unzip \
        zlib1g-dev \
        pkg-config \
        apt-utils \
        libbz2-dev \
        liblzma-dev \
        software-properties-common \
        libarchive-tools \
        liblz4-tool \
        lbzip2 \
        python3-dev \
        python3-pip \
        python3-setuptools \
        python3-wheel \
        python3-requests \
        python3-yaml \
        python3-dateutil \
        python3-psutil \
        python3-cutadapt \
        python3-scipy \
        samtools \
        fastx-toolkit \
        seqtk \
        bedtools \
        dh-autoreconf \
        nasm \
        build-essential 

# The following packages pull in python2.7
RUN apt-get -q install -y \
        bowtie2 \
        spades \
        ncbi-blast+

RUN pip3 install boto3==1.23.10 marisa-trie==0.7.7 pytest
RUN pip3 install miniwdl==${MINIWDL_VERSION} miniwdl-s3parcp==0.0.5 miniwdl-s3upload==0.0.4
RUN pip3 install https://github.com/chanzuckerberg/miniwdl-plugins/archive/f0465b0.zip#subdirectory=sfn-wdl
RUN pip3 install https://github.com/chanzuckerberg/s3mi/archive/v0.8.0.tar.gz

ADD https://raw.githubusercontent.com/chanzuckerberg/miniwdl/v${MINIWDL_VERSION}/examples/clean_download_cache.sh /usr/local/bin
RUN chmod +x /usr/local/bin/clean_download_cache.sh

RUN curl -L -o /usr/bin/czid-dedup https://github.com/chanzuckerberg/czid-dedup/releases/download/v0.1.2/czid-dedup-Linux; chmod +x /usr/bin/czid-dedup

RUN curl -Ls https://github.com/chanzuckerberg/s3parcp/releases/download/v0.2.0-alpha/s3parcp_0.2.0-alpha_Linux_x86_64.tar.gz | tar -C /usr/bin -xz s3parcp

# FIXME: check if use of pandas, pysam is necessary
RUN pip3 install pysam==0.14.1 pandas==1.1.5

# Workaround for srst2 refusing to work with upstream bowtie2 and samtools
# FIXME: replace srst2 with a more appropriate tool
RUN apt-get -q install -y srst2
RUN sed -i '/Incorrect version/ s/\w.*/return "1.7"/' /usr/bin/srst2

# Picard for average fragment size https://github.com/broadinstitute/picard
# r-base is a dependency of collecting input size metrics https://github.com/bioconda/bioconda-recipes/pull/16398
RUN apt-get install -y r-base
RUN curl -L -o /usr/local/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/2.21.2/picard.jar
# Create a single executable so we can use SingleCommand
RUN printf '#!/bin/bash\njava -jar /usr/local/bin/picard.jar "$@"\n' > /usr/local/bin/picard
RUN chmod +x /usr/local/bin/picard

# install gsnap and rapsearch2 for local alignment
RUN apt-get install -y \
        g++ \
        libperl4-corelibs-perl \
        make \
        # Runtime dependencies
        # On versions higher than 1.62 compilation errors because of an overloading ambiguity with the `advance` method
        libboost1.62-all-dev


RUN apt-get -y update && apt-get install -y libz-dev git python3-pip cmake

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

# Host filtering (2022 version) dependencies
# fastp (libdeflate libisal (dh-autoreconf nasm))
# hisat2
# bowtie2 [already installed]
# kallisto + python gtfparse
WORKDIR /tmp
RUN wget -nv -O - https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz | tar zx
RUN cd isa-l-* && ./autogen.sh && ./configure && make -j8 && make install
RUN wget -nv -O - https://github.com/ebiggers/libdeflate/archive/refs/tags/v1.12.tar.gz | tar zx
RUN cd libdeflate-* && make -j8 && make install
RUN ldconfig
RUN git clone https://github.com/mlin/fastp.git && git -C fastp checkout c427e57
RUN cd fastp && make -j8 && ./fastp test && cp fastp /usr/local/bin
WORKDIR /
RUN wget -nv -O /tmp/HISAT2.zip https://czid-public-references.s3.us-west-2.amazonaws.com/test/hisat2/hisat2.zip \
        && unzip /tmp/HISAT2.zip && rm /tmp/HISAT2.zip
RUN curl -L https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz | tar xz -C /
RUN pip3 install gtfparse==1.2.1

COPY --from=lib idseq-dag /tmp/idseq-dag
RUN pip3 install /tmp/idseq-dag && rm -rf /tmp/idseq-dag

COPY --from=lib idseq_utils /tmp/idseq_utils
RUN pip3 install /tmp/idseq_utils && rm -rf /tmp/idseq_utils

COPY --from=lib /bin/raise_error /usr/local/bin/raise_error

COPY --from=s3quilt /usr/lib/python3/dist-packages/s3quilt/ /usr/lib/python3/dist-packages/s3quilt/