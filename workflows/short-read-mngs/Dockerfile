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

RUN go get golang.org/x/tools/cmd/goimports@v0.2.0
RUN go install golang.org/x/tools/cmd/goimports
RUN go get github.com/go-python/gopy@v0.4.4
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
        docker.io \
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

RUN pip3 install boto3==1.23.10 marisa-trie==0.7.7 pytest requests==2.27.1
RUN pip3 install miniwdl==${MINIWDL_VERSION} miniwdl-s3parcp==0.0.5 miniwdl-s3upload==0.0.4
RUN pip3 install https://github.com/chanzuckerberg/miniwdl-plugins/archive/f0465b0.zip#subdirectory=sfn-wdl
RUN pip3 install https://github.com/chanzuckerberg/s3mi/archive/v0.8.0.tar.gz

ADD https://raw.githubusercontent.com/chanzuckerberg/miniwdl/v${MINIWDL_VERSION}/examples/clean_download_cache.sh /usr/local/bin
RUN chmod +x /usr/local/bin/clean_download_cache.sh

# docker.io is the largest package at 250MB+ / half of all package disk space usage.
# The docker daemons never run inside the container - removing them saves 150MB+
RUN rm -f /usr/bin/dockerd /usr/bin/containerd*

RUN cd /usr/bin; curl -O https://amazon-ecr-credential-helper-releases.s3.amazonaws.com/0.4.0/linux-amd64/docker-credential-ecr-login
RUN chmod +x /usr/bin/docker-credential-ecr-login
RUN mkdir -p /root/.docker
RUN jq -n '.credsStore="ecr-login"' > /root/.docker/config.json

RUN curl -L -o /usr/bin/czid-dedup https://github.com/chanzuckerberg/czid-dedup/releases/download/v0.1.2/czid-dedup-Linux; chmod +x /usr/bin/czid-dedup

# Note: bsdtar is available in libarchive-tools
# Note: python3-scipy pulls in gcc (fixed in Ubuntu 19.10)
# TODO: kSNP3 (separate phylotree image?)

# Note: the NonHostAlignment stage uses a different version of gmap custom to CZ ID, installed here:
# https://github.com/chanzuckerberg/czid/blob/master/workflows/docker/gsnap/Dockerfile#L16-L20
# TODO: migrate both to https://packages.ubuntu.com/focal/gmap (updates to gmap require revalidation)
RUN apt-get -q install -y gmap

# FIXME: replace trimmomatic with cutadapt (trimmomatic pulls in too many deps)
RUN apt-get -q install -y trimmomatic
RUN ln -sf /usr/share/java/trimmomatic-0.36.jar /usr/local/bin/trimmomatic-0.38.jar

# FIXME: replace PriceSeqFilter with cutadapt quality/N-fraction cutoff
RUN curl -s https://idseq-prod-pipeline-public-assets-us-west-2.s3-us-west-2.amazonaws.com/PriceSource140408/PriceSeqFilter > /usr/bin/PriceSeqFilter
RUN chmod +x /usr/bin/PriceSeqFilter

RUN curl -Ls https://github.com/chanzuckerberg/s3parcp/releases/download/v0.2.0-alpha/s3parcp_0.2.0-alpha_Linux_x86_64.tar.gz | tar -C /usr/bin -xz s3parcp

# FIXME: check if use of pandas, pysam is necessary
RUN pip3 install pysam==0.14.1 pandas==1.1.5

# Picard for average fragment size https://github.com/broadinstitute/picard
# r-base is a dependency of collecting input size metrics https://github.com/bioconda/bioconda-recipes/pull/16398
RUN apt-get install -y r-base
RUN curl -L -o /usr/local/bin/picard.jar https://github.com/broadinstitute/picard/releases/download/2.21.2/picard.jar
# Create a single executable so we can use SingleCommand
RUN printf '#!/bin/bash\njava -jar /usr/local/bin/picard.jar "$@"\n' > /usr/local/bin/picard
RUN chmod +x /usr/local/bin/picard

# install STAR, the package rna-star does not include STARlong
RUN curl -L https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz | tar xz
RUN mv STAR-2.5.3a/bin/Linux_x86_64_static/* /usr/local/bin
RUN rm -rf STAR-2.5.3a

# install gsnap and rapsearch2 for local alignment
RUN apt-get install -y \
        g++ \
        libperl4-corelibs-perl \
        make \
        # Runtime dependencies
        # On versions higher than 1.62 compilation errors because of an overloading ambiguity with the `advance` method
        libboost1.62-all-dev


# install gsnap, we are using a pinned version for alignment
#  we are also using a different version for host filtering so we append -208-10-26 to each binary.
WORKDIR /gmap-gsnap-2018-10-26
RUN curl -L  https://idseq-gsnap-prerelease.s3-us-west-2.amazonaws.com/gmap-2018-10-26.tar.gz | tar xz -C /gmap-gsnap-2018-10-26 --strip-components 1
RUN CFLAGS="-O3 -mavx2 -mssse3" ./configure --prefix=/gmap-gsnap-2018-10-26 --with-simd-level=avx2 && make -j 8 && make install
RUN ls bin | xargs -I % mv bin/% bin/%-2018-10-26
ENV PATH="${PATH}:/gmap-gsnap-2018-10-26/bin/"

WORKDIR /rapsearch2/Src
RUN curl -L  https://idseq-rapsearch2.s3-us-west-2.amazonaws.com/RAPSearch2.24_64bits.tar.gz  | tar xz -C /rapsearch2 --strip-components 1
# Use the version of boost we installed
RUN sed -i -e 's|^INC.*|INC := -I /usr/include/boost|' -e 's|^LIB.*|LIB :=|' Makefile
RUN make
ENV PATH="${PATH}:/rapsearch2/Src/"

RUN apt-get -y update && apt-get install -y build-essential libz-dev git python3-pip cmake

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

RUN curl -Ls https://github.com/chanzuckerberg/s3parcp/releases/download/v0.2.0-alpha/s3parcp_0.2.0-alpha_Linux_x86_64.tar.gz | tar -C /usr/bin -xz s3parcp

# Host filtering (2022 version) dependencies
# fastp (libdeflate libisal (dh-autoreconf nasm))
# hisat2
# bowtie2 [already installed]
# kallisto
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

# Uninstall build only dependencies
RUN apt-get purge -y g++ libperl4-corelibs-perl make

COPY --from=lib idseq-dag /tmp/idseq-dag
RUN pip3 install /tmp/idseq-dag && rm -rf /tmp/idseq-dag

COPY --from=lib idseq_utils /tmp/idseq_utils
RUN pip3 install /tmp/idseq_utils && rm -rf /tmp/idseq_utils

COPY --from=lib /bin/raise_error /usr/local/bin/raise_error
COPY --from=lib /bin/log_assembly_fail.py /usr/local/bin/log_assembly_fail.py

COPY --from=s3quilt /usr/lib/python3/dist-packages/s3quilt/ /usr/lib/python3/dist-packages/s3quilt/
