FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
ARG ARIA2_VERSION=1.36.0

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        # Build only dependencies
        ## diamond
        g++ \
        gdb \
        libclang-common-6.0-dev \
        cmake\
        g++\
        zlib1g-dev \
        ## minimap2
        build-essential \
        libz-dev \
        # Runtime dependencies
        make \
        curl \
        git \
        wget \
        liblz4-tool \
        pigz \
        python-pandas \
        python3-gdbm \
        python3-pip \
        software-properties-common

# install aria2
WORKDIR /tmp
RUN curl -L https://github.com/aria2/aria2/releases/download/release-${ARIA2_VERSION}/aria2-${ARIA2_VERSION}.tar.gz -o aria2-${ARIA2_VERSION}.tar.gz && \
    tar xzvf aria2-${ARIA2_VERSION}.tar.gz
RUN cd /tmp/aria2-${ARIA2_VERSION} && ./configure && make && mv ./src/aria2c /usr/local/bin/

# install minimap2 and seqkit
RUN git clone --single-branch --branch distributed-mapping  https://github.com/mlin/minimap2.git \
&& cd minimap2 && make && mv ./minimap2 /usr/local/bin/
RUN curl -L https://github.com/shenwei356/seqkit/releases/download/v2.0.0/seqkit_linux_amd64.tar.gz | tar xz && mv seqkit /usr/local/bin/

# install diamond 
WORKDIR /tmp
RUN git clone https://github.com/morsecodist/diamond
WORKDIR /tmp/diamond
RUN git checkout scatter-gather
WORKDIR /tmp/diamond/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j6
RUN mv diamond /usr/local/bin
WORKDIR /workspace

# install boto 3 and awscli
RUN python3 -m pip install boto3==1.12.27 botocore==1.15.39
RUN python3 -m pip install awscli==1.18.39 --upgrade
ENV PATH="${PATH}:/usr/local/bin"

COPY *.py /usr/local/bin/

COPY ncbi_download /usr/local/bin/

WORKDIR /workspace
