# syntax=docker/dockerfile:1.4
ARG ARIA2_VERSION=1.36.0
ARG DEBIAN_FRONTEND=noninteractive

FROM ubuntu:22.04 AS builder
ARG ARIA2_VERSION
ARG DEBIAN_FRONTEND
# Install  build dependencies
RUN apt-get update && \
    apt-get install -y \
    curl \
    git \
    ## aria2
    libssl-dev \
    liblz4-tool \
    software-properties-common \
    pkg-config \
    ## diamond
    g++ \
    gdb \
    libclang-common-12-dev \
    cmake\
    g++\
    zlib1g-dev \
    ## minimap2
    build-essential \
    libz-dev \
    zip

# install aria2
RUN curl -L https://github.com/aria2/aria2/releases/download/release-${ARIA2_VERSION}/aria2-${ARIA2_VERSION}.tar.gz -o aria2-${ARIA2_VERSION}.tar.gz && \
    tar xzvf aria2-${ARIA2_VERSION}.tar.gz
RUN cd /aria2-${ARIA2_VERSION} && ./configure --without-gnutls --with-openssl && make

# install minimap2
RUN git clone --single-branch --branch distributed-mapping https://github.com/mlin/minimap2.git \
    && cd minimap2 && make

# install seqkit
RUN curl -L https://github.com/shenwei356/seqkit/releases/download/v2.0.0/seqkit_linux_amd64.tar.gz | tar xz

# install diamond 
RUN git clone --single-branch --branch scatter-gather https://github.com/morsecodist/diamond \
    && mkdir diamond/build && cd diamond/build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j6

# install aws cli
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip && ./aws/install -i /usr/local/aws-cli -b /usr/local/bin

FROM ubuntu:22.04
ARG ARIA2_VERSION
ARG DEBIAN_FRONTEND
# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        wget \
        pigz \
        python3-pip \
        mysql-client \
        kraken2 \
        jq

COPY --from=builder /diamond/build/diamond /usr/local/bin
COPY --from=builder /minimap2/minimap2 /usr/local/bin/
COPY --from=builder /aria2-${ARIA2_VERSION}/src/aria2c /usr/local/bin/
COPY --from=builder /seqkit /usr/local/bin/
COPY --from=builder /usr/local/aws-cli /usr/local/aws-cli
COPY --from=builder /usr/local/bin/aws /usr/local/bin

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY *.py /usr/local/bin/
COPY ncbi_download /usr/local/bin/

WORKDIR /workspace
