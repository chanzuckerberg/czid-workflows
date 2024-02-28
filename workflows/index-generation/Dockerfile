# syntax=docker/dockerfile:1.4
ARG DEBIAN_FRONTEND=noninteractive

FROM ubuntu:22.04 AS builder
ARG DEBIAN_FRONTEND
# Install  build dependencies
RUN apt-get update && \
    apt-get install -y \
    curl \
    git \
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

# install minimap2
RUN git clone --single-branch --branch distributed-mapping https://github.com/mlin/minimap2.git \
    && cd minimap2 && make

# install seqkit
RUN curl -L https://github.com/shenwei356/seqkit/releases/download/v2.0.0/seqkit_linux_amd64.tar.gz | tar xz

# install diamond 
RUN git clone --single-branch --branch scatter-gather https://github.com/morsecodist/diamond \
    && mkdir diamond/build && cd diamond/build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j6


FROM ubuntu:22.04 AS ncbi_compress_builder

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    ## rocksdb
    libgflags-dev \
    libsnappy-dev \
    zlib1g-dev \
    libbz2-dev \
    liblz4-dev \
    libzstd-dev \
    clang \
    # needed to download ncbi-blast
    alien

# based on https://github.com/rust-lang/docker-rust/blob/dcb74d779e8a74263dc8b91d58d8ce7f3c0c805b/1.70.0/buster/Dockerfile
ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH \
    RUST_VERSION=1.70.0

RUN set -eux; \
    dpkgArch="$(dpkg --print-architecture)"; \
    case "${dpkgArch##*-}" in \
        amd64) rustArch='x86_64-unknown-linux-gnu'; rustupSha256='0b2f6c8f85a3d02fde2efc0ced4657869d73fccfce59defb4e8d29233116e6db' ;; \
        armhf) rustArch='armv7-unknown-linux-gnueabihf'; rustupSha256='f21c44b01678c645d8fbba1e55e4180a01ac5af2d38bcbd14aa665e0d96ed69a' ;; \
        arm64) rustArch='aarch64-unknown-linux-gnu'; rustupSha256='673e336c81c65e6b16dcdede33f4cc9ed0f08bde1dbe7a935f113605292dc800' ;; \
        i386) rustArch='i686-unknown-linux-gnu'; rustupSha256='e7b0f47557c1afcd86939b118cbcf7fb95a5d1d917bdd355157b63ca00fc4333' ;; \
        *) echo >&2 "unsupported architecture: ${dpkgArch}"; exit 1 ;; \
    esac; \
    url="https://static.rust-lang.org/rustup/archive/1.26.0/${rustArch}/rustup-init"; \
    wget "$url"; \
    echo "${rustupSha256} *rustup-init" | sha256sum -c -; \
    chmod +x rustup-init; \
    ./rustup-init -y --no-modify-path --profile minimal --default-toolchain $RUST_VERSION --default-host ${rustArch}; \
    rm rustup-init; \
    chmod -R a+w $RUSTUP_HOME $CARGO_HOME; \
    rustup --version; \
    cargo --version; \
    rustc --version;

COPY ncbi-compress /ncbi-compress
RUN cd /ncbi-compress && cargo build --release

FROM ubuntu:22.04
ARG DEBIAN_FRONTEND
# Install dependencies

RUN apt-get update && \
    apt-get install -y \
        wget \
        pigz \
        python3-pip \
        mysql-client \
        kraken2 \
        curl \
        jq

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast-2.15.0+-x64-linux.tar.gz && \
    cp ncbi-blast-2.15.0+/bin/* /usr/local/bin/

COPY --from=builder /diamond/build/diamond /usr/local/bin
COPY --from=builder /minimap2/minimap2 /usr/local/bin/
COPY --from=builder /seqkit /usr/local/bin/
COPY --from=ncbi_compress_builder /ncbi-compress/target/release/ncbi-compress /usr/local/bin

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY *.py /usr/local/bin/

WORKDIR /workspace
