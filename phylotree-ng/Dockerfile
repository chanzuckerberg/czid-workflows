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
    iqtree \
    build-essential \
    fonts-open-sans \
    && locale-gen en_US.UTF-8

RUN curl -L https://github.com/simonrharris/SKA/archive/refs/tags/v1.0.tar.gz | tar -xvz && \
    make -C SKA-1.0 && \
    make -C SKA-1.0 install

COPY requirements.txt requirements.txt

RUN pip3 install -r requirements.txt

COPY python_steps/* /bin/
