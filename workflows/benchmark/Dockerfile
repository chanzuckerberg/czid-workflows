# syntax=docker/dockerfile:1.4
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="CZID Team idseq-tech@chanzuckerberg.com"
LABEL description = "Image for benchmarking CZ ID workflows"

RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list; \
    echo 'APT::Install-Recommends "false";' > /etc/apt/apt.conf.d/98czid; \
    echo 'APT::Install-Suggests "false";' > /etc/apt/apt.conf.d/99czid

RUN apt-get -qq update && apt-get -qq -y install \
    python3-pip

COPY requirements-dev.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY benchmark_helpers benchmark_helpers
RUN pip3 install benchmark_helpers/