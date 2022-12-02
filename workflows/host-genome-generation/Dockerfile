# syntax=docker/dockerfile:1.4
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

# this brings in python2.7
RUN apt-get update && apt-get install -y bowtie2 curl minimap2

# Install STAR, the package rna-star does not include STARlong
RUN curl -L https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz | tar xz
RUN mv STAR-2.5.3a/bin/Linux_x86_64_static/* /usr/local/bin
RUN rm -rf STAR-2.5.3a
