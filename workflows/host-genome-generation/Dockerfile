# syntax=docker/dockerfile:1.4
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y wget bowtie2 curl minimap2 pigz dh-autoreconf nasm make git g++ unzip python3-pip

# nescessary for hisat2
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install STAR, the package rna-star does not include STARlong
RUN curl -L https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz | tar xz
RUN mv STAR-2.5.3a/bin/Linux_x86_64_static/* /usr/local/bin
RUN rm -rf STAR-2.5.3a

# Install fastp (libdeflate libisal (dh-autoreconf nasm))
WORKDIR /tmp
RUN wget -nv -O - https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz | tar zx
RUN cd isa-l-* && ./autogen.sh && ./configure && make -j8 && make install
RUN wget -nv -O - https://github.com/ebiggers/libdeflate/archive/refs/tags/v1.12.tar.gz | tar zx
RUN cd libdeflate-* && make -j8 && make install
RUN ldconfig
RUN git clone https://github.com/mlin/fastp.git && git -C fastp checkout 37edd60
RUN cd fastp && make -j8 && ./fastp test && cp fastp /usr/local/bin

# Install hisat2
WORKDIR /hisat2
RUN wget -nv -O /tmp/HISAT2.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download \
    && unzip /tmp/HISAT2.zip && mv hisat2-*/* . && rm /tmp/HISAT2.zip

# Install kallisto + python gtfparse
RUN curl -L https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz | tar xz -C /
RUN pip3 install gtfparse==1.2.1

WORKDIR /