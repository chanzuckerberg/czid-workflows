# syntax=docker/dockerfile:1.4
FROM ubuntu:20.04

LABEL maintainer="CZ ID Team <idseq-tech@chanzuckerberg.com>"

RUN apt-get update
RUN apt-get -y install \
    curl \
    git \
    python3 \
    python3-pip \
    gcc \
    make \
    libz-dev \
    libncurses-dev \
    libbz2-dev \
    liblzma-dev \
    g++ \
    zip \
    bedtools \
    dh-autoreconf \
    nasm

RUN ln -s /usr/bin/python3 /usr/bin/python

# Host filtering (2022 version) dependencies
# fastp (libdeflate libisal (dh-autoreconf nasm))
# hisat2
# bowtie2 [already installed]
WORKDIR /tmp
RUN wget -nv -O - https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz | tar zx
RUN cd isa-l-* && ./autogen.sh && ./configure && make -j8 && make install
RUN wget -nv -O - https://github.com/ebiggers/libdeflate/archive/refs/tags/v1.12.tar.gz | tar zx
RUN cd libdeflate-* && make -j8 && make install
RUN ldconfig
RUN git clone https://github.com/mlin/fastp.git && git -C fastp checkout 37edd60
RUN cd fastp && make -j8 && ./fastp test && cp fastp /usr/local/bin
WORKDIR /hisat2
RUN wget -nv -O /tmp/HISAT2.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download \
        && unzip /tmp/HISAT2.zip && mv hisat2-*/* . && rm /tmp/HISAT2.zip
RUN curl -L https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz | tar xz -C /

# Install rgi dependencies
RUN curl -L https://github.com/bbuchfink/diamond/releases/download/v0.8.36/diamond-linux64.tar.gz | tar xz && mv diamond /usr/local/bin
RUN curl -L http://www.unafold.org/download/oligoarrayaux-3.8.tar.gz | tar xz && cd oligoarrayaux-3.8/ && ./configure && make install 
RUN curl -L https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 | tar xj && cd samtools-1.9/ && make && mv samtools /usr/local/bin/
RUN  curl -L https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar xj && cd bwa-0.7.17/ && make &&  mv bwa /usr/local/bin/
RUN curl -L  https://bitbucket.org/genomicepidemiology/kma/get/1.4.3.tar.gz | tar xz && cd genomicepidemiology-kma-47505a7dd4eb/ && make && mv kma kma_index kma_shm kma_update /usr/local/bin/

RUN apt-get update && apt-get -y install ncbi-blast+=2.9.0-2 prodigal=1:2.6.3-4 bamtools=2.5.1+dfsg-5build1 bowtie2=2.3.5.1-6build1 bedtools=2.27.1+dfsg-4ubuntu1


COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt


RUN git clone https://github.com/rzlim08/rgi.git
RUN mkdir -p rgi/inputs/heatmap_inputs
RUN cd rgi/ && git checkout 20b22dab2676e6cdfd16245dd15b9af3d266cb1b && python setup.py build && python setup.py test && python setup.py install

# load databases
RUN curl -LO https://czid-public-references.s3.us-west-2.amazonaws.com/test/AMRv2/card.json && \
    rgi card_annotation -i card.json && \
    rgi load \
        --card_annotation card_database_v3.2.3.fasta \
        -i card.json 

COPY empty-main-header.txt /tmp/empty-main-header.txt 
