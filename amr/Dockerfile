FROM ubuntu:20.04

LABEL maintainer="CZ ID Team <idseq-tech@chanzuckerberg.com>"

RUN apt-get update && apt-get -y install curl git python3 python3-pip gcc make libz-dev libncurses-dev libbz2-dev liblzma-dev g++ 

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN curl -L https://github.com/bbuchfink/diamond/releases/download/v0.8.36/diamond-linux64.tar.gz | tar xz && mv diamond /usr/local/bin
RUN curl -L http://www.unafold.org/download/oligoarrayaux-3.8.tar.gz | tar xz && cd oligoarrayaux-3.8/ && ./configure && make install 
RUN curl -L https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 | tar xj && cd samtools-1.9/ && make && mv samtools /usr/local/bin/
RUN curl -L https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz | tar xz && cd bedtools2 && make && mv bin/* /usr/local/bin/
RUN  curl -L https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 | tar xj && cd bwa-0.7.17/ && make &&  mv bwa /usr/local/bin/
RUN curl -L  https://bitbucket.org/genomicepidemiology/kma/get/1.4.3.tar.gz | tar xz && cd genomicepidemiology-kma-47505a7dd4eb/ && make && mv kma kma_index kma_shm kma_update /usr/local/bin/

RUN apt-get update && apt-get -y install  ncbi-blast+=2.9.0-2 prodigal=1:2.6.3-4 bamtools=2.5.1+dfsg-5build1 bowtie2=2.3.5.1-6build1


COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt


RUN git clone -b rlim/fix-header-parse-space https://github.com/rzlim08/rgi.git
RUN mkdir -p rgi/inputs/heatmap_inputs
RUN cd rgi/ && python setup.py build && python setup.py test && python setup.py install

RUN mkdir -p wildcard/

RUN curl -LO https://czid-public-references.s3.us-west-2.amazonaws.com/test/AMRv2/card.json && rgi card_annotation -i card.json && curl -LO https://czid-public-references.s3.us-west-2.amazonaws.com/test/AMRv2/wildcard_database_v3.1.0.fasta && curl -L https://czid-public-references.s3.us-west-2.amazonaws.com/test/AMRv2/61_kmer_db.json -o wildcard/61_kmer_db.json && curl -L https://czid-public-references.s3.us-west-2.amazonaws.com/test/AMRv2/all_amr_61mers.txt -o wildcard/all_amr_61mers.txt && curl -L https://czid-public-references.s3.us-west-2.amazonaws.com/test/AMRv2/index-for-model-sequences.txt -o wildcard/index-for-model-sequences.txt && rgi load --wildcard_annotation wildcard_database_v3.1.0.fasta --wildcard_version 3.1.0 --wildcard_index wildcard/index-for-model-sequences.txt --card_annotation card_database_v3.2.3.fasta -i card.json --kmer_database wildcard/61_kmer_db.json --amr_kmers wildcard/all_amr_61mers.txt --kmer_size 61 && rm -r wildcard/ wildcard_database_v3.1.0.fasta

