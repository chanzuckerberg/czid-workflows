# Set the base image
FROM ubuntu:16.04
ENV DEBIAN_FRONTEND noninteractive

# File Author / Maintainer
MAINTAINER IDseq Team idseq-tech@chanzuckerberg.com

# Add packages, update image, and clear cache
RUN apt-get update && apt-get install -y build-essential curl wget python-pip python-dev python-scipy python-redis gdebi-core zip unzip g++ zlib1g-dev gcc pkg-config apt-utils make perl cmake libbz2-dev liblzma-dev

RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update && apt-get install -y python3.7

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.5 1
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 2

RUN pip install --upgrade pip
RUN pip install redis biopython pysam
RUN pip install htseq==0.6.1p1


WORKDIR /tmp
# install STAR
RUN curl -L https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz | tar xz
RUN mv STAR-2.5.3a/bin/Linux_x86_64_static/* /usr/local/bin
RUN rm -rf STAR-2.5.3a

# Compile and install bowtie
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-source.zip
RUN unzip bowtie2-2.3.2-source.zip
WORKDIR /tmp/bowtie2-2.3.2
RUN make -j 16 NO_TBB=1 # we may want tbb (for "superior thread scaling") later, so will need to find a way to compile without "NO_TBB=1"
RUN mv -t /usr/local/bin/ bowtie2 bowtie2-align-l bowtie2-align-s bowtie2-build bowtie2-build-l bowtie2-build-s bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-s

WORKDIR /tmp
# install samtools
RUN apt-get install -y libncurses-dev libbz2-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2
RUN tar -jxf samtools-1.5.tar.bz2
WORKDIR /tmp/samtools-1.5
RUN ./configure --disable-lzma
RUN make -j 16
RUN mv samtools /usr/local/bin/

WORKDIR /tmp
# Compile and install PriceSeqFilter
RUN wget http://derisilab.ucsf.edu/software/price/PriceSource140408.tar.gz
RUN tar -xzf PriceSource140408.tar.gz
WORKDIR /tmp/PriceSource140408
RUN make -j 16
RUN mv PriceSeqFilter /usr/local/bin/

WORKDIR /tmp
# Compile and install cdhit-dup tools
RUN wget https://github.com/weizhongli/cdhit/archive/V4.6.8.zip
RUN unzip V4.6.8.zip
WORKDIR /tmp/cdhit-4.6.8
RUN make -j 16
WORKDIR /tmp/cdhit-4.6.8/cd-hit-auxtools
RUN make -j 16
RUN mv cd-hit-dup /usr/local/bin/

WORKDIR /tmp
# Compile and install Fastax tools
#RUN wget http://launchpadlibrarian.net/161878011/libgtextutils0_0.7-1_amd64.deb
#RUN wget http://launchpadlibrarian.net/162265652/fastx-toolkit_0.0.14-1_amd64.deb
#RUN gdebi --non-interactive libgtextutils0_0.7-1_amd64.deb
#RUN gdebi --non-interactive fastx-toolkit_0.0.14-1_amd64.deb

# For aegea
RUN apt-get install -y python3-pip cython python3.7-dev
RUN pip3 install awscli-cwlogs==1.4.0 keymaker==0.2.1 boto3==1.12.39 awscli==1.11.44 dynamoq==0.0.5 tractorbeam==0.1.3
# Avoid cython error in 0.15.4
RUN pip3 install pysam==0.15.3 biopython
#RUN echo iptables-persistent iptables-persistent/autosave_v4 boolean true | debconf-set-selections
RUN apt-get update && apt-get install -y iptables-persistent debian-goodies bridge-utils pixz cryptsetup-bin mdadm btrfs-tools libffi-dev libssl-dev libxml2-dev libxslt1-dev libyaml-dev libcurl4-openssl-dev libjemalloc-dev libzip-dev libsnappy-dev liblz4-dev libgmp-dev libmpfr-dev libhts-dev libncurses5-dev htop pydf jq httpie python-dev python-cffi python-pip python-setuptools python-wheel python-virtualenv python-requests python-yaml python3-dev python3-cffi python3-pip python3-setuptools python3-wheel python3-requests python3-yaml nfs-common unzip build-essential cmake libtool autoconf ruby sysstat dstat numactl gdebi-core stunnel moreutils curl wget git aria2 sift
# Strange that we have to do this, but if we don't, aegea tries to do it, and it fails then with some urllib3 bug.
RUN apt-get -y install awscli

RUN apt-get install -y bsdtar alien

# For de-novo assembly
WORKDIR /tmp/spades_build
RUN git clone https://github.com/ablab/spades.git
WORKDIR /tmp/spades_build/spades
RUN git checkout spades_3.11.0
WORKDIR /tmp/spades_build/spades/assembler
RUN PREFIX=/usr/local ./spades_compile.sh
RUN /usr/local/bin/spades.py --test

# For nonhost fastq filtering
WORKDIR /tmp/seqtk_build
RUN git clone https://github.com/lh3/seqtk.git
WORKDIR /tmp/seqtk_build/seqtk
RUN make && make install

WORKDIR /tmp

# Compile and install gmap/gsnap
RUN wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-11-15.tar.gz
RUN mkdir gmap-gsnap && tar xf gmap-gsnap-2017-11-15.tar.gz -C gmap-gsnap --strip-components 1
WORKDIR /tmp/gmap-gsnap
RUN ./configure --prefix=/usr/local
RUN make -j 16 && make check && make install
RUN rm -rf /tmp/gmap-gsnap /tmp/gmap-gsnap-2017-11-15.tar.gz
RUN gsnapl --version

# For srst2, install forked srst2 python 3 compatible repo
RUN pip install scipy
RUN pip install git+https://github.com/chanzuckerberg/srst2
# TODO: Test both pip installations, consider keeping pip use consistent
RUN pip3 install pandas
RUN pip3 install psutil==5.6.3
# Bedtools for obtaining total reads for each gene and calculating reads per million
RUN apt-get install -y bedtools

# Blast command line
RUN wget -N ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-1.x86_64.rpm
RUN alien -i ncbi-blast-2.6.0+-1.x86_64.rpm
RUN rm -rf ncbi-blast-2.6.0+-1.x86_64.rpm

# For adapter trimming
RUN apt install python-cutadapt
WORKDIR /tmp
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
RUN unzip Trimmomatic-0.38.zip
RUN mv Trimmomatic-0.38/trimmomatic-0.38.jar /usr/local/bin/
RUN apt-get update && apt-get install -y default-jre

# For phylogenetic trees
WORKDIR /tmp
RUN wget https://sourceforge.net/projects/ksnp/files/kSNP3.1_Linux_package.zip
RUN unzip -o kSNP3.1_Linux_package.zip
WORKDIR /tmp/kSNP3.1_Linux_package/kSNP3
RUN cp -r * /usr/local/bin/
RUN sed -i 's:set kSNP=/usr/local/kSNP3:set kSNP=/usr/local/bin:g' /usr/local/bin/kSNP3 # edit the kSNP3 executable so it looks for other executables in /usr/local/bin
RUN apt install tcsh
RUN kSNP3

RUN apt-get -y install liblz4-tool
RUN apt-get -y install lbzip2

# Picard for average fragment size https://github.com/broadinstitute/picard
WORKDIR /tmp
# r-base is a dependency of collecting input size metrics https://github.com/bioconda/bioconda-recipes/pull/16398
RUN apt-get -y install r-base
RUN wget https://github.com/broadinstitute/picard/releases/download/2.21.2/picard.jar
RUN mv picard.jar /usr/local/bin/
# Create a single executable so we can use SingleCommand
RUN printf '#!/bin/bash\njava -jar /usr/local/bin/picard.jar "$@"\n' > /usr/local/bin/picard
RUN chmod +x /usr/local/bin/picard

# Install s3parcp
WORKDIR /tmp
RUN curl -L https://github.com/chanzuckerberg/s3parcp/releases/download/v0.1.2-alpha/s3parcp_0.1.2-alpha_Linux_x86_64.tar.gz | tar zx
RUN mv s3parcp /usr/local/bin

# Cleanup
RUN rm -rf /tmp/*

WORKDIR /
