FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="CZ ID Team <idseq-tech@chanzuckerberg.com>"
LABEL description = "Image for consensus genome by metagenomic sequencing with spiked primer enrichment or amplicon sequencing"

RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list; \
  echo 'APT::Install-Recommends "false";' > /etc/apt/apt.conf.d/98czid; \
  echo 'APT::Install-Suggests "false";' > /etc/apt/apt.conf.d/99czid

# script to raise error as json
RUN echo  'set +x \nexport error=$1 cause=$2 \njq -nc ".wdl_error_message=true | .error=env.error | .cause=env.cause" > /dev/stderr \nexit 1' >> /usr/local/bin/raise_error; \
  chmod +x /usr/local/bin/raise_error

RUN apt-get -qq update && apt-get -qq -y install \
  jq \
  moreutils \
  curl \
  locales \
  zip \
  unzip \
  hmmer \
  httpie \
  zlib1g-dev \
  libfreetype6-dev \
  libhts-dev \
  pkg-config \
  apt-utils \
  libinline-c-perl \
  libwww-perl \
  bcftools \
  freebayes \
  kraken2 \
  minimap2 \
  samtools \
  seqtk \
  trim-galore \
  python3-dev \
  python3-pip \
  python3-setuptools \
  python3-wheel \
  python3-yaml \
  python3-dateutil \
  python3-biopython \
  python3-pysam \
  python3-seaborn \
  seqtk \
  bedtools \
  build-essential \
  automake \
  tabix \
  fasta3 \
  wget \
  && locale-gen en_US.UTF-8

# These newer versions of infernal and ncbi-blast+ are required by VADR
RUN curl -O -L https://ftp.osuosl.org/pub/ubuntu/pool/universe/i/infernal/infernal_1.1.4-1_amd64.deb && \
  dpkg -i infernal_1.1.4-1_amd64.deb && \
  curl -O -L https://ftp.osuosl.org/pub/ubuntu/pool/universe/n/ncbi-blast+/ncbi-blast+_2.10.1+ds-1_amd64.deb && \
  dpkg -i ncbi-blast+_2.10.1+ds-1_amd64.deb

# See https://github.com/ablab/quast/issues/157
RUN pip3 install multiqc==1.8 quast==5.0.2 && \
  sed -i 's/cgi/html/' /usr/local/lib/python3.8/dist-packages/quast_libs/site_packages/jsontemplate/jsontemplate.py

RUN /bin/bash -c "set -e; mkdir ivar; pushd ivar; \
  curl -L https://github.com/andersen-lab/ivar/archive/v1.3.1.tar.gz | tar -xvz --strip-components 1; \
  ./autogen.sh; \
  ./configure; \
  make -j8; \
  make install; \
  popd; rm -rf ivar"

RUN /bin/bash -c "set -e; mkdir -p /usr/local/lib/site_perl; \
  curl -L https://github.com/nawrockie/sequip/archive/sequip-0.08.tar.gz | tar -xz --strip-components 1 -C /usr/local/lib/site_perl; \
  rm -rf Bio-Easel vadr; \
  mkdir -p Bio-Easel; pushd Bio-Easel; \
  curl -L https://github.com/nawrockie/Bio-Easel/archive/Bio-Easel-0.14.tar.gz | tar -xz --strip-components 1; \
  mkdir -p src/easel; pushd src/easel; \
  curl -L https://github.com/EddyRivasLab/easel/archive/easel-0.48.tar.gz | tar -xz --strip-components 1; \
  autoconf; \
  ./configure --enable-pic --enable-sse4 --enable-avx512 --enable-threads; \
  make -j8; \
  make install; \
  popd; \
  perl Makefile.PL; \
  make -j8; \
  make install; \
  popd; \
  rm -rf Bio-Easel; \
  mkdir -p vadr; \
  curl -L https://github.com/ncbi/vadr/archive/vadr-1.2.tar.gz | tar -xz --strip-components 1 -C vadr; \
  sed -i -e 's|/scripts/esl-ssplit.pl|/esl-ssplit.pl|' vadr/v-annotate.pl; \
  mv vadr/*.pl /usr/local/bin/; \
  mv vadr/*.pm /usr/local/lib/site_perl/; \
  rm -rf vadr; \
  echo export VADRSCRIPTSDIR=/usr/local/bin \
  VADRMODELDIR=/usr/local/share/vadr/models \
  VADRBLASTDIR=/usr/bin \
  VADRFASTADIR=/usr/bin \
  VADRINFERNALDIR=/usr/bin \
  VADRHMMERDIR=/usr/bin \
  VADREASELDIR=/usr/local/bin \
  VADRBIOEASELDIR=/usr/local/bin > /etc/profile.d/vadr.sh"

RUN wget https://drive5.com/muscle/muscle_src_3.8.1551.tar.gz \
&& mkdir muscle && tar zxvf muscle_src_3.8.1551.tar.gz -C muscle && cd muscle && sed -i 's/Str\[16\]/Str[32]/' progress.cpp && make && mv muscle /usr/local/bin/


# ARTIC, Medaka and dependencies
RUN apt-get install -y python3-cffi python3-h5py python3-intervaltree python3-edlib git
RUN pip3 install ont-fast5-api parasail mappy pyspoa protobuf==3.19.* tensorflow==2.7.1 https://github.com/artic-network/fieldbioinformatics/archive/1.2.1.tar.gz
RUN pip3 install medaka==1.6.0 --no-deps
RUN pip3 install marisa-trie
RUN pip3 install boto3
RUN pip3 install git+https://github.com/rzlim08/PyVCF.git


RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.0.0/seqkit_linux_amd64.tar.gz && tar zxvf seqkit_linux_amd64.tar.gz \
&& mv seqkit /usr/local/bin/ && rm seqkit_linux_amd64.tar.gz

