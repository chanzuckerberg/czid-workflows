FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="IDseq Team idseq-tech@chanzuckerberg.com"
LABEL description = "Image for consensus genome by metagenomic sequencing with spiked primer enrichment or amplicon sequencing"

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
  hmmer \
  httpie \
  zlib1g-dev \
  libfreetype6-dev \
  pkg-config \
  apt-utils \
  ncbi-blast+ \
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
  autoconf \
  && locale-gen en_US.UTF-8

RUN curl -O -L http://ftp.osuosl.org/pub/ubuntu/pool/universe/i/ivar/ivar_1.2.2+dfsg-3_amd64.deb && \
  dpkg -i ivar_1.2.2+dfsg-3_amd64.deb && \
 curl -O -L http://ftp.osuosl.org/pub/ubuntu/pool/universe/i/infernal/infernal_1.1.4-1_amd64.deb && \
 dpkg -i infernal_1.1.4-1_amd64.deb

RUN pip3 install multiqc==1.8 quast==5.0.2

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
  curl -L https://github.com/ncbi/vadr/archive/vadr-1.1.3.tar.gz | tar -xz --strip-components 1 -C vadr; \
  sed -i -e 's|/scripts/esl-ssplit.pl|/esl-ssplit.pl|' vadr/v-annotate.pl; \
  mv vadr/*.pl /usr/local/bin/; \
  mv vadr/*.pm /usr/local/lib/site_perl/; \
  rm -rf vadr; \
  echo export VADRSCRIPTSDIR=/usr/local/bin \
  VADRMODELDIR=/usr/local/share/vadr/models \
  VADRBLASTDIR=/usr/bin \
  VADRINFERNALDIR=/usr/bin \
  VADRHMMERDIR=/usr/bin \
  VADREASELDIR=/usr/local/bin \
  VADRBIOEASELDIR=/usr/local/bin > /etc/profile.d/vadr.sh"

RUN apt-get remove -y autoconf build-essential python3-dev && apt-get autoremove -y