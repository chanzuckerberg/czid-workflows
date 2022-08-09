FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
ARG MINIWDL_VERSION=1.1.5

LABEL maintainer="CZ ID Team <idseq-tech@chanzuckerberg.com>"

RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list; \
    echo 'APT::Install-Recommends "false";' > /etc/apt/apt.conf.d/98czid; \
    echo 'APT::Install-Suggests "false";' > /etc/apt/apt.conf.d/99czid

RUN apt-get -q update
RUN apt-get -q install -y curl python3-pip git build-essential libz-dev

# Match versions from Ubuntu 20.04 LTS (please do not upgrade except to fix bugs)
RUN pip3 install boto3==1.9.253 marisa-trie==0.7.7 requests biopython

ADD https://raw.githubusercontent.com/chanzuckerberg/miniwdl/v${MINIWDL_VERSION}/examples/clean_download_cache.sh /usr/local/bin
RUN chmod +x /usr/local/bin/clean_download_cache.sh

RUN curl -Ls https://github.com/chanzuckerberg/s3parcp/releases/download/v0.2.0-alpha/s3parcp_0.2.0-alpha_Linux_x86_64.tar.gz | tar -C /usr/bin -xz s3parcp

WORKDIR /tmp
RUN git clone --recursive https://github.com/mlin/minimap2-scatter.git
WORKDIR /tmp/minimap2-scatter
RUN make minimap2
RUN mv /tmp/minimap2-scatter/minimap2/minimap2 /usr/local/bin/minimap2

RUN curl -Ls https://github.com/chanzuckerberg/s3parcp/releases/download/v0.2.0-alpha/s3parcp_0.2.0-alpha_Linux_x86_64.tar.gz | tar -C /usr/bin -xz s3parcp

COPY idseq-dag /tmp/idseq-dag
RUN pip3 install /tmp/idseq-dag && rm -rf /tmp/idseq-dag

COPY idseq_utils /tmp/idseq_utils
RUN pip3 install /tmp/idseq_utils && rm -rf /tmp/idseq_utils
