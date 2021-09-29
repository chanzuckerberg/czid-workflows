FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list

RUN apt-get update && apt-get install -y curl build-essential libz-dev git python3-pip

LABEL maintainer="IDseq Team <idseq-tech@chanzuckerberg.com>"

WORKDIR /tmp
RUN git clone --recursive https://github.com/mlin/minimap2-scatter.git

RUN curl -Ls https://github.com/chanzuckerberg/s3parcp/releases/download/v0.2.0-alpha/s3parcp_0.2.0-alpha_Linux_x86_64.tar.gz | tar -C /usr/bin -xz s3parcp

WORKDIR /tmp/minimap2-scatter
RUN make minimap2
RUN mv /tmp/minimap2-scatter/minimap2/minimap2 /usr/local/bin/minimap2

COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY idseq_utils /tmp/idseq_utils
RUN pip3 install /tmp/idseq_utils && rm -rf /tmp/idseq_utils