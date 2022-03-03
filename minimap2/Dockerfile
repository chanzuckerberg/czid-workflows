FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list

RUN apt-get update && apt-get install -y build-essential libz-dev git python3-pip

WORKDIR /tmp
RUN git clone --recursive https://github.com/mlin/minimap2-scatter.git

WORKDIR /tmp/minimap2-scatter
RUN make minimap2
RUN mv /tmp/minimap2-scatter/minimap2/minimap2 /usr/local/bin/minimap2