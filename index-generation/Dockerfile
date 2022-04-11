FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        # Build only dependencies
        g++ \
        libperl4-corelibs-perl \
        # Runtime dependencies
        make \
        # On versions higher than 1.62 compilation errors because of an overloading ambiguity with the `advance` method
        libboost1.62-all-dev \
        curl \
        git \
        wget \
        liblz4-tool \
        pigz \
        python-pandas \
        python3-gdbm \
        python3-pip \
        software-properties-common

# Install s3parcp
RUN curl -Ls https://github.com/chanzuckerberg/s3parcp/releases/download/v0.2.0-alpha/s3parcp_0.2.0-alpha_Linux_x86_64.tar.gz | tar -C /usr/bin -xz s3parcp

# install minimap2 and seqkit
RUN git clone --single-branch --branch distributed-mapping  https://github.com/mlin/minimap2.git \
&& cd minimap2 && make && mv ./minimap2 /usr/local/bin/
RUN curl -L https://github.com/shenwei356/seqkit/releases/download/v2.0.0/seqkit_linux_amd64.tar.gz | tar xz && mv seqkit /usr/local/bin/

# install diamond 
RUN git clone --single-branch --branch scatter-gather https://github.com/morsecodist/diamond && cmake -H. -B diamond/build diamond

# install boto 3 and awscli
RUN python3 -m pip install boto3==1.12.27 botocore==1.15.39
RUN python3 -m pip install awscli==1.18.39 --upgrade
ENV PATH="${PATH}:/usr/local/bin"

COPY *.py /usr/local/bin/

WORKDIR /workspace
