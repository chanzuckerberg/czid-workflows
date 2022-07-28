FROM ubuntu:20.04

LABEL maintainer="CZ ID Team <idseq-tech@chanzuckerberg.com>"

SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get -y install curl git python3 python3-pip

RUN pip3 install pandas

RUN curl -sO https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
    bash Miniconda2-latest-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda2-latest-Linux-x86_64.sh  && \
    echo 'PATH=/usr/local/miniconda/bin:$PATH' >> /etc/bash.bashrc

ENV PATH=/usr/local/miniconda/bin:$PATH

RUN conda init bash 

RUN git clone -b rlim/parse-kma https://github.com/rzlim08/rgi.git && conda env create -f rgi/conda_env.yml

RUN mkdir -p rgi/inputs/heatmap_inputs
RUN . /usr/local/miniconda/etc/profile.d/conda.sh && conda activate rgi && cd rgi/ && python setup.py build && python setup.py test && python setup.py install

#RUN mkdir /opt/CARD_DBs/ && cd /opt/CARD_DBs && curl -L https://card.mcmaster.ca/latest/data  | tar -xvj
