# syntax=docker/dockerfile:1.4
FROM ubuntu:20.04

LABEL maintainer="CZ ID Team <idseq-tech@chanzuckerberg.com>"

RUN apt-get update && apt-get -y install curl git python3 python3-pip make zip
RUN ln -s /usr/bin/python3 /usr/bin/python

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
