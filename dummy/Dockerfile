FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
RUN sed -i s/archive.ubuntu.com/us-west-2.ec2.archive.ubuntu.com/ /etc/apt/sources.list

RUN apt-get update && apt-get -y install python3-pip

RUN pip3 install mysql-connector-python boto3

COPY dummy.py /usr/local/bin/ 
