# test pip install of starseqr
FROM ubuntu:16.04

MAINTAINER Jeff Jasper jasper1918@gmail.com

RUN apt-get update
RUN apt-get -y install wget git tar python python-pip g++ make zlib1g-dev

RUN pip install --pre starseqr
