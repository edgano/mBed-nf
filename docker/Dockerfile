FROM debian:jessie
MAINTAINER Edgar Garriga <docker@edgargarriga.com>

## 
## Docker image for mBed distances
##

RUN echo "deb http://ftp.us.debian.org/debian unstable main contrib non-free" >> /etc/apt/sources.list.d/unstable.list

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
                ed \
                less \
                wget \
                git \
                python \
                build-essential \
                cmake \
                curl \
                libargtable2-0 \
                ca-certificates \
                libpng12-0 \
                python-biopython \
                python-numpy \
                python-setuptools \
                default-jdk \
                libpng-dev \
                vim
                
## download the repo
RUN git clone https://github.com/edgano/mBed-nf && \

##compile SeqFile
    cd mBed-nf/bin/SeqFilter && \
    make

##compile mBed
RUN cd mBed-nf/bin/mBed && \
    make

