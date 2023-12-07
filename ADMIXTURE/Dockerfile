FROM --platform=linux/amd64 ubuntu:20.04
MAINTAINER Zan Koenig


RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        autoconf \
        automake \
        bzip2 \
        build-essential \
        ca-certificates \
        curl \
        g++ \
        gcc \
        gzip \
        libboost-all-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        sudo \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN cd /opt && \
    wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz && \
    tar -xvf admixture_linux-1.3.0.tar.gz && \
    mv /opt/dist/admixture_linux-1.3.0/admixture /usr/local/bin







