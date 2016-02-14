FROM debian:8.3
MAINTAINER Kevin Murray <spam@kdmurray.id.au>

RUN apt-get update && \
    apt-get -yy upgrade && \
    apt-get -yy install build-essential \
                    zlib1g-dev \
                    cmake \
                    git \
                    libncurses5-dev \
                    python-dev \
                    python3-dev \
                    python-pip \
                    python3-pip \
                    && \
    apt-get clean && \
    apt-get autoclean && \
    apt-get -yy autoremove && \
    pip install screed==0.9 docopt==0.6.2 numpy==1.10.2 && \
    cd /usr/local/src && \
    git clone https://github.com/kdmurray91/libqes && \
    cd libqes && \
    cmake . && \
    make && make test && make install && \
    cd .. && \
    git clone --recursive https://github.com/nh13/DWGSIM/ dwgsim && \
    cd dwgsim && \
    make && \
    mv dwgsim /usr/local/bin && \
    rm -rf /usr/local/src/*
