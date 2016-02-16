FROM debian:8.3
MAINTAINER Kevin Murray <spam@kdmurray.id.au>

RUN apt-get update && \
    apt-get -yy upgrade && \
    apt-get -yy install build-essential \
                    zlib1g-dev \
                    cmake \
                    git \
                    python3-dev \
                    python3-pip \
                    python3-docopt \
                    python3-numpy \
                    && \
    apt-get clean && \
    apt-get autoclean && \
    apt-get -yy autoremove

RUN pip3 install screed==0.9 \
                 snakemake==3.5.5

RUN cd /usr/local/src && \
    git clone https://github.com/kdmurray91/libqes && \
    cd libqes && \
    git checkout 0.1.21 && \
    cmake . && \
    make && make test && make install && \
    rm -rf /usr/local/src/*

ADD http://packages.seqan.de/mason2/mason2-2.0.1-Linux-x86_64.tar.bz2 /usr/local/src
RUN cd /usr/local/src && \
    tar xvf mason2-2.0.1-Linux-x86_64.tar.bz2 && \
    mv mason2-2.0.1-Linux-x86_64/bin/* /usr/local/bin && \
    rm -rf /usr/local/src/*

ADD keyfiles /
ADD Snakefile /
