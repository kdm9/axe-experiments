#/bin/bash
set -xe

AXE_VERSION=0.3.2
MASON_VERSION=2.0.5
SEQAN_VERSION=2.1.1
FLEXBAR_VERSION=2.7.0

################################################################################
#                                     APT                                      #
################################################################################

apt-get update
apt-get -yy upgrade
apt-get -yy install build-essential \
                cmake               \
                libbz2-dev          \
                libtbb-dev          \
                python              \
                python3-dev         \
                python3-docopt      \
                python3-numpy       \
                python3-pip         \
                wget                \
                xz-utils            \
                zlib1g-dev
apt-get clean
apt-get autoclean
apt-get -yy autoremove
rm -rf /var/lib/apt/lists/*

pip3 install screed==0.9 snakemake==3.8.2

cd /usr/local/src


################################################################################
#                                     AXE                                      #
################################################################################

wget -O axe-${AXE_VERSION}.tar.gz \
    https://github.com/kdmurray91/axe/archive/${AXE_VERSION}.tar.gz

tar xvf axe-${AXE_VERSION}.tar.gz

pushd axe-${AXE_VERSION}
cmake .
make && make test && make install
popd
rm -rf /usr/local/src/*


################################################################################
#                                    MASON                                     #
################################################################################

wget http://packages.seqan.de/mason2/mason2-${MASON_VERSION}-Linux-x86_64.tar.xz
tar xvf mason2*.tar*
mv mason2-*-Linux-x86_64/bin/* /usr/local/bin
rm -rf /usr/local/src/*


################################################################################
#                                   FLEXBAR                                    #
################################################################################

wget https://github.com/seqan/seqan/releases/download/seqan-v${SEQAN_VERSION}/seqan-library-${SEQAN_VERSION}.tar.xz
tar xvf seqan-library*.tar*

wget -O flexbar-${FLEXBAR_VERSION}.tar.gz \
    https://github.com/seqan/flexbar/archive/v${FLEXBAR_VERSION}.tar.gz
tar xvf flexbar-${FLEXBAR_VERSION}.tar*

mv seqan-library*/include flexbar-${FLEXBAR_VERSION}/

pushd flexbar-${FLEXBAR_VERSION}
cmake .
make && mv flexbar /usr/local/bin/flexbar
popd
rm -rf /usr/local/src/*
