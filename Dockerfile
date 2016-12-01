FROM ubuntu:16.04
MAINTAINER Kevin Murray

ADD . /experiments
WORKDIR /experiments
RUN bash /experiments/docker/build.sh
CMD snakemake
