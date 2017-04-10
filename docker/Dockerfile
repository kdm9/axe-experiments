FROM ubuntu:16.04
MAINTAINER Kevin Murray

ADD loadpkgs.jl /root/loadpkgs.jl
ADD build.sh /root/build_image.sh
RUN bash /root/build_image.sh
WORKDIR /experiments
