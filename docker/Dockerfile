FROM debian:stable
MAINTAINER crusoe@ucdavis.edu

ENV PACKAGES python-dev zlib1g-dev libbz2-dev gcc git python-setuptools g++ make ca-certificates python-pip 

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

WORKDIR /home

RUN git clone --branch v2.0 https://github.com/dib-lab/khmer.git && \
    cd khmer && \
    make install-dep && \
    python setup.py build_ext --include-dirs third-party/seqan/core/include:third-party/smhasher:lib --libraries z,bz2 install
