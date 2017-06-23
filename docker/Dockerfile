FROM debian:stable
MAINTAINER khmer-project@idyll.org
ARG branch=v2.0

ENV PACKAGES python-dev zlib1g-dev libbz2-dev gcc git python-setuptools g++ make ca-certificates python-pip python-wheel 

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

WORKDIR /home

RUN git clone --branch $branch https://github.com/dib-lab/khmer.git && \
    cd khmer && \
    make install-dep && \
    python setup.py build_ext --include-dirs third-party/seqan/core/include:third-party/smhasher:lib --libraries z,bz2 install
