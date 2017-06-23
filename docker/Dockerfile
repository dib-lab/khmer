FROM debian:stable
MAINTAINER khmer-project@idyll.org

ENV PACKAGES python-dev zlib1g-dev libbz2-dev gcc git python-setuptools g++ \
             make ca-certificates python-pip python-wheel python-enchant

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

WORKDIR /home

ARG branch
ENV branch ${branch:-v2.0}

RUN git clone https://github.com/dib-lab/khmer.git && \
    cd khmer && \
    git checkout ${branch} && \
    make install-dep && \
    python setup.py build_ext --include-dirs third-party/seqan/core/include:third-party/smhasher:lib --libraries z,bz2 install
