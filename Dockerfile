FROM ubuntu:15.10
MAINTAINER ecs260@luizirber.org

ENV PACKAGES python-dev zlib1g git python-setuptools g++ make ca-certificates wget lcov
ENV KHMER_VERSION pull/1324/head:feature/hypothesis
ENV TESTATTR=hypothesis

WORKDIR /home

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
	wget https://bootstrap.pypa.io/get-pip.py && \
	python get-pip.py

RUN git clone https://github.com/dib-lab/khmer.git && \
    cd khmer && \
	git fetch origin ${KHMER_VERSION} && \
	git checkout feature/hypothesis && \
    make install-dep && \
    python setup.py develop

WORKDIR /home/khmer
CMD make test
