FROM debian:stable
MAINTAINER khmer-project@idyll.org

ENV PACKAGES python3-dev zlib1g-dev libbz2-dev gcc git python3-setuptools g++ \
             make ca-certificates python3-pip python3-wheel python3-enchant

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

WORKDIR /home

ARG branch
ENV branch ${branch:-v2.1.2}
ARG slug
ENV slug ${slug:-"dib-lab/khmer"}

RUN git clone https://github.com/dib-lab/khmer.git && \
    cd khmer && \
    git checkout ${branch} && \
    python3 setup.py build_ext --include-dirs include:third-party/seqan/core/include:third-party/smhasher:third-party/cqf:third-party/rollinghash --libraries z,bz2 install
