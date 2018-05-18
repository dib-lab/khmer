FROM debian:stable
MAINTAINER khmer-project@idyll.org

ENV PACKAGES python3-dev zlib1g-dev libbz2-dev gcc git python3-setuptools g++ \
             make ca-certificates python3-pip python3-wheel

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

WORKDIR /home

ARG branch
ENV branch ${branch:-v2.1.2}
ARG slug
ENV slug ${slug:-dib-lab/khmer}

RUN pip3 install git+https://github.com/${slug}.git@${branch}
