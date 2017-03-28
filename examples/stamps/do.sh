#!/bin/bash
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
#
set -e # exit as soon as one command fails
set -x # echo commands before executing them
load-into-counting.py -x 1e8 -k 20 stamps-reads.ct \
	../../data/stamps-reads.fa.gz
abundance-dist.py stamps-reads.ct ../../data/stamps-reads.fa.gz \
	stamps-reads.hist
normalize-by-median.py -k 20 -C 10 -x 1e8 ../../data/stamps-reads.fa.gz \
	--savegraph stamps-dn.ct
abundance-dist.py stamps-dn.ct stamps-reads.fa.gz.keep stamps-dn.hist
do-partition.py -k 32 -x 1e8 -s 1e4 -T 8 stamps-part \
	../../data/stamps-reads.fa.gz
../../sandbox/error-correct-pass2.py --trusted-cov 10 stamps-dn.ct \
	../../data/stamps-reads.fa.gz
load-into-counting.py -x 1e8 -k 20 stamps-corr.ct stamps-reads.fa.gz.corr
abundance-dist.py stamps-corr.ct stamps-reads.fa.gz.corr stamps-corr.hist
extract-partitions.py stamps-part stamps-reads.fa.gz.part
extract-partitions.py -X 1 stamps-part stamps-reads.fa.gz.part
load-into-counting.py -x 1e8 -k 20 stamps-part.g0.ct stamps-part.group0000.fa
load-into-counting.py -x 1e8 -k 20 stamps-part.g1.ct stamps-part.group0001.fa
abundance-dist.py stamps-part.g0.ct stamps-part.group0000.fa stamps-part.g0.hist
abundance-dist.py stamps-part.g1.ct stamps-part.group0001.fa stamps-part.g1.hist

filter-abund.py stamps-dn.ct stamps-reads.fa.gz.keep
normalize-by-median.py -x 1e8 -k 20 -C 10 stamps-reads.fa.gz.keep.abundfilt \
	--savegraph stamps-dn3.ct

abundance-dist.py stamps-dn3.ct stamps-reads.fa.gz.keep.abundfilt.keep \
	stamps-dn3.hist
