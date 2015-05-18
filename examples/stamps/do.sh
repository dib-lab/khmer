#!/bin/bash
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
load-into-counting.py -x 1e8 -k 20 stamps-reads.ct \
	../../data/stamps-reads.fa.gz
abundance-dist.py stamps-reads.ct ../../data/stamps-reads.fa.gz \
	stamps-reads.hist
normalize-by-median.py -k 20 -C 10 -x 1e8 ../../data/stamps-reads.fa.gz \
	--savetable stamps-dn.ct
abundance-dist.py stamps-dn.ct stamps-reads.fa.gz.keep stamps-dn.hist
do-partition.py -k 32 -x 1e8 -s 1e4 -T 8 stamps-part \
	../../data/stamps-reads.fa.gz
../../sandbox/error-correct-pass2.py -C 10 stamps-dn.ct \
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
	--savetable stamps-dn3.ct

abundance-dist.py stamps-dn3.ct stamps-reads.fa.gz.keep.abundfilt.keep \
	stamps-dn3.hist
