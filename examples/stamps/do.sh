#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
python ../../scripts/load-into-counting.py -x 1e8 -k 20 stamps-reads.kh ../../data/stamps-reads.fa.gz 
python ../../scripts/abundance-dist.py stamps-reads.kh ../../data/stamps-reads.fa.gz stamps-reads.hist
python ../../scripts/normalize-by-median.py -k 20 -C 10 -x 1e8 ../../data/stamps-reads.fa.gz --savehash stamps-dn.kh
python ../../scripts/abundance-dist.py stamps-dn.kh stamps-reads.fa.gz.keep stamps-dn.hist
python ../../scripts/do-partition.py -k 32 -x 1e8 -s 1e4 -T 8 stamps-part ../../data/stamps-reads.fa.gz 
python ../../sandbox/error-correct-pass2.py -C 10 stamps-dn.kh ../../data/stamps-reads.fa.gz 
python ../../scripts/load-into-counting.py -x 1e8 -k 20 stamps-corr.kh stamps-reads.fa.gz.corr
python ../../scripts/abundance-dist.py stamps-corr.kh stamps-reads.fa.gz.corr stamps-corr.hist
python ../../scripts/extract-partitions.py stamps-part stamps-reads.fa.gz.part
python ../../scripts/extract-partitions.py -X 1 stamps-part stamps-reads.fa.gz.part
python ../../scripts/load-into-counting.py -x 1e8 -k 20 stamps-part.g0.kh stamps-part.group0000.fa 
python ../../scripts/load-into-counting.py -x 1e8 -k 20 stamps-part.g1.kh stamps-part.group0001.fa 
python ../../scripts/abundance-dist.py stamps-part.g0.kh stamps-part.group0000.fa stamps-part.g0.hist
python ../../scripts/abundance-dist.py stamps-part.g1.kh stamps-part.group0001.fa stamps-part.g1.hist

python ../../scripts/filter-abund.py stamps-dn.kh stamps-reads.fa.gz.keep
python ../../scripts/normalize-by-median.py -x 1e8 -k 20 -C 10 stamps-reads.fa.gz.keep.abundfilt --savehash stamps-dn3.kh

python ../../scripts/abundance-dist.py stamps-dn3.kh stamps-reads.fa.gz.keep.abundfilt.keep stamps-dn3.hist
