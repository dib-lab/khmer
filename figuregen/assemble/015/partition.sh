#!/bin/sh

memusg python ../graph-size.py 31 1776873010 3 1 400 ../../data/msb2.fq msb2.fa
memusg python ../graph-size.py 31 1776873010 3 1 200 msb2.fa msb2a.fa
rm msb2.fa
mv msb2a.fa msb2.fa
memusg python ../graph-size.py 31 1776873010 3 1 200 msb2.fa msb2a.fa
rm msb2.fa
mv msb2a.fa msb2.fa
memusg python ../load-graph.py -k 31 --n_hashes 3 --hashsize 1776873010 msb2 msb2.fa
memusg python ../partition-graph.py --threads 1 msb2
memusg python ../merge-partitions.py -k 31 msb2
memusg python ../annotate-partitions.py -k 31 msb2 msb2.fa
memusg python ../extract-partitions.py msb2 msb2.fa.part
memusg ./assemble.sh
