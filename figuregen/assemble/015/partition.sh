#!/bin/sh

python ../graph-size.py 31 1776873010 3 1 500 ../../data/msb2.fq msb2.fa
python ../load-graph.py -k 31 --n_hashes 3 --hashsize 1776873010 msb2 msb2.fa
python ../partition-graph.py --threads 1 msb2
python ../merge-partitions.py -k 31 msb2
python ../annotate-partitions.py -k 31 msb2 msb2.fa
python ../extract-partitions.py msb2 msb2.fa.part
./assemble.sh
