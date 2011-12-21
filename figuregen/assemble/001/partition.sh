#!/bin/sh

python ../graph-size.py 31 1848546972 7 1 200 ../../data/msb2.fq msb2.fa
python ../load-graph.py -k 31 --n_hashes 7 --hashsize 1848546972 msb2 msb2.fa
python ../partition-graph.py --threads 1 msb2
python ../merge-partitions.py -k 31 msb2
python ../annotate-partitions.py -k 31 msb2 msb2.fa
python ../extract-partitions.py msb2 msb2.fa.part

./assemble.sh
