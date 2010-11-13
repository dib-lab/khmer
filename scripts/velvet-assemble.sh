#! /bin/bash
filename=$1
K=$2

python ~/khmer/scripts/strip-partition.py $filename | \
    velveth $filename.ass.$K $K -fasta -short - && \
velvetg $filename.ass.$K -read_trkg yes -exp_cov 3 -cov_cutoff 0
