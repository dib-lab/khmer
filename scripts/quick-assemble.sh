#! /bin/bash
K=33
if [ "$2" ]; then
   K=$2
fi
echo assembling $1 with K = $K

python ~/dev/khmer/scripts/strip-partition.py $1 | velveth $1.ass.$K $K -short -
velvetg $1.ass.$K -exp_cov auto -read_trkg yes -cov_cutoff 0
