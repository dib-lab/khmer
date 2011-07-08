#! /bin/bash
scriptpath=$1
filename=$2
K=$3

BASE=`basename $filename`

if [ \! -f $BASE.se -o \! -f $BASE.pe \]; then
   python $scriptpath/strip-and-split-for-assembly.py $filename $BASE
fi

velveth $BASE.ass.$K $K -fasta -short ${BASE}.se -shortPaired ${BASE}.pe && \
velvetg $BASE.ass.$K -read_trkg yes -exp_cov auto -cov_cutoff 0

rm ${BASE}.ass.$K/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
