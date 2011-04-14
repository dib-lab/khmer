#! /bin/bash
scriptpath=$1
filename=$2
K=$3

BASE=`basename $filename`

if [ \! -f $BASE.strip \]; then
   python $scriptpath/strip-partition.py $filename > $BASE.strip
fi

velveth $BASE.ass.$K.single $K -fasta -short ${BASE}.strip && \
velvetg $BASE.ass.$K.single -read_trkg yes
