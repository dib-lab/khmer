#! /bin/bash
scriptpath=$1
filename=$2
K=$3

BASE=`basename $filename`

if [ \! -f $BASE.strip \]; then
   python $scriptpath/strip-partition.py $filename > $BASE.strip
fi

echo running velveth && \
velveth $BASE.ass.$K.oases $K -fasta -short ${BASE}.strip && \
echo running velvetg && \
velvetg $BASE.ass.$K.oases -read_trkg yes  && \
echo running oases && \
oases $BASE.ass.$K.oases
