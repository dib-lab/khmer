#! /bin/bash
filename=$1
K=$2

scriptpath=`dirname $0`
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
