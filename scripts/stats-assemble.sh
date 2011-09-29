#! /bin/bash
scriptpath=$1
filename=$2


BASE=`basename $filename`

python $scriptpath/assemstats3.py 500 $BASE/contigs.fa > $BASE.stat 
