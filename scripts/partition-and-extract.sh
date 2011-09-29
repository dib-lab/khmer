#! /bin/bash
SCRIPTPATH=`dirname $0`
python $SCRIPTPATH/do-partition.py $1 && \
python $SCRIPTPATH/extract-partitions.py $1.part && \
tail $1.part.dist
