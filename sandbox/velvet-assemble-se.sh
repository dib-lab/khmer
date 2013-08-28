#! /bin/bash
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
scriptpath=$1
filename=$2
K=$3

BASE=`basename $filename`

if [ \! -f $BASE.strip \]; then
   python $scriptpath/strip-partition.py $filename > $BASE.strip
fi

velveth $BASE.ass.$K.single $K -fasta -short ${BASE}.strip && \
velvetg $BASE.ass.$K.single -read_trkg yes
