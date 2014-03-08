#! /bin/bash
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
filename=$1
K=$2
scriptpath=`dirname $0`

BASE=`basename $filename`

if [ \! -f $BASE.se -o \! -f $BASE.pe \]; then
   python $scriptpath/strip-and-split-for-assembly.py $filename $BASE
fi

if [ \! -s $BASE.se -o \! -s $BASE.pe ]; then
   echo 'WARNING -- one or more sequences files are EMPTY; may fail'
fi

velveth $BASE.ass.$K $K -fasta -short ${BASE}.se -shortPaired ${BASE}.pe && \
velvetg $BASE.ass.$K -read_trkg yes -exp_cov auto -cov_cutoff auto -scaffolding no

rm ${BASE}.ass.$K/{Graph2,LastGraph,PreGraph,Roadmaps,Sequences}
