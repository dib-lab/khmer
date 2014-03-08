#! /bin/bash
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#

# this is a hack for my laptop, which requires a special Python env. --titus
if [ -d env ]; then
   . env/bin/activate
fi

SEQS=$1
SCRIPTPATH=`dirname $0`

BASENAME=`basename $SEQS .gz`
BASENAME=`basename $BASENAME .fa`
BASENAME=`basename $BASENAME .fasta`
BASENAME=`basename $BASENAME .fq`

K=20
HASHBITS_SIZE=1e7
N_TABLES=4
SUBSET_SIZE=1e2
KEEP_SUBSETS=    # --keep-subsets turns on
OUTPUT_GROUPS=-n # -n prevents groups from being output

$SCRIPTPATH/load-graph.py -k $K -x $HASHBITS_SIZE -N $N_TABLES $BASENAME $SEQS
$SCRIPTPATH/partition-graph.py --subset-size $SUBSET_SIZE $BASENAME
$SCRIPTPATH/merge-partitions.py -k $K $BASENAME
$SCRIPTPATH/annotate-partitions.py -k $K $BASENAME $SEQS
$SCRIPTPATH/extract-partitions.py $BASENAME `basename ${SEQS}`.part
