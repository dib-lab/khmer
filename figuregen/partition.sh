#!/bin/sh
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#

python simgenome.py genome.fa 10000000 1000 seq

for p in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15
do
   mkdir $p
   cp do-partition.py $p
   cp genome.fa $p
   cd $p
   python do-partition.py genome.fa $p
   cd ..
done
