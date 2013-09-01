#! /bin/bash
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
K=33
if [ "$2" ]; then
   K=$2
fi
echo assembling $1 with K = $K

python ~/dev/khmer/scripts/strip-partition.py $1 | velveth $1.ass.$K $K -short -
velvetg $1.ass.$K -read_trkg yes
#-scaffolding no
