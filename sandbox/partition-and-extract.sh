#! /bin/bash
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
SCRIPTPATH=`dirname $0`
python $SCRIPTPATH/do-partition.py $1 && \
python $SCRIPTPATH/extract-partitions.py $1.part && \
tail $1.part.dist
