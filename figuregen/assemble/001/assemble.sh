#!/bin/sh
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#

for i in msb2.group*.fa
do
   ABYSS -k33 $i -o contigs.$i
done
