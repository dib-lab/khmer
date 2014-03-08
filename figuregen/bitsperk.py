#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import math

for p in [0.001, 0.01, 0.05, 0.1, 0.15, 0.20]:
   print p, 1 / math.log(2, 2) * math.log(1 / p, 2)
