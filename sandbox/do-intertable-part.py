#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import khmer

infile = sys.argv[1]
khmer.do_intersection_partition(19, 250000013,  # 1000000007,
                                infile, infile + '.ipart')
