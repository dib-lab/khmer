#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
from screed.fasta import fasta_iter

for record in fasta_iter(open(sys.argv[1]), parse_description=False):
    area = record['name'].split()[-1]
    assert area[:2] == 'r='
    area = float(area[2:])

    if area < 500:
        print '>%s\n%s' % (record['name'], record['sequence'],)
