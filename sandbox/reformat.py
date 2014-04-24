#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import sys
import screed

for record in screed.open(sys.argv[1], parse_description=False):
    name, descr = record.name.split(" ", 1)
    name = name.lstrip('>')
    name += "\t".join(descr.split())

    print '>%s\n%s' % (name, record.sequence,)
