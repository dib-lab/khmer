#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import screed
import sys

for record in screed.open(sys.argv[1]):
    name = record['name']
    sequence = record['sequence']

    name = name.split()[0]

    print '>%s\n%s' % (name, sequence,)
