#! /usr/bin/env python
#
# This script is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
import screed
import os.path

for filename in sys.argv[1:]:
    assert 'R1' in filename
    filename2 = filename.replace('R1', 'R2')

    r1 = iter(screed.open(filename)).next()
    r2 = iter(screed.open(filename2)).next()

    assert r1.name == r2.name, (r1.name, r2.name)

    final = filename.replace('R1', '')
    print 'python /root/khmer/sandbox/interleave.py %s %s | gzip -9c > %s' % (
        filename, filename2, final)

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
