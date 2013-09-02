#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import sys
from screed.fasta import fasta_iter

n = 0
for filename in sys.argv[1:]:
    sys.stderr.write('... %s %d\n' % (filename, n))
    idx = filename.find('group')
    assert idx != -1, filename

    group_num = int(filename[idx + 5:].split('.')[0])

    for record in fasta_iter(open(filename + '/contigs.fa')):
        print n, group_num, len(record['sequence'])
        n += 1
