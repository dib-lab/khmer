#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# *obligatory line change*
import sys
import screed

total_bp = 0
total_seqs = 0

output = []
for filename in sys.argv[1:]:
    bp = 0
    seqs = 0
    for record in screed.open(filename):
        if seqs % 100000 == 0:
            print >>sys.stderr, '...', filename, seqs
        bp += len(record.sequence)
        seqs += 1

    if seqs == 0:
        print >>sys.stderr, 'No sequences found in %s' % filename
    else:
        avg_len = bp / float(seqs)
        s = '%d bp / %d seqs; %.1f average length -- %s' % (bp,
                                                            seqs,
                                                            avg_len,
                                                            filename)
        print >>sys.stderr, '... found', s
        output.append(s)

        total_bp += bp
        total_seqs += seqs

if total_seqs == 0:
    print 'No sequences found in %d files' % len(sys.argv[1:])
else:
    print '---------------'
    print "\n".join(output)
    print '---------------'
    print '%d bp / %d seqs; %.1f average length -- total'.format(
        total_bp, total_seqs, total_bp / float(total_seqs))
