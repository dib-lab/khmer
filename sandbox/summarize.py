#! /usr/bin/env python
import screed
import sys

fp = open('summary.txt', 'w')

total_bp = 0
total_seqs = 0
for m, filename in enumerate(sys.argv[1:]):
    print 'opening', filename
    file_bp = 0

    for n, record in enumerate(screed.open(filename)):
        if n % 100000 == 0 and n > 0:
            print '... %d, in %s -- %d of %d files' % \
                (n, filename, m, len(sys.argv) - 1)
        file_bp += len(record.sequence)

    file_seqs = n

    print >>fp, '%d %d %s' % (file_seqs, file_bp, filename)

    total_bp += file_bp
    total_seqs += file_seqs

print >>fp, '%d %d %s' % (total_seqs, total_bp, 'total')

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
