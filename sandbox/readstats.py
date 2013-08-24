#! /usr/bin/env python
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

    s = '%d bp / %d seqs; %.1f average length -- %s' % (bp,
                                                        seqs,
                                                        bp / float(seqs),
                                                        filename)
    print >>sys.stderr, '... found', s
    output.append(s)

    total_bp += bp
    total_seqs += seqs

print '---------------'
print "\n".join(output)
print '---------------'
print '%d bp / %d seqs; %.1f average length -- total' % (total_bp,
                                                total_seqs,
                                                total_bp / float(total_seqs))
