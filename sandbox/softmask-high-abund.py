#! /usr/bin/env python
"""
Lowercase high-copy-number k-mers.
"""
import sys
import khmer
import screed
import gzip

THRESHOLD = 5

hashfile = sys.argv[1]
filename = sys.argv[2]
output = sys.argv[3]

print 'loading ht'
ht = khmer.load_counting_hash(hashfile)
print '...done!'

K = ht.ksize()
print 'loaded ht; K is %d, n_ht is %d, size ~ %g' % (K,
                                                     len(ht.hashsizes()),
                                                     ht.hashsizes()[0])

outfp = gzip.open(output, 'w')

total = 0
total_masked = 0

for n, record in enumerate(screed.open(filename)):
    if n % 1000 == 0:
        print '...', n

    x = []
    seq = record.sequence

    total += len(seq) - K + 1

    pos = 0
    while pos < len(seq) - K + 1:
        kmer = seq[pos:pos + K]
        if 'N' in kmer.upper():
            x.extend(kmer)
            pos += K
            continue

        count = ht.get(kmer)

        if count < THRESHOLD:
            x.append(seq[pos].upper())
            pos += 1
        else:
            x.extend(seq[pos:pos + K].lower())
            total_masked += K
            pos += K

    leftover = len(seq) - len(x)
    if leftover:
        x.extend(seq[-leftover:])

    assert len(x) == len(seq), (record.name, len(x), len(seq))

    outfp.write(">%s\n%s\n" % (record.name, "".join(x)))

print 'total k-mers seen:', total
print 'n softmasked:', total_masked
