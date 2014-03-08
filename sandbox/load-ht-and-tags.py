#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. 
# Contact: khmer-project@idyll.org
#
import khmer
import sys
import os.path

K = 32
# HASHTABLE_SIZE=int(1e3)
HASHTABLE_SIZE = int(64e9)

ht = khmer.new_hashbits(K, HASHTABLE_SIZE, 2)


def main(basename, filenames):
    print 'K', K
    print 'HASHTABLE SIZE %g' % HASHTABLE_SIZE
    print '--'

    # populate the hash table and tag set
    for filename in filenames:
        print 'reading sequences and loading tagset from %s...' % (filename,)
        ht.consume_fasta_and_tag(filename)

    # save to a file (optional)
    print 'saving...'
    ht.save(basename + '.ht')
    print 'saving tagset...'
    ht.save_tagset(basename + '.tagset')

    # calculate the hashtable occupancy
    print '---'
    print 'hashtable occupancy:', ht.n_occupied() / float(HASHTABLE_SIZE)
    print '---'

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
