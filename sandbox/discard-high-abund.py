#! /usr/bin/env python
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
"""
Eliminate reads with median k-mer abundance higher than DESIRED_COVERAGE.

Do so stochastically, to simplify matters.

Example usage:

% python scripts/load-into-counting.py mda-test.ht data/mda-test.fa
% python scripts/discard-high-abund.py mda-test.ht data/mda-test.fa > new.fa

Parameters to adjust: DESIRED_COVERAGE
"""

import sys
import screed
import os
import khmer
import random

DESIRED_COVERAGE = 5


def main():
    ht_filename = sys.argv[1]
    contig_filename = sys.argv[2]

    print>>sys.stderr, 'loading ht from', ht_filename
    ht = khmer.new_counting_hash(1, 1, 1)
    ht.load(ht_filename)

    for record in screed.open(contig_filename):
        seq = record.sequence.upper()
        if 'N' in seq:
            seq = seq.replace('N', 'G')

        med, _, _ = ht.get_median_count(seq)

        if med > DESIRED_COVERAGE:
            ratio = float(med) / float(DESIRED_COVERAGE)
            ratio = ratio * 100.
            if random.randrange(int(ratio)) > 100:
                sys.stderr.write('SKIPPING %s\n' % record.name)
                continue                 # ELIMINATE sequence

        print '>%s\n%s' % (record.name, record.sequence)

if __name__ == '__main__':
    main()
