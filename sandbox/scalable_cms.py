#! /usr/bin/env python2
#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#

import khmer, math

INITIAL_SIZE = 512
GROWTH_RATE = 2
N_TABLES = 4
ERROR=.01
ERROR_TIGHTENING_RATIO=0.5

class ScalableCounter(object):
    """
    A simple trial implementation of a dynamically scaling CountMin Sketch,
    based on a scalable Bloom filter.
    """
    tables = []
    table_counts = []
    table_capacity = []

    def __init__(self, k):
        self.k = k
        self.add_table()

    def add_table(self):
        n_tables = len(self.tables)
        size = INITIAL_SIZE * (GROWTH_RATE**n_tables)
        kh = khmer.new_counting_hash(self.k, size, N_TABLES)
        self.tables.append(kh)
        self.table_counts.append(0)

        error_tightening = ERROR_TIGHTENING_RATIO**n_tables

        capacity = int(N_TABLES * size * math.log(2) * math.log(2) / \
                       abs(math.log(ERROR * error_tightening)))
        
        self.table_capacity.append(capacity)

        print 'added new table of size %d/capacity %d (now %d tables total)' %\
              (size, capacity, n_tables + 1)

    def count(self, dna):
        n = len(self.tables) - 1

        # iterate over the tables
        while n >= 0:
            kh = self.tables[n]
            if kh.get(dna):
                kh.count(dna)
#                if self.table_counts[n] >= self.table_capacity[n]:
#                    kh.decrement_random()
                return
            
            n -= 1

        # do we need to expand the table capacity?
        if self.table_counts[-1] >= self.table_capacity[-1]:
            kh = self.tables[-1]
            old_fpr = khmer.calc_expected_collisions(kh)
            print 'last table is full! %d counts, FP rate %.3f' % (\
                self.table_counts[-1], old_fpr)
            self.add_table()
            
        kh = self.tables[-1]
        kh.count(dna)
        self.table_counts[-1] += 1

    def get(self, dna):
        n = len(self.tables) - 1
        while n >= 0:
            kh = self.tables[n]
            count = kh.get(dna)
            if count:
                return count

            n -= 1

        return 0
            
if __name__ == '__main__':
    import random
    thenums = []
    
    K = 20
    limit = 4**K

    print 'Creating new ScalableCounter: growth rate %d, error ratio %.2f, bound %.3f' % \
          (GROWTH_RATE, ERROR_TIGHTENING_RATIO, ERROR)
    kh = ScalableCounter(K)

    NN = 10000

    for i in range(NN):
        if i % 1000 == 0:
            print '...added', i
        r = random.randint(0, limit)
        dna = khmer.reverse_hash(r, K)
        kh.count(dna)
        thenums.append(dna)

    #####

    total = 0
    #print 'shuffling items'
    #random.shuffle(thenums)
    for i, dna in enumerate(thenums):
        count = kh.get(dna)
        total += count - 1

    for i in range(10):
        print thenums[i], kh.get(thenums[i])

    total_size = sum([ sum(t.hashsizes()) for t in kh.tables ])
    print 'total memory used: %.1fk' % (total_size / 1000.)
    print 'average miscount:', float(total) / float(NN)

