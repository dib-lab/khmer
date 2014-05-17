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
ERROR=.1
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

        print 'new table:', n_tables, size, capacity

    def count(self, dna):
        n = len(self.tables) - 1
        while n >= 0:
            kh = self.tables[n]
            if kh.get(dna):
                kh.count(dna)
                if self.table_counts[n] >= self.table_capacity[n]:
                    kh.decrement_random()
                    pass
                return
            
            n -= 1

        if self.table_counts[-1] >= self.table_capacity[-1]:
            kh = self.tables[-1]
            fpr = khmer.calc_expected_collisions(kh)
            print 'adding new table', self.table_counts[-1], fpr
            self.add_table()
            
        kh = self.tables[-1]
        kh.count(dna)
        self.table_counts[-1] += 1

    def get(self, dna):
        n = len(self.tables) - 1
        while n >= 0:
            kh = self.tables[n]
            return kh.get(dna)

            n -= 1

        return 0
            
if __name__ == '__main__':
    import random
    thenums = []
    
    K = 20
    limit = 4**(K-1)

    kh = ScalableCounter(K)

    NN = 10000

    for i in range(NN):
        r = random.randint(0, limit)
        dna = khmer.reverse_hash(r, K-1)
        dna = 'A'+dna
        kh.count(dna)
        thenums.append(dna)

    #####

    total = 0
    #random.shuffle(thenums)
    for i, dna in enumerate(thenums):
        count = kh.get(dna)
        total += count - 1

    for i in range(10):
        print thenums[i], kh.get(thenums[i])

    print 'average miscount:', float(total) / float(NN)
