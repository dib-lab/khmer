import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer

####

# from http://www.rsok.com/~jrm/printprimes.html
PRIMES_1m = [1000003, 1009837]
PRIMES_100m = [100009979, 100000007]
PRIMES_1b = [1000000007, 1000000919]
PRIMES_2b = [1999999973, 1999999943]
PRIMES_4b = [4000000007, 4000000009]
PRIMES_8b = [8000000011, 8000000051]

DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

class Test_HashtableIntersect(object):
    def setup(self):
        self.hi = khmer._new_counting_hash(12, PRIMES_1m)

    def test_basic(self):
        hi = self.hi
        hi.consume(DNA)

        x = hi.get_min_count(DNA)
        assert x == 1, x

        x = hi.get_max_count(DNA)
        assert x == 1, x

        hi.consume('ATTCTGACTGCA')

        x = hi.get_min_count(DNA)
        assert x == 1, x

        x = hi.get_max_count(DNA)
        assert x == 2, x

    def test_collision_1(self):
        kt = khmer.new_ktable(12)
        
        GG = 'G' * 12                   # forward_hash: 11184810
        assert kt.forward_hash(GG) == 11184810

        collision_1 = 'AAACGTATGACT'      # forward_hash: 184777L
        assert kt.forward_hash(collision_1) == 184777L

        collision_2 = 'AAATACCGAGCG'      # forward_hash: 76603L
        assert kt.forward_hash(collision_2) == 76603L

        # note, hash(GG) % 1000003 == hash(collision_1)
        # note, hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_1)

        assert hi.get_min_count(GG) == 1
        assert hi.get_max_count(GG) == 1

    def test_collision_2(self):
        kt = khmer.new_ktable(12)
        
        GG = 'G' * 12                   # forward_hash: 11184810
        assert kt.forward_hash(GG) == 11184810

        collision_1 = 'AAACGTATGACT'      # forward_hash: 184777L
        assert kt.forward_hash(collision_1) == 184777L

        collision_2 = 'AAATACCGAGCG'      # forward_hash: 76603L
        assert kt.forward_hash(collision_2) == 76603L

        # note, hash(GG) % 1000003 == hash(collision_1)
        # note, hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_2)

        assert hi.get_min_count(GG) == 1
        assert hi.get_max_count(GG) == 1

    def test_collision_3(self):
        kt = khmer.new_ktable(12)
        
        GG = 'G' * 12                   # forward_hash: 11184810
        assert kt.forward_hash(GG) == 11184810

        collision_1 = 'AAACGTATGACT'      # forward_hash: 184777L
        assert kt.forward_hash(collision_1) == 184777L

        collision_2 = 'AAATACCGAGCG'      # forward_hash: 76603L
        assert kt.forward_hash(collision_2) == 76603L

        # note, hash(GG) % 1000003 == hash(collision_1)
        # note, hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_1)
        hi.consume(collision_2)

        assert hi.get_min_count(GG) == 2
        assert hi.get_max_count(GG) == 2
