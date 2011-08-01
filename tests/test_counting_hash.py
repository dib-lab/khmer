import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer
MAX_COUNT=255

####

# from http://www.rsok.com/~jrm/printprimes.html
PRIMES_1m = [1000003, 1009837]
PRIMES_100m = [100009979, 100000007]
PRIMES_1b = [1000000007, 1000000919]
PRIMES_2b = [1999999973, 1999999943]
PRIMES_4b = [4000000007, 4000000009]
PRIMES_8b = [8000000011, 8000000051]

DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

class Test_CountingHash(object):
    def setup(self):
        self.hi = khmer._new_counting_hash(12, PRIMES_1m)

    def test_collision_1(self):
        kt = khmer.new_ktable(12)
        
        GG = 'G' * 12                   # forward_hash: 11184810
        assert khmer.forward_hash(GG, 12) == 11184810

        collision_1 = 'AAACGTATGACT'
        assert khmer.forward_hash(collision_1, 12) == 184777L

        collision_2 = 'AAATACCGAGCG'
        assert khmer.forward_hash(collision_2, 12) == 76603L

        # note, hash(GG) % 1000003 == hash(collision_1)
        # note, hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_1)

        assert hi.get(GG) == 1

    def test_collision_2(self):
        kt = khmer.new_ktable(12)
        
        GG = 'G' * 12                   # forward_hash: 11184810
        assert khmer.forward_hash(GG, 12) == 11184810

        collision_1 = 'AAACGTATGACT'
        assert khmer.forward_hash(collision_1, 12) == 184777L

        collision_2 = 'AAATACCGAGCG'
        assert khmer.forward_hash(collision_2, 12) == 76603L

        # hash(GG) % 1000003 == hash(collision_1)
        # hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_2)

        assert hi.get(GG) == 1

    def test_collision_3(self):
        kt = khmer.new_ktable(12)
        
        GG = 'G' * 12                   # forward_hash: 11184810
        assert khmer.forward_hash(GG, 12) == 11184810

        collision_1 = 'AAACGTATGACT'
        assert khmer.forward_hash(collision_1, 12) == 184777L

        collision_2 = 'AAATACCGAGCG'
        assert khmer.forward_hash(collision_2, 12) == 76603L

        # hash(GG) % 1000003 == hash(collision_1)
        # hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_1)
        hi.consume(collision_2)

        assert hi.get(GG) == 2

def test_3_tables():
    x = list(PRIMES_1m)
    x.append(1000005)
    
    hi = khmer._new_counting_hash(12, x)
      
    GG = 'G' * 12                   # forward_hash: 11184810
    assert khmer.forward_hash(GG, 12) == 11184810

    collision_1 = 'AAACGTATGACT'
    assert khmer.forward_hash(collision_1, 12) == 184777L

    collision_2 = 'AAATACCGAGCG'
    assert khmer.forward_hash(collision_2, 12) == 76603L

    collision_3 = 'AAACGTATCGAG'
    assert khmer.forward_hash(collision_3, 12) == 184755L

    # hash(GG) % 1000003 == hash(collision_1)
    # hash(GG) % 1009837 == hash(collision_2)
    # hash(GG) % 1000005 == hash(collision_3)
    hi.consume(GG)
    assert hi.get(GG) == 1
    
    hi.consume(collision_1)
    assert hi.get(GG) == 1
    
    hi.consume(collision_2)
    assert hi.get(GG) == 1
    
    hi.consume(collision_3)
    assert hi.get(GG) == 2

def test_simple_median():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume("AAAAAA")
    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print median, average, stddev
    assert median == 1
    assert average == 1.0
    assert stddev == 0.0

    hi.consume("AAAAAA")
    (median, average, stddev) = hi.get_median_count("AAAAAA")
    print median, average, stddev
    assert median == 2
    assert average == 2.0
    assert stddev == 0.0

    hi.consume("AAAAAT")
    (median, average, stddev) = hi.get_median_count("AAAAAAT")
    print median, average, stddev
    assert median == 2
    assert average == 1.5
    assert int(stddev*100) == 50        # .5
    
    hi.consume("AAAAAT")
    (median, average, stddev) = hi.get_median_count("AAAAAAT")
    print median, average, stddev
    assert median == 2
    assert average == 2.0
    assert stddev == 0.0

    hi.consume("AAAAAT")
    (median, average, stddev) = hi.get_median_count("AAAAAAT")
    print median, average, stddev
    assert median == 3
    assert average == 2.5
    assert int(stddev*100) == 50        # .5

def test_save_load():
    thisdir = os.path.dirname(__file__)
    inpath = os.path.join(thisdir, 'test-data/random-20-a.fa')    
    savepath = os.path.join(thisdir, 'tempcountingsave.ht')
    
    x = list(PRIMES_1m)
    x.append(1000005)
    
    hi = khmer._new_counting_hash(12, x)
    hi.consume_fasta(inpath)
    hi.save(savepath)

    ht = khmer._new_counting_hash(12, x)
    ht.load(savepath)

    x = hi.abundance_distribution(inpath)
    y = ht.abundance_distribution(inpath)

    assert sum(x) == 3966
    assert x == y, (x,y)

def test_trim_full():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA)

    seq, pos = hi.trim_on_abundance(DNA, 2)
    assert DNA == seq, seq
    
def test_trim_short():
    hi = khmer.new_counting_hash(6, 1e6, 2)

    hi.consume(DNA)
    hi.consume(DNA[:50])

    seq, pos = hi.trim_on_abundance(DNA, 2)
    assert DNA[:50] == seq, (seq, pos)
    assert hi.get(seq[-6:]) == 2
    assert hi.get(DNA[:51][-6:]) == 1
    

def test_maxcount():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.new_counting_hash(4, 4**4, 4)
    
    last_count = None
    for i in range(0, 10000):
        kh.count('AAAA')
        c = kh.get('AAAA')
        
        print last_count, c
        if c == last_count:
            break
        last_count = c

    assert c != 10000, "should not be able to count to 10000: %d" % c
    assert c == MAX_COUNT, c       # this will depend on HashcountType...

def test_maxcount_with_bigcount():
    # hashtable should not saturate, if use_bigcount is set.
    kh = khmer.new_counting_hash(4, 4**4, 4)
    kh.set_use_bigcount(True)
    
    last_count = None
    for i in range(0, 10000):
        kh.count('AAAA')
        c = kh.get('AAAA')
        
        print last_count, c
        if c == last_count:
            break
        last_count = c

    assert c == 10000, "should be able to count to 10000: %d" % c
    assert c != MAX_COUNT, c
