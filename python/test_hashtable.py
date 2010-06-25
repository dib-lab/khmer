import khmer

def test_no_collision():
    kh = khmer.new_hashtable(4, 86)

    kh.count('AAAA')
    assert kh.get('AAAA') == 1

    kh.count('TTTT')                    # reverse complement
    assert kh.get('TTTT') == 2

def test_collision():
    kh = khmer.new_hashtable(4, 85)

    kh.count('AAAA')
    assert kh.get('AAAA') == 1

    kh.count('TTTT')
    assert kh.get('TTTT') == 2

def test_complete_no_collision():
    kh = khmer.new_hashtable(4, 4**4)
    kt = khmer.new_ktable(4)

    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        kh.count(s)

    n_palindromes = 0
    n_rc_filled = 0
    n_fwd_filled = 0
    
    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        if kh.get(s):                   # string hashing is rc aware
            n_rc_filled += 1
        if kh.get(s) == 1:              # palindromes are singular
            n_palindromes += 1
        if kh.get(i):                   # int hashing is not rc aware
            n_fwd_filled += 1

    assert n_rc_filled == kt.n_entries(),  n_rc_filled
    assert n_palindromes == 16, n_palindromes # @CTB check this
    assert n_fwd_filled == kt.n_entries() / 2 + n_palindromes / 2, \
           n_fwd_filled

def test_complete_2_collision():
    kh = khmer.new_hashtable(4, 4**4 / 2)
    kt = khmer.new_ktable(4)

    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        kh.count(s)

    n_rc_filled = 0
    n_fwd_filled = 0
    
    for i in range(0, 128):
        s = kt.reverse_hash(i)
        if kh.get(s):                   # string hashing is rc aware
            n_rc_filled += 1
        if kh.get(i):                   # int hashing is not rc aware
            n_fwd_filled += 1

    assert n_rc_filled == 128,  n_rc_filled
    # @CTB assert n_fwd_filled == 100 # kt.n_entries() / 2, n_fwd_filled

def test_complete_4_collision():
    kh = khmer.new_hashtable(4, 4**4 / 4)
    kt = khmer.new_ktable(4)

    for i in range(0, kt.n_entries()):
        s = kt.reverse_hash(i)
        kh.count(s)

    n_rc_filled = 0
    n_fwd_filled = 0
    
    for i in range(0, 64):
        s = kt.reverse_hash(i)
        if kh.get(s):                   # string hashing is rc aware
            n_rc_filled += 1
        if kh.get(i):                   # int hashing is not rc aware
            n_fwd_filled += 1

    assert n_rc_filled == 64,  n_rc_filled
    # @CTB assert n_fwd_filled == kt.n_entries() / 2, n_fwd_filled

def test_maxcount():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.new_hashtable(4, 4**4)
    
    last_count = None
    for i in range(0, 10000):
        kh.count('AAAA')
        c = kh.get('AAAA')
        
        print last_count, c
        if c == last_count:
            break
        last_count = c

    assert c != 10000, "should not be able to count to 10000"
    assert c == 255                     # this will depend on HashcountType...

def test_maxcount_consume():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.new_hashtable(4, 4**4)

    s = "A"*10000
    kh.consume(s)

    c = kh.get('AAAA')
    assert c == 255, c          # this will depend on HashcountType...

def test_get_mincount():
    kh = khmer.new_hashtable(4, 4**4)

    s = "AAAAACGT"
    kh.consume(s)

    x = kh.get_min_count(s)
    assert x == 1
    
    kh.consume(s)
    x = kh.get_min_count(s)
    assert x == 2

def test_get_maxcount():
    kh = khmer.new_hashtable(4, 4**4)

    s = "AAAAACGT"
    kh.consume(s)

    x = kh.get_max_count(s)
    assert x == 2
    
    kh.consume(s)
    x = kh.get_max_count(s)
    assert x == 4

def test_get_maxcount_rc():
   kh = khmer.new_hashtable(4, 4**4)

   s = "AAAAACGT"
   src = "ACGTTTTT"
   kh.consume(s)

   x = kh.get_max_count(s)
   assert x == 2

   kh.consume(src)
   x = kh.get_max_count(s)
   assert x == 4

def test_get_mincount_rc():
   kh = khmer.new_hashtable(4, 4**4)

   s = "AAAAACGT"
   src = "ACGTTTTT"

   kh.consume(s)
   x = kh.get_min_count(s)
   assert x == 1

   kh.consume(src)
   x = kh.get_min_count(s)
   assert x == 2

def test_consume_fasta():
   kh = khmer.new_hashtable(25, 4**15)
   kh.consume_fasta("1000.fa")

   minCount = kh.get_min_count("GCGCCTGGCTCAAGAATGCGGCACTAACCGGCACTGCCGTGGTAAACAATCCGTTTTGGTGGAGCGTGGATGAAAAATTCTTCACTAATGCGCTGGCCACAAAACTCGCCGTCG")
   
   assert minCount > 0

DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

class Test_HashtableIntersect:
    def setup(self):
        self.hi = khmer.HashtableIntersect(10, *khmer.PRIMES_1m)

    def test_basic(self):
        hi = self.hi
        hi.consume(DNA)

        x = hi.get_min_count(DNA)
        assert x == 1, x

        x = hi.get_max_count(DNA)
        assert x == 1, x

        hi.consume('ATTCTGACTG')

        x = hi.get_min_count(DNA)
        assert x == 1, x

        x = hi.get_max_count(DNA)
        assert x == 2, x

    def test_collision_1(self):
        kt = khmer.new_ktable(10)
        
        GG = 'G' * 10                   # forward_hash: 1048575
        assert kt.forward_hash(GG) == 1048575

        collision_1 = 'AACGGTCGGA'      # forward_hash: 48572
        assert kt.forward_hash(collision_1) == 48572

        collision_2 = 'AACTTGTTAC'      # forward_hash: 38738
        assert kt.forward_hash(collision_2) == 38738

        # note, hash(GG) % 1000003 == hash(collision_1)
        # note, hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_1)

        assert hi._kh1.get(GG) == 2
        assert hi._kh2.get(GG) == 1

        assert hi.get_min_count(GG) == 1
        assert hi.get_max_count(GG) == 1

    def test_collision_2(self):
        kt = khmer.new_ktable(10)
        
        GG = 'G' * 10                   # forward_hash: 1048575
        assert kt.forward_hash(GG) == 1048575

        collision_1 = 'AACGGTCGGA'      # forward_hash: 48572
        assert kt.forward_hash(collision_1) == 48572

        collision_2 = 'AACTTGTTAC'      # forward_hash: 38738
        assert kt.forward_hash(collision_2) == 38738

        # note, hash(GG) % 1000003 == hash(collision_1)
        # note, hash(GG) % 1009837 == hash(collision_2)

        hi = self.hi
        hi.consume(GG)
        hi.consume(collision_2)

        assert hi._kh1.get(GG) == 1
        assert hi._kh2.get(GG) == 2

        assert hi.get_min_count(GG) == 1
        assert hi.get_max_count(GG) == 1
