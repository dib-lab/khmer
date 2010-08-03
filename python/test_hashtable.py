import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

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
    assert c == 127                     # this will depend on HashcountType...

def test_maxcount_consume():
    # hashtable should saturate at some point so as not to overflow counter
    kh = khmer.new_hashtable(4, 4**4)

    s = "A"*10000
    kh.consume(s)

    c = kh.get('AAAA')
    assert c == 127, c          # this will depend on HashcountType...

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

def test_64bitshift():
   kh = khmer.new_hashtable(25, 4**10)
   fullstr = "GTATGCCAGCTCCAACTGGGCCGGTACGAGCAGGCCATTGCCTCTTGCCGCGATGCGTCGGCG"
   substr =    "ATGCCAGCTCCAACTGGGCCGGTACGAGCAGGCCATTGCCTCTTGC"
   
   kh.consume(fullstr)
   assert 0 < kh.get_min_count(substr), kh.get_min_count(substr)
   
def test_64bitshift_2():
   kh = khmer.new_hashtable(25, 4**10)
   fullstr = "GTATGCCAGCTCCAACTGGGCCGGTACGAGCAGGCCATTGCCTCTTGCCGCGATGCGTCGGCG"
   
   kh.consume(fullstr)
   for i in range(len(fullstr) - 25 + 1):
       substr = fullstr[i:i+25]
       assert kh.get(substr) > 0

class Test_ConsumeString(object):
    def setup(self):
        self.kh = khmer.new_hashtable(4, 4**4)

    def test_n_occupied(self):
        assert self.kh.n_occupied() == 0
        n = self.kh.consume('AAAA')
        assert self.kh.n_occupied() == 1
        n = self.kh.consume('AACT')
        assert self.kh.n_occupied() == 2

    def test_abundance_by_pos(self):
        kh = self.kh

        for i in range(0, 300):
            kh.count('ATCG')

        for i in range(0, 10):
            kh.count('ATGG')

        short_filename = os.path.join(thisdir, 'test-short.fa')
        dist = kh.fasta_count_kmers_by_position(short_filename, 6, 10)
        assert dist[4] == 1
        assert sum(dist) == 1
        
        dist = kh.fasta_count_kmers_by_position(short_filename, 6, 127)
        assert dist[0] == 1, dist[0]
        assert dist[2] == 1
        assert sum(dist) == 2

    def test_abundance_dist(self):
        dist = self.kh.abundance_distribution()
        assert dist[0] == 4**4
        assert sum(dist[1:]) == 0
        
        n = self.kh.consume('AAAA')
        n = self.kh.consume('AACT')
        
        dist = self.kh.abundance_distribution()
        assert sum(dist[1:]) == 2
        assert dist[1] == 2

        n = self.kh.consume('AAAA')
        n = self.kh.consume('AACT')

        dist = self.kh.abundance_distribution()
        assert sum(dist[1:]) == 2
        assert dist[2] == 2, dist

    def test_n_occupied_args(self):
        assert self.kh.n_occupied() == 0
        n = self.kh.consume('AAAA')
        assert self.kh.n_occupied(0, 1) == 1
        assert self.kh.n_occupied(1, 4**4) == 0

        hash = khmer.forward_hash('AACT', 4)
        n = self.kh.consume('AACT')
        assert self.kh.n_occupied(0, hash + 1) == 2
        assert self.kh.n_occupied(hash + 1, 4**4) == 0
        assert self.kh.n_occupied(hash, hash + 1) == 1

    def test_simple(self):
        n = self.kh.consume('AAAA')
        assert n == 1
        assert self.kh.get(0) == 1
        
    def test_simple_2(self):
        n = self.kh.consume('AAAAA')
        assert n == 2
        assert self.kh.get(0) == 2
        
    def test_simple_rc(self):
        n = self.kh.consume('TTTTT')
        assert n == 2
        assert self.kh.get(0) == 2
        
    def test_bounded(self):
        n = self.kh.consume('AAAAA', 1, 4**4)
        assert n == 0, n
        assert self.kh.get(0) == 0

    def test_bounded_2(self):
        n = self.kh.consume('AAAAA', 0, 1)
        assert n == 2, n
        assert self.kh.get(0) == 2
        
    def test_bounded_rc(self):
        n = self.kh.consume('TTTTT', 1, 4**4)
        assert n == 0, n
        assert self.kh.get(0) == 0
        
    def test_bounded_2_rc(self):
        n = self.kh.consume('TTTTT', 0, 1)
        assert n == 2, n
        assert self.kh.get(0) == 2

    def test_min_count(self):
        self.kh.consume('AAAA')

        count = self.kh.get_min_count('AAAA')
        assert count == 1
        
    def test_min_count_in_bound(self):
        self.kh.consume('AAAA')

        count = self.kh.get_min_count('AAAAA', 0, 1)
        assert count == 1
        
    def test_min_count_out_bound(self):
        self.kh.consume('AAAA')

        count = self.kh.get_min_count('AAAAA', 1, 4**4)
        assert count == 127

    def test_max_count(self):
        self.kh.consume('AAAA')

        count = self.kh.get_max_count('AAAA')
        assert count == 1
        
    def test_max_count_in_bound(self):
        self.kh.consume('AAAA')

        count = self.kh.get_max_count('AAAAA', 0, 1)
        assert count == 1
        
    def test_max_count_out_bound(self):
        self.kh.consume('AAAA')

        count = self.kh.get_max_count('AAAAA', 1, 4**4)
        assert count == 0

class Test_ExactGraphFu(object):
    def setup(self):
        self.ht = khmer.new_hashtable(12, 4**12)

    def test_counts(self):
        ht = self.ht
        ht.consume_fasta(os.path.join(thisdir, 'test-graph.fa'))

        kmer = "TTAGGACTGCAC"
        x = ht.calc_connected_graph_size(kmer)
        assert x == 69, x
        
        kmer = "TGCGTTTCAATC"
        x = ht.calc_connected_graph_size(kmer)
        assert x == 68, x

        kmer = "ATACTGTAAATA"
        x = ht.calc_connected_graph_size(kmer)
        assert x == 36, x

    def test_trim(self):
        ht = self.ht
        
        filename = os.path.join(thisdir, 'test-graph.fa')
        outfile = os.path.join(thisdir, 'test-graph.fa.out')
        ht.consume_fasta(filename)
        ht.trim_graphs(filename, 40, outfile)

        ht = khmer.new_hashtable(12, 4**12)
        ht.consume_fasta(outfile)

        x = ht.calc_connected_graph_size("TTAGGACTGCAC")
        assert x == 69, x
        
        x = ht.calc_connected_graph_size("TGCGTTTCAATC")
        assert x == 68, x
        
        x = ht.calc_connected_graph_size("ATACTGTAAATA")
        assert x == 0, x

    def test_graphsize_distrib(self):
        ht = self.ht
        ht.consume_fasta(os.path.join(thisdir, 'test-graph.fa'))
        x = ht.graphsize_distribution(200)

        assert sum(x) == 3, x
        assert x[69] == 1
        assert x[68] == 1
        assert x[36] == 1

    def test_graph_links_next_a(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "A")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_next_c(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "C")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_next_g(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "G")

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_next_t(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume(word[1:] + "T")

        x = ht.calc_connected_graph_size(word)
        assert x == 2
        
    def test_graph_links_prev_a(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("A" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_c(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("C" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_g(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("G" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

    def test_graph_links_prev_t(self):
        ht = self.ht
        word = "TGCGTTTCAATC"
        ht.consume(word)
        ht.consume("T" + word[:-1])

        x = ht.calc_connected_graph_size(word)
        assert x == 2

class Test_InexactGraphFu(object):
    def setup(self):
        self.ht = khmer.new_hashtable(12, 4**8+1)

    def test_trim(self):
        ht = self.ht
        filename = os.path.join(thisdir, 'test-graph.fa')
        outfile = os.path.join(thisdir, 'test-graph.fa.out')
        
        ht.consume_fasta(filename)
        ht.trim_graphs(filename, 40, outfile)

        ht = khmer.new_hashtable(12, 4**12)
        ht.consume_fasta(outfile)

        x = ht.calc_connected_graph_size("TTAGGACTGCAC")
        assert x >= 69, x
        
        x = ht.calc_connected_graph_size("TGCGTTTCAATC")
        assert x >= 68, x               # @CTB why 69??
        
        x = ht.calc_connected_graph_size("ATACTGTAAATA")
        assert x == 0, x

####

DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

class Test_HashtableIntersect(object):
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
        return                          # @CTB
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
        return                          # @CTB
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
