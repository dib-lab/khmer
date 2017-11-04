"""
Tests for common features of *all* tables, including 1-bit presence/absence
tables.

This tests:

* method implementations in common between all "table types", including
  error handling code;
* basic API compatibility between all "table types", for methods with
  their own implementations.

Error handling & other class-specific code (e.g. in methods like
'load' and 'save' that are implemented differently for each type) will
be tested separately.  We can use code coverage to identify that
code...
"""

import math
import sys
import pytest

from . import khmer_tst_utils as utils
import khmer

from khmer import Countgraph, SmallCountgraph, Nodegraph
from khmer import Nodetable, Counttable, SmallCounttable, QFCounttable
from khmer import CyclicCounttable

from khmer import ReadParser

from .table_fixtures import (AnyTabletype, Tabletype, params_1m, PRIMES_1m,
                             QF_SIZE)
import screed


# For map(long, [list of ints]) cross-version hackery
if sys.version_info.major > 2:
    long = int  # pylint: disable=redefined-builtin
    unicode = str


def test_presence(AnyTabletype):
    # basic get/add test
    tt = AnyTabletype(12)

    kmer = 'G' * 12
    hashval = tt.hash('G' * 12)

    assert tt.get(kmer) == 0
    assert tt.get(hashval) == 0

    tt.add(kmer)
    assert tt.get(kmer) == 1
    assert tt.get(hashval) == 1

    tt.add(kmer)
    # Node* types can only tell presence/absence
    if 'Node' in tt.__class__.__name__:
        assert tt.get(kmer) == 1
        assert tt.get(hashval) == 1
    else:
        assert tt.get(kmer) == 2
        assert tt.get(hashval) == 2


def test_n_occupied(AnyTabletype):
    # basic get/add test
    tt = AnyTabletype(12)

    kmer = 'G' * 12

    assert tt.n_occupied() == 0
    assert tt.n_unique_kmers() == 0

    tt.add(kmer)
    assert tt.n_occupied() == 1
    assert tt.n_unique_kmers() == 1

    tt.add(kmer)
    # the CQF implementation we use can use more than one slot to represent
    # counts for a single kmer
    if not tt.__class__.__name__.startswith("QF"):
        assert tt.n_occupied() == 1
    else:
        assert tt.n_occupied() == 2
    assert tt.n_unique_kmers() == 1


def test_bad_create(Tabletype):
    # creation should fail w/bad parameters
    try:
        tt = Tabletype(5, [])
    except ValueError as err:
        assert 'tablesizes needs to be one or more numbers' in str(err)


def test_get_ksize(AnyTabletype):
    # ksize() function.
    kh = AnyTabletype(22)
    assert kh.ksize() == 22


def test_hash(AnyTabletype):
    # hashing of strings -> numbers.
    kh = AnyTabletype(5)
    x = kh.hash("ATGGC")
    assert type(x) == long


def test_hash_bad_dna(AnyTabletype):
    # hashing of bad dna -> succeeds w/o complaint
    kh = AnyTabletype(5)

    x = kh.hash("ATGYC")


def test_hash_bad_length(AnyTabletype):
    # hashing of bad dna length -> error
    kh = AnyTabletype(5)

    with pytest.raises(ValueError):
        x = kh.hash("ATGGGC")

    with pytest.raises(ValueError):
        x = kh.hash("ATGG")


def test_reverse_hash(AnyTabletype):
    # hashing of strings -> numbers.
    kh = AnyTabletype(5)

    try:
        x = kh.reverse_hash(15)
    except ValueError:
        pytest.skip("reverse_hash not implemented on this table type")

    assert isinstance(x, (unicode, str))


def test_hashsizes(AnyTabletype):
    # hashsizes method.
    kh = AnyTabletype(5)
    assert (kh.hashsizes() == PRIMES_1m or
            # CQF allocates some extra slots beyond what you request
            # exactly how many extra is an implementation detail
            kh.hashsizes()[0] >= QF_SIZE)


def test_add_hashval(AnyTabletype):
    # test add(hashval)
    kh = AnyTabletype(5)
    x = kh.hash("ATGGC")
    y = kh.add(x)
    assert y

    z = kh.get(x)
    assert z == 1


def test_add_dna_kmer(AnyTabletype):
    # test add(dna)
    kh = AnyTabletype(5)
    x = kh.add("ATGGC")
    assert x

    z = kh.get("ATGGC")
    assert z == 1


def test_add_bad_dna_kmer(AnyTabletype):
    # even with 'bad' dna, should succeed.
    kh = AnyTabletype(5)

    x = kh.add("ATYGC")


def test_get_hashval(AnyTabletype):
    # test get(hashval)
    kh = AnyTabletype(5)
    hashval = kh.hash("ATGGC")
    kh.add(hashval)

    z = kh.get(hashval)
    assert z == 1


def test_get_hashval_rc(AnyTabletype):
    # test get(hashval)
    kh = AnyTabletype(4)
    hashval = kh.hash("ATGC")
    rc = kh.hash("GCAT")

    assert hashval == rc


def test_get_dna_kmer(AnyTabletype):
    # test get(dna)
    kh = AnyTabletype(5)
    hashval = kh.hash("ATGGC")
    kh.add(hashval)

    z = kh.get("ATGGC")
    assert z == 1


def test_get_bad_dna_kmer(AnyTabletype):
    # test get(dna) with bad dna; should be fine.
    kh = AnyTabletype(5)

    kh.hash("ATYGC")


def test_consume_and_count(AnyTabletype):
    tt = AnyTabletype(6)

    x = "ATGCCGATGCA"
    num_kmers = tt.consume(x)
    assert num_kmers == len(x) - tt.ksize() + 1   # num k-mers consumed

    for start in range(len(x) - 6 + 1):
        assert tt.get(x[start:start + 6]) == 1


def test_consume_and_count_bad_dna(AnyTabletype):
    # while we don't specifically handle bad DNA, we should at least be
    # consistent...
    tt = AnyTabletype(6)

    x = "ATGCCGNTGCA"
    num_kmers = tt.consume(x)

    for start in range(len(x) - 6 + 1):
        assert tt.get(x[start:start + 6]) == 1


def test_consume_short(AnyTabletype):
    # raise error on too short when consume is run
    tt = AnyTabletype(6)

    x = "ATGCA"
    with pytest.raises(ValueError):
        tt.consume(x)


def test_get_kmer_counts(AnyTabletype):
    hi = AnyTabletype(6)

    hi.consume("AAAAAA")
    counts = hi.get_kmer_counts("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] == 1

    hi.consume("AAAAAA")
    counts = hi.get_kmer_counts("AAAAAA")
    print(counts)
    assert len(counts) == 1
    assert counts[0] >= 1

    hi.consume("AAAAAT")
    counts = hi.get_kmer_counts("AAAAAAT")
    print(counts)
    assert len(counts) == 2
    assert counts[0] >= 1
    assert counts[1] == 1


def test_get_kmer_hashes(AnyTabletype):
    hi = AnyTabletype(6)

    hashes = hi.get_kmer_hashes("ACGTGCGT")
    print(hashes)
    assert len(hashes) == 3
    assert hashes[0] == hi.hash("ACGTGC")
    assert hashes[1] == hi.hash("CGTGCG")
    assert hashes[2] == hi.hash("GTGCGT")


def test_get_min_count(AnyTabletype):
    hi = AnyTabletype(6)

    # master string, 3 k-mers
    x = "ACGTGCGT"

    hi.add("ACGTGC")  # 3
    hi.add("ACGTGC")
    hi.add("ACGTGC")

    hi.add("CGTGCG")  # 1

    hi.add("GTGCGT")  # 2
    hi.add("GTGCGT")

    counts = hi.get_kmer_counts(x)
    assert hi.get_min_count(x) == min(counts)
    assert hi.get_max_count(x) == max(counts)
    med, _, _ = hi.get_median_count(x)
    assert med == list(sorted(counts))[len(counts) // 2]


def test_get_kmers(AnyTabletype):
    hi = AnyTabletype(6)

    kmers = hi.get_kmers("AAAAAA")
    assert kmers == ["AAAAAA"]

    kmers = hi.get_kmers("AAAAAAT")
    assert kmers == ["AAAAAA", "AAAAAT"]

    kmers = hi.get_kmers("AGCTTTTC")
    assert kmers == ['AGCTTT', 'GCTTTT', 'CTTTTC']


def test_trim_on_abundance(AnyTabletype):
    hi = AnyTabletype(6)

    x = "ATGGCAGTAGCAGTGAGC"
    hi.consume(x[:10])

    (y, pos) = hi.trim_on_abundance(x, 1)
    assert pos == 10
    assert x[:pos] == y


def test_trim_below_abundance(AnyTabletype):
    hi = AnyTabletype(6)

    x = "ATGGCAGTAGCAGTGAGC"
    x_rc = screed.rc(x)
    hi.consume(x_rc[:10])

    print(len(x))

    (y, pos) = hi.trim_below_abundance(x, 0)
    assert pos == len(x) - hi.ksize() + 1
    assert x[:pos] == y


DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"


def test_find_spectral_error_positions(AnyTabletype):
    kh = AnyTabletype(8)

    kh.consume(DNA[:30])

    for n in range(len(DNA) - 8 + 1):
        print(n, kh.get(DNA[n:n + 8]))

    posns = kh.find_spectral_error_positions(DNA, 0)
    assert posns == [30], posns


def test_find_spectral_error_positions_6(AnyTabletype):
    kh = AnyTabletype(8)

    kh.consume(DNA[1:])

    for n in range(len(DNA) - 8 + 1):
        print(n, kh.get(DNA[n:n + 8]))

    posns = kh.find_spectral_error_positions(DNA, 0)
    assert posns == [0], posns


def test_find_spectral_error_positions_5(AnyTabletype):
    kh = AnyTabletype(8)

    kh.consume(DNA[:10])
    kh.consume(DNA[11:])

    posns = kh.find_spectral_error_positions(DNA, 0)
    assert posns == [10], posns


def test_consume_seqfile_reads_parser(AnyTabletype):
    kh = AnyTabletype(5)
    rparser = ReadParser(utils.get_test_data('test-fastq-reads.fq'))

    kh.consume_seqfile(rparser)

    kh2 = AnyTabletype(5)
    for record in screed.open(utils.get_test_data('test-fastq-reads.fq')):
        kh2.consume(record.sequence)

    assert kh.get('CCGGC') == kh2.get('CCGGC')


def test_consume_seqfile(AnyTabletype):
    kh = AnyTabletype(5)
    kh.consume_seqfile(utils.get_test_data('test-fastq-reads.fq'))

    kh2 = AnyTabletype(5)
    for record in screed.open(utils.get_test_data('test-fastq-reads.fq')):
        kh2.consume(record.sequence)

    assert kh.get('CCGGC') == kh2.get('CCGGC')


def test_save_load(Tabletype):
    kh = Tabletype(5)
    ttype = type(kh)
    savefile = utils.get_temp_filename('tablesave.out')

    # test add(dna)
    x = kh.add("ATGGC")
    z = kh.get("ATGGC")
    assert z == 1

    kh.save(savefile)

    # should we provide a single load function here? yes, probably. @CTB
    loaded = ttype.load(savefile)

    z = loaded.get('ATGGC')
    assert z == 1


def test_get_bigcount(Tabletype):
    # get_bigcount should return false by default
    tt = Tabletype(12)

    assert not tt.get_use_bigcount()


def test_set_bigcount(Tabletype):
    supports_bigcount = [Countgraph, Counttable, CyclicCounttable]
    tt = Tabletype(12)

    if type(tt) in supports_bigcount:
        tt.set_use_bigcount(True)

        for i in range(300):
            tt.add('G' * 12)
        assert tt.get('G' * 12) == 300

    else:
        with pytest.raises(ValueError):
            tt.set_use_bigcount(True)


def test_abund_dist_A(AnyTabletype):
    A_filename = utils.get_test_data('all-A.fa')

    kh = AnyTabletype(4)
    tracking = Nodegraph(4, 1, 1, primes=PRIMES_1m)

    kh.consume_seqfile(A_filename)
    dist = kh.abundance_distribution(A_filename, tracking)

    print(dist[:10])
    assert sum(dist) == 1
    assert dist[0] == 0


def test_abund_dist_A_readparser(AnyTabletype):
    A_filename = utils.get_test_data('all-A.fa')
    rparser = ReadParser(A_filename)

    kh = AnyTabletype(4)
    tracking = Nodegraph(4, 1, 1, primes=PRIMES_1m)

    kh.consume_seqfile(A_filename)
    dist = kh.abundance_distribution(rparser, tracking)

    print(dist[:10])
    assert sum(dist) == 1
    assert dist[0] == 0
