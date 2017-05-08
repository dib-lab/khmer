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

import sys
import pytest


from . import khmer_tst_utils as utils
import khmer
from khmer import _Countgraph, _Counttable, _SmallCountgraph, _SmallCounttable
from khmer import _Nodegraph, _Nodetable
from khmer import ReadParser
import screed


PRIMES_1m = [1000003, 1009837]


# all the table types!
@pytest.fixture(params=[_Countgraph, _Counttable, _SmallCountgraph,
                        _SmallCounttable, _Nodegraph, _Nodetable])
def tabletype(request):
    return request.param

# For map(long, [list of ints]) cross-version hackery
if sys.version_info.major > 2:
    long = int  # pylint: disable=redefined-builtin
    unicode = str


def test_presence(tabletype):
    # basic get/add test
    tt = tabletype(12, PRIMES_1m)

    kmer = 'G' * 12
    hashval = tt.hash('G' * 12)

    assert tt.get(kmer) == 0
    assert tt.get(hashval) == 0

    tt.add(kmer)
    assert tt.get(kmer) == 1
    assert tt.get(hashval) == 1


def test_bad_create(tabletype):
    # creation should fail w/bad parameters
    try:
        tt = tabletype(5, [])
    except ValueError as err:
        assert 'tablesizes needs to be one or more numbers' in str(err)


def test_get_ksize(tabletype):
    # ksize() function.
    kh = tabletype(22, PRIMES_1m)
    assert kh.ksize() == 22


def test_hash(tabletype):
    # hashing of strings -> numbers.
    kh = tabletype(5, PRIMES_1m)
    x = kh.hash("ATGGC")
    assert type(x) == long


def test_hash_bad_dna(tabletype):
    # hashing of bad dna -> succeeds w/o complaint
    kh = tabletype(5, PRIMES_1m)

    x = kh.hash("ATGYC")


def test_hash_bad_length(tabletype):
    # hashing of bad dna length -> error
    kh = tabletype(5, PRIMES_1m)

    with pytest.raises(ValueError):
        x = kh.hash("ATGGGC")

    with pytest.raises(ValueError):
        x = kh.hash("ATGG")


def test_reverse_hash(tabletype):
    # hashing of strings -> numbers.
    kh = tabletype(5, PRIMES_1m)

    try:
        x = kh.reverse_hash(15)
    except ValueError:
        pytest.skip("reverse_hash not implemented on this table type")

    assert isinstance(x, (unicode, str))


def test_hashsizes(tabletype):
    # hashsizes method.
    kh = tabletype(5, PRIMES_1m)
    assert kh.hashsizes() == PRIMES_1m


def test_add_hashval(tabletype):
    # test add(hashval)
    kh = tabletype(5, PRIMES_1m)
    x = kh.hash("ATGGC")
    y = kh.add(x)
    assert y

    z = kh.get(x)
    assert z == 1


def test_add_dna_kmer(tabletype):
    # test add(dna)
    kh = tabletype(5, PRIMES_1m)
    x = kh.add("ATGGC")
    assert x

    z = kh.get("ATGGC")
    assert z == 1


def test_add_bad_dna_kmer(tabletype):
    # even with 'bad' dna, should succeed.
    kh = tabletype(5, PRIMES_1m)

    x = kh.add("ATYGC")


def test_get_hashval(tabletype):
    # test get(hashval)
    kh = tabletype(5, PRIMES_1m)
    hashval = kh.hash("ATGGC")
    kh.add(hashval)

    z = kh.get(hashval)
    assert z == 1


def test_get_hashval_rc(tabletype):
    # test get(hashval)
    kh = tabletype(4, PRIMES_1m)
    hashval = kh.hash("ATGC")
    rc = kh.hash("GCAT")

    assert hashval == rc


def test_get_dna_kmer(tabletype):
    # test get(dna)
    kh = tabletype(5, PRIMES_1m)
    hashval = kh.hash("ATGGC")
    kh.add(hashval)

    z = kh.get("ATGGC")
    assert z == 1


def test_get_bad_dna_kmer(tabletype):
    # test get(dna) with bad dna; should be fine.
    kh = tabletype(5, PRIMES_1m)

    kh.hash("ATYGC")


def test_consume_and_count(tabletype):
    tt = tabletype(6, PRIMES_1m)

    x = "ATGCCGATGCA"
    num_kmers = tt.consume(x)
    assert num_kmers == len(x) - tt.ksize() + 1   # num k-mers consumed

    for start in range(len(x) - 6 + 1):
        assert tt.get(x[start:start + 6]) == 1


def test_consume_and_count_bad_dna(tabletype):
    # while we don't specifically handle bad DNA, we should at least be
    # consistent...
    tt = tabletype(6, PRIMES_1m)

    x = "ATGCCGNTGCA"
    num_kmers = tt.consume(x)

    for start in range(len(x) - 6 + 1):
        assert tt.get(x[start:start + 6]) == 1


def test_consume_short(tabletype):
    # raise error on too short when consume is run
    tt = tabletype(6, PRIMES_1m)

    x = "ATGCA"
    with pytest.raises(ValueError):
        tt.consume(x)


def test_get_kmer_counts(tabletype):
    hi = tabletype(6, PRIMES_1m)

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


def test_get_kmer_hashes(tabletype):
    hi = tabletype(6, PRIMES_1m)

    hashes = hi.get_kmer_hashes("ACGTGCGT")
    print(hashes)
    assert len(hashes) == 3
    assert hashes[0] == hi.hash("ACGTGC")
    assert hashes[1] == hi.hash("CGTGCG")
    assert hashes[2] == hi.hash("GTGCGT")


def test_get_min_count(tabletype):
    hi = tabletype(6, PRIMES_1m)

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


def test_get_kmers(tabletype):
    hi = tabletype(6, PRIMES_1m)

    kmers = hi.get_kmers("AAAAAA")
    assert kmers == ["AAAAAA"]

    kmers = hi.get_kmers("AAAAAAT")
    assert kmers == ["AAAAAA", "AAAAAT"]

    kmers = hi.get_kmers("AGCTTTTC")
    assert kmers == ['AGCTTT', 'GCTTTT', 'CTTTTC']


def test_trim_on_abundance(tabletype):
    hi = tabletype(6, PRIMES_1m)

    x = "ATGGCAGTAGCAGTGAGC"
    hi.consume(x[:10])

    (y, pos) = hi.trim_on_abundance(x, 1)
    assert pos == 10
    assert x[:pos] == y


def test_trim_below_abundance(tabletype):
    hi = tabletype(6, PRIMES_1m)

    x = "ATGGCAGTAGCAGTGAGC"
    x_rc = screed.rc(x)
    hi.consume(x_rc[:10])

    print(len(x))

    (y, pos) = hi.trim_below_abundance(x, 0)
    assert pos == len(x) - hi.ksize() + 1
    assert x[:pos] == y


DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"


def test_find_spectral_error_positions(tabletype):
    kh = tabletype(8, PRIMES_1m)

    kh.consume(DNA[:30])

    for n in range(len(DNA) - 8 + 1):
        print(n, kh.get(DNA[n:n + 8]))

    posns = kh.find_spectral_error_positions(DNA, 0)
    assert posns == [30], posns


def test_find_spectral_error_positions_6(tabletype):
    kh = tabletype(8, PRIMES_1m)

    kh.consume(DNA[1:])

    for n in range(len(DNA) - 8 + 1):
        print(n, kh.get(DNA[n:n + 8]))

    posns = kh.find_spectral_error_positions(DNA, 0)
    assert posns == [0], posns


def test_find_spectral_error_positions_5(tabletype):
    kh = tabletype(8, PRIMES_1m)

    kh.consume(DNA[:10])
    kh.consume(DNA[11:])

    posns = kh.find_spectral_error_positions(DNA, 0)
    assert posns == [10], posns


def test_consume_seqfile_reads_parser(tabletype):
    kh = tabletype(5, PRIMES_1m)
    rparser = ReadParser(utils.get_test_data('test-fastq-reads.fq'))

    kh.consume_seqfile_with_reads_parser(rparser)

    kh2 = tabletype(5, PRIMES_1m)
    for record in screed.open(utils.get_test_data('test-fastq-reads.fq')):
        kh2.consume(record.sequence)

    assert kh.get('CCGGC') == kh2.get('CCGGC')


def test_consume_seqfile(tabletype):
    kh = tabletype(5, PRIMES_1m)
    kh.consume_seqfile(utils.get_test_data('test-fastq-reads.fq'))

    kh2 = tabletype(5, PRIMES_1m)
    for record in screed.open(utils.get_test_data('test-fastq-reads.fq')):
        kh2.consume(record.sequence)

    assert kh.get('CCGGC') == kh2.get('CCGGC')


def test_save_load(tabletype):
    kh = tabletype(5, PRIMES_1m)
    savefile = utils.get_temp_filename('tablesave.out')

    # test add(dna)
    x = kh.add("ATGGC")
    z = kh.get("ATGGC")
    assert z == 1

    kh.save(savefile)

    # should we provide a single load function here? yes, probably. @CTB
    if tabletype == _Countgraph:
        loaded = khmer.load_countgraph(savefile)
    elif tabletype == _Counttable:
        loaded = khmer.load_counttable(savefile)
    elif tabletype == _SmallCountgraph:
        loaded = khmer.load_countgraph(savefile, small=True)
    elif tabletype == _SmallCounttable:
        loaded = khmer.load_counttable(savefile, small=True)
    elif tabletype == _Nodegraph:
        loaded = khmer.load_nodegraph(savefile)
    elif tabletype == _Nodetable:
        loaded = khmer.load_nodetable(savefile)
    else:
        raise Exception("unknown tabletype")

    z = loaded.get('ATGGC')
    assert z == 1


def test_get_bigcount(tabletype):
    # get_bigcount should return false by default
    tt = tabletype(12, PRIMES_1m)

    assert not tt.get_use_bigcount()


def test_set_bigcount(tabletype):
    supports_bigcount = [_Countgraph, _Counttable]
    tt = tabletype(12, PRIMES_1m)

    if tabletype in supports_bigcount:
        tt.set_use_bigcount(True)

        for i in range(300):
            tt.add('G'*12)
        assert tt.get('G'*12) == 300

    else:
        with pytest.raises(ValueError):
            tt.set_use_bigcount(True)


def test_abund_dist_A(tabletype):
    A_filename = utils.get_test_data('all-A.fa')

    kh = tabletype(4, PRIMES_1m)
    tracking = khmer._Nodegraph(4, PRIMES_1m)

    kh.consume_seqfile(A_filename)
    dist = kh.abundance_distribution(A_filename, tracking)

    print(dist[:10])
    assert sum(dist) == 1
    assert dist[0] == 0


def test_abund_dist_A_readparser(tabletype):
    A_filename = utils.get_test_data('all-A.fa')
    rparser = ReadParser(A_filename)

    kh = tabletype(4, PRIMES_1m)
    tracking = khmer._Nodetable(4, PRIMES_1m)

    kh.consume_seqfile(A_filename)
    dist = kh.abundance_distribution(A_filename, tracking)

    print(dist[:10])
    assert sum(dist) == 1
    assert dist[0] == 0
