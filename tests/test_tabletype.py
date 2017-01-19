"""
Tests for common features of *all* tables, including 1-bit presence/absence
tables.
"""

import pytest


from khmer import _Countgraph, _Counttable, _SmallCountgraph, _SmallCounttable
from khmer import _Nodegraph, _Nodetable

PRIMES_1m = [1000003, 1009837]


@pytest.fixture(params=[_Countgraph, _Counttable, _SmallCountgraph,
                        _SmallCounttable, _Nodegraph, _Nodetable])
def tabletype(request):
    return request.param


def test_presence(tabletype):
    tt = tabletype(12, PRIMES_1m)

    kmer = 'G' * 12
    hashval = tt.hash('G' * 12)

    assert tt.get(kmer) == 0
    assert tt.get(hashval) == 0

    tt.add(kmer)
    assert tt.get(kmer) == 1
    assert tt.get(hashval) == 1


def test_bad_create(tabletype):
    try:
        tt = tabletype(5, [])
    except ValueError as err:
        assert 'tablesizes needs to be one or more numbers' in str(err)


def test_get_ksize(tabletype):
    kh = tabletype(22, PRIMES_1m)
    assert kh.ksize() == 22


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

    hi.consume("AAAAAA")
    hashes = hi.get_kmer_hashes("AAAAAA")
    print(hashes)
    assert len(hashes) == 1
    assert hi.get(hashes[0]) == 1

    hi.consume("AAAAAT")
    hashes = hi.get_kmer_hashes("AAAAAAT")
    print(hashes)
    assert len(hashes) == 2
    assert hi.get(hashes[0]) >= 1
    assert hi.get(hashes[1]) == 1


def test_get_kmers(tabletype):
    hi = tabletype(6, PRIMES_1m)

    kmers = hi.get_kmers("AAAAAA")
    assert kmers == ["AAAAAA"]

    kmers = hi.get_kmers("AAAAAAT")
    assert kmers == ["AAAAAA", "AAAAAT"]

    kmers = hi.get_kmers("AGCTTTTC")
    assert kmers == ['AGCTTT', 'GCTTTT', 'CTTTTC']
