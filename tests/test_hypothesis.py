#! /usr/bin/env python

"""
Hypothesis tests for khmer.

Hypothesis is a library implementing a mix of property-based testing,
unit testing and fuzzing. Tests are similar to unit tests,
but specify strategies for generating input.
Hypothesis takes these strategies and tries to find inputs to falsify the tests.
"""

from __future__ import division, unicode_literals

from collections import Counter
from tempfile import NamedTemporaryFile

from hypothesis import given, example, strategies as st
from nose.plugins.attrib import attr

import screed
import khmer
from khmer import reverse_hash, forward_hash, forward_hash_no_rc


# TODO: we are only testing with fixed k, table size and number of tables for now.
KSIZE = 13
N_TABLES = 4
TABLE_SIZE = 1e6

# strategy for creating kmers. Alphabet is derived from nucleotides.
st_kmer = st.text("ACGT", min_size=KSIZE, max_size=KSIZE)

st_sequence = st.text("ACGT", min_size=KSIZE, max_size=1000)

# strategy for creating valid FASTA records
st_record = st.fixed_dictionaries({
             'name': st.characters(min_codepoint=32, max_codepoint=126),
             'sequence': st_sequence}
            )

st_invalid_record = st.fixed_dictionaries({
             'name': st.characters(),
             'sequence': st.characters()}
            )

# Reverse complement utilities.
TRANSLATE = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
def revcomp(kmer):
    return "".join(TRANSLATE[c] for c in kmer[::-1])


def lex_rc(kmer):
    """ Returns the smaller string between kmer and its reverse complement,
    in lexicographical order """
    # We abuse one property from reverse_hash here,
    # since its implementation always return the smaller between kmer and rc.
    return reverse_hash(forward_hash(kmer, len(kmer)), len(kmer))


# FASTA utilities
def fasta_build(records):
    return "".join(">{name}\n{sequence}\n".format(**r) for r in records)


@attr('hypothesis')
@given(st_kmer)
def test_forward_hash(kmer):
    ksize = len(kmer)

    rh = reverse_hash(forward_hash(kmer, ksize), ksize)

    # We need to check for either kmer or reverse complement,
    # due to how the hash functions are implemented.
    assert rh == kmer or rh == revcomp(kmer)


@attr('hypothesis')
@given(st_kmer)
def test_forward_hash_no_rc(kmer):
    ksize = len(kmer)

    rh = reverse_hash(forward_hash_no_rc(kmer, ksize), ksize)

    # forward_hash_no_rc only generates the hash for the kmer,
    # not the reverse complement. In this case we just need to check if
    # the kmer is the same.
    assert rh == kmer


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1))
def test_countgraph_undercounting(kmers):
    """Testing countgraph undercounting.

    A collections.Counter serves as an oracle for Count-Min sketches,
    since both implement a frequency counter interface."""

    oracle = Counter()
    countgraph = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)

    for kmer in kmers:
        oracle.update([kmer])
        countgraph.count(kmer)

    for kmer in oracle:
        # Our implementation only counts to 255 by default,
        # so we need to check for at most 255 in the comparison.
        assert countgraph.get(kmer) >= min(oracle[kmer], 255)


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1))
def test_countgraph_undercounting_consume(kmers):
    """Testing countgraph undercounting, using consume instead of count.

        A collections.Counter serves as an oracle for Count-Min sketches,
        since both implement a frequency counter interface."""

    oracle = Counter()
    countgraph = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)

    for kmer in kmers:
        oracle.update([kmer])
        countgraph.consume(kmer)

    for kmer in oracle:
        # Our implementation only counts to 255 by default,
        # so we need to check for at most 255 in the comparison.
        assert countgraph.get(kmer) >= min(oracle[kmer], 255)


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1))
def test_countgraph_undercounting_bigcounts(kmers):
    """Testing countgraph undercounting, using the bigcount feature.

    A collections.Counter serves as an oracle for Count-Min sketches,
    since both implement a frequency counter interface."""

    oracle = Counter()
    countgraph = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)
    countgraph.set_use_bigcount(True)

    for kmer in kmers:
        oracle.update([kmer])
        countgraph.count(kmer)

    for kmer in oracle:
        # The bigcount implementation counts to 65535,
        # so we need to check for this at most in the comparison.
        assert countgraph.get(kmer) >= min(oracle[kmer], 65535)


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1))
def test_countgraph_undercounting_bigcounts_consume(kmers):
    """Testing countgraph undercounting, using the bigcount feature and using consume.

    A collections.Counter serves as an oracle for Count-Min sketches,
    since both implement a frequency counter interface."""

    oracle = Counter()
    countgraph = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)
    countgraph.set_use_bigcount(True)

    for kmer in kmers:
        oracle.update([kmer])
        countgraph.consume(kmer)

    for kmer in oracle:
        # The bigcount implementation counts to 65535,
        # so we need to check for this at most in the comparison.
        assert countgraph.get(kmer) >= min(oracle[kmer], 65535)


@attr('hypothesis')
@given(st.one_of(st.lists(st_record),
                 st.lists(st_invalid_record)))
def test_countgraph_consume_fasta(records):
    """
    """
    cg_fasta = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)
    cg_string = khmer.Countgraph(KSIZE, TABLE_SIZE, N_TABLES)

    fasta_error = None

    fasta = fasta_build(records)
    with NamedTemporaryFile() as temp:
        temp.write(fasta.encode('utf-8'))
        temp.flush()

        try:
            cg_fasta.consume_fasta(temp.name)
        except OSError:
            assert fasta == ''
        except Exception as e:
            fasta_error = e

    for record in records:
        try:
            cg_string.consume(record['sequence'])
        except ValueError:
            assert fasta_error is not None

    assert cg_fasta.n_unique_kmers() == cg_string.n_unique_kmers()
    assert cg_fasta.n_occupied() == cg_string.n_occupied()
    assert cg_fasta.hashsizes() == cg_string.hashsizes()


@attr('hypothesis')
@attr('known_failing')
@given(st.lists(st_kmer, min_size=10, unique_by=lex_rc))
def test_hll_cardinality(kmers):
    """Testing HyperLogLog cardinality estimation.

    len(set) serves as an oracle for HyperLogLog.
    """

    # We use the 'unique_by' argument to guarantee we have only one of
    # {kmer, revcomp(kmer)} in the set: We would overestimate otherwise.

    oracle = set()
    hll = khmer.HLLCounter(0.01, KSIZE)

    for kmer in kmers:
        oracle.update([kmer])
        hll.consume_string(kmer)

    if len(oracle) < 100:
        # This is tricky: if the cardinality is smaller than 100,
        # HyperLogLog might over or underestimate by 1 or 2,
        # which is greater than the default error rate (1%).
        # So, we just check to see if the absolute difference is small.
        assert abs(len(oracle) - len(hll)) <= 2
    else:
        # HyperLogLog error is calculated based on an uniform hash function,
        # but every hash function has some bias, even if small.
        # We set the error rate previously to 1%,
        # but we check for 2% here.
        error = round(abs(len(hll) - len(oracle)) / len(oracle), 2)
        assert error <= 0.02


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=1), st.lists(st_kmer, min_size=1))
def test_hll_merge_commutativity(kmers_1, kmers_2):
    """
    """
    hll1 = khmer.HLLCounter(0.01, KSIZE)
    hll2 = khmer.HLLCounter(0.01, KSIZE)

    for kmer in kmers_1:
        hll1.consume_string(kmer)

    for kmer in kmers_2:
        hll2.consume_string(kmer)

    hll2.merge(hll1)
    hll1.merge(hll2)

    assert hll1.counters == hll2.counters
    assert len(hll1) == len(hll2)
    assert hll1.error_rate == hll2.error_rate
    assert hll1.ksize == hll2.ksize


@attr('hypothesis')
@given(st.lists(st_record))
def test_hll_consume_fasta(records):
    """
    """
    hll_fasta = khmer.HLLCounter(0.01, KSIZE)
    hll_string = khmer.HLLCounter(0.01, KSIZE)

    fasta = fasta_build(records)
    with NamedTemporaryFile() as temp:
        temp.write(fasta.encode('utf-8'))
        temp.flush()

        try:
            hll_fasta.consume_fasta(temp.name)
        except OSError:
            assert fasta == ''

    for record in records:
        hll_string.consume_string(record['sequence'])

    assert len(hll_fasta) == len(hll_string)


@attr('hypothesis')
@given(st.sets(st_kmer, min_size=1))
def test_nodegraph_presence(kmers):
    """Testing nodegraph for presence checking.

    A set serves as an oracle for Bloom Filters,
    since both implement set membership operations."""

    oracle = set()
    nodegraph = khmer.Nodegraph(KSIZE, TABLE_SIZE, N_TABLES)

    for kmer in kmers:
        oracle.update([kmer])
        nodegraph.count(kmer)

    for kmer in oracle:
        # Nodegraph and Countgraph have similar API,
        # but Nodegraph.get returns only {0,1}
        assert nodegraph.get(kmer) == 1


@attr('hypothesis')
@given(st.sets(st_kmer, min_size=1))
def test_nodegraph_update(kmers):
    """

    Modeled after test_nodegraph:test_update_from_2
    """
    ng1 = khmer.Nodegraph(KSIZE, TABLE_SIZE, N_TABLES)
    ng2 = khmer.Nodegraph(KSIZE, TABLE_SIZE, N_TABLES)

    for kmer in kmers:
        ng1.count(kmer)
        ng2.count(kmer)

    assert ng1.n_occupied() == ng2.n_occupied()

    ng1.update(ng2)
    assert ng1.n_occupied() == ng2.n_occupied()


@attr('hypothesis')
@given(st.sets(st_kmer, min_size=1), st.sets(st_kmer, min_size=1))
def test_nodegraph_update_commutativity(kmers_1, kmers_2):
    """
    """
    ng1 = khmer.Nodegraph(KSIZE, TABLE_SIZE, N_TABLES)
    ng2 = khmer.Nodegraph(KSIZE, TABLE_SIZE, N_TABLES)

    for kmer in kmers_1:
        ng1.count(kmer)

    for kmer in kmers_2:
        ng2.count(kmer)

    ng2.update(ng1)
    ng1.update(ng2)
    assert ng1.n_occupied() == ng2.n_occupied()


@attr('hypothesis')
@given(st.lists(st_record))
def test_nodegraph_consume_fasta(records):
    """
    """
    ng_fasta = khmer.Nodegraph(KSIZE, TABLE_SIZE, N_TABLES)
    ng_string = khmer.Nodegraph(KSIZE, TABLE_SIZE, N_TABLES)

    fasta = fasta_build(records)
    with NamedTemporaryFile() as temp:
        temp.write(fasta.encode('utf-8'))
        temp.flush()

        try:
            ng_fasta.consume_fasta(temp.name)
        except OSError:
            assert fasta == ''

    for record in records:
        ng_string.consume(record['sequence'])

    assert ng_fasta.n_unique_kmers() == ng_string.n_unique_kmers()
    assert ng_fasta.n_occupied() == ng_string.n_occupied()
    assert ng_fasta.hashsizes() == ng_string.hashsizes()


@attr('hypothesis')
@given(st.lists(st_kmer, min_size=10, unique_by=lex_rc),
       st.integers(min_value=int(TABLE_SIZE)))
def test_n_unique(kmers, table_size):
    oracle = set()
    ng = khmer.Nodegraph(KSIZE, table_size, N_TABLES)
    cg = khmer.Countgraph(KSIZE, table_size, N_TABLES)

    for kmer in kmers:
        oracle.update([kmer])
        ng.consume(kmer)
        cg.consume(kmer)

    assert len(oracle) == ng.n_unique_kmers(), (len(oracle), ng.n_unique_kmers())
    assert len(oracle) == cg.n_unique_kmers(), (len(oracle), cg.n_unique_kmers())

    assert ng.n_occupied() == cg.n_occupied(), (ng.n_occupied(), cg.n_occupied())

    assert ng.hashsizes() == cg.hashsizes(), (ng.hashsizes(), cg.hashsizes())


@attr('hypothesis')
@given(st.lists(st_record), st.lists(st_record))
@example([{"name": "a", "sequence": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"}],
         [{"name": "1", "sequence": "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"},
          {"name": "2", "sequence": "GAGATCAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"},
          {"name": "3", "sequence": "AGAGATACACAAGATAGAGAGACCCAGGAGGGGG"}])
def test_nodegraph_filter_if_present(mask, records):
    """Testing nodegraph.filter_if_present

    filter_if_present removes all reads from an input file already covered by
    any kmer in the countgraph.
    """
    ng = khmer.Nodegraph(KSIZE, TABLE_SIZE, N_TABLES)

    mask_fasta = fasta_build(mask)
    with NamedTemporaryFile() as maskfile:
        maskfile.write(mask_fasta.encode('utf-8'))
        maskfile.flush()

        try:
            ng.consume_fasta(maskfile.name)
        except OSError:
            assert mask_fasta == ''

    input_fasta = fasta_build(records)
    with NamedTemporaryFile() as outfile:
        with NamedTemporaryFile() as inputfile:
            inputfile.write(mask_fasta.encode('utf-8'))
            inputfile.write(input_fasta.encode('utf-8'))
            inputfile.flush()

            try:
                ng.filter_if_present(inputfile.name, outfile.name)
            except OSError:
                assert input_fasta == ''

        with screed.open(outfile.name) as filtered:
            filtered_records = [r for r in filtered]

    assert all(m['sequence'] not in r['sequence']
               for m in mask
               for r in filtered_records)


#def test_nodegraph_combine_pe(kmers):
#  - count_partitions()
#  - get_partition_od
#  - join_partitions

#def test_nodegraph_count_within_radius
#  - count_kmers_within_radius

#def test_nodegraph_count_kmer_degree
#  - kmer_degree

#def test_nodegraph_save_load_tagset
#  - add_tag
#  - save_tagset
#  - load_tagset (clear_tags=False)

#def test_nodegraph_stop_traverse
#  - add_stop_tag
#  - consume_fasta_and_tag
#  - do_subset_partition
#  - merge_subset

#def test_nodegraph_find_unpart
#  - consume_fasta_and_tag
#  - do_subset_partition
#  - merge_subset
#  - count_partitions
#  - find_unpart

#def test_nodegraph_find_stoptags
#  - identify_stoptags_by_position

#def test_nodegraph_ksize
#  - ksize

#def test_nodegraph_get_hashsizes
#  - get_hashsizes

#def test_nodegraph_extract_unique_paths
#  - extract_unique_paths

#def test_nodegraph_get_raw_tables
#  - get_raw_tables

#def test_nodegraph_get_median_count
#  - get_median_count

#def test_countgraph_median_at_least
#  - median_at_least

#def test_countgraph_get_kmer_counts
#  - get_kmer_counts

#def test_countgraph_get_kmer_hashes
#  - get_kmer_hashes

#def test_countgraph_get_kmer
#  - get_kmers

#def test_countgraph_abundance_distribution
#  - abundance_distribution

#def test_countgraph_trim_on_abundance
#  - trim_on_abundance

#def test_countgraph_find_spectral_error_positions
#  - find_spectral_error_positions

#def test_countgraph_get_min_count
#  - get_min_count

#def test_countgraph_get_max_count
#  - get_max_count

#def test_countgraph_get_median_count
#  - get_median_count

#def test_countgraph_get_tags_and_positions
#  - consume_and_tag
#  - get_tags_and_positions
#  - find_all_tags_list

#def test_countgraph_consume_and_retrieve_tags_empty
#  - consume
#  - get_tags_and_positions
#  - find_all_tags_list
