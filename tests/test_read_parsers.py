from __future__ import print_function
from __future__ import absolute_import
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#

# Tests for the ReadParser and Read classes.


import khmer
from khmer import ReadParser
from . import khmer_tst_utils as utils
from nose.plugins.attrib import attr
from functools import reduce


def test_read_properties():

    # Note: Using a data file with only one read.
    rparser = ReadParser(utils.get_test_data("single-read.fq"))

    # Check the properties of all one reads in data set.
    for read in rparser:
        assert read.name == "895:1:1:1246:14654 1:N:0:NNNNN"
        assert read.sequence == "CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT"
        assert read.annotations == ""
        assert read.quality == """][aaX__aa[`ZUZ[NONNFNNNNNO_____^RQ_"""


def test_with_default_arguments():

    read_names = []
    # Note: Using a data file where read names are just integers on [0,99).
    rparser = ReadParser(utils.get_test_data("random-20-a.fa"))

    for read in rparser:
        read_names.append(int(read.name))

    # "Derandomize".
    read_names.sort()

    # Each read number should match the corresponding name.
    for m, n in enumerate(read_names):
        assert m == n


def test_num_reads():
    """Test ReadParser.num_reads"""
    reads_count = 0
    rparser = ReadParser(utils.get_test_data("100-reads.fq.gz"))
    for _ in rparser:
        reads_count += 1

    assert reads_count == 100
    assert rparser.num_reads == 100


@attr('multithread')
def test_num_reads_threads():
    """Test threadsaftey of ReadParser's read counting"""
    import threading

    def count_reads(rparser):
        for _ in rparser:
            pass

    n_threads = 4
    threads = []
    rparser = ReadParser(utils.get_test_data("100-reads.fq.gz"))
    for _ in range(n_threads):
        thr = threading.Thread(target=count_reads, args=[rparser, ])
        threads.append(thr)
        thr.start()
    for thr in threads:
        thr.join()

    assert rparser.num_reads == 100


def test_num_reads_truncated():

    n_reads = 0
    rparser = ReadParser(utils.get_test_data("truncated.fq"))
    try:
        for read in rparser:
            n_reads += 1
    except ValueError as err:
        assert "Sequence is empty" in str(err), str(err)
    assert rparser.num_reads == 1, "%d valid reads in file, got %d" % (
        n_reads, rparser.num_reads)


def test_gzip_decompression():
    reads_count = 0
    rparser = ReadParser(utils.get_test_data("100-reads.fq.gz"))
    for read in rparser:
        reads_count += 1

    assert 100 == reads_count


def test_gzip_decompression_truncated():

    rparser = ReadParser(utils.get_test_data("100-reads.fq.truncated.gz"))
    try:
        for read in rparser:
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))


def test_gzip_decompression_truncated_pairiter():

    rparser = ReadParser(utils.get_test_data("100-reads.fq.truncated.gz"))
    try:
        for read in rparser.iter_read_pairs():
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))


def test_bzip2_decompression():

    reads_count = 0
    rparser = ReadParser(utils.get_test_data("100-reads.fq.bz2"))
    for read in rparser:
        reads_count += 1

    assert 100 == reads_count


def test_bzip2_decompression_truncated():

    rparser = ReadParser(utils.get_test_data("100-reads.fq.truncated.bz2"))
    try:
        for read in rparser:
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))


def test_bzip2_decompression_truncated_pairiter():

    rparser = ReadParser(utils.get_test_data("100-reads.fq.truncated.bz2"))
    try:
        for read in rparser.iter_read_pairs():
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))


def test_badbzip2():
    try:
        rparser = ReadParser(utils.get_test_data("test-empty.fa.bz2"))
        for read in rparser:
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))
    except ValueError as err:
        print(str(err))


@attr('multithread')
def test_with_multiple_threads(testfile="test-reads.fq.bz2"):

    import operator
    import threading

    reads_count_1thr = 0
    rparser = ReadParser(utils.get_test_data(testfile))
    for read in rparser:
        reads_count_1thr += 1

    def count_reads(rparser, counters, tnum):
        counters[tnum] = reduce(operator.add, (1 for read in rparser))

    N_THREADS = 4
    threads = []
    reads_counts_per_thread = [0] * N_THREADS
    rparser = ReadParser(utils.get_test_data(testfile))
    for tnum in range(N_THREADS):
        t = \
            threading.Thread(
                target=count_reads,
                args=[rparser, reads_counts_per_thread, tnum]
            )
        threads.append(t)
        t.start()
    for t in threads:
        t.join()

    assert reads_count_1thr == sum(reads_counts_per_thread), \
        reads_counts_per_thread


@attr('multithread')
def test_with_multiple_threads_big():
    test_with_multiple_threads(testfile="test-large.fa")


@attr('multithread')
def test_old_illumina_pair_mating():

    import threading

    rparser = ReadParser(utils.get_test_data("test-reads.fa"))

    def thread_1_runtime(rparser):
        for read in rparser:
            pass

    def thread_2_runtime(rparser):
        for readnum, read in enumerate(rparser):
            if 0 == readnum:
                pass

    t1 = threading.Thread(target=thread_1_runtime, args=[rparser])
    t2 = threading.Thread(target=thread_2_runtime, args=[rparser])

    t1.start()
    t2.start()

    t1.join()
    t2.join()


@attr('multithread')
def test_casava_1_8_pair_mating():

    import threading

    # Note: This file, when used in conjunction with a 64 KiB per-thread
    #       prefetch buffer, tests the paired read mating logic with the
    #       Casava >= 1.8 read name format.
    rparser = ReadParser(utils.get_test_data("test-reads.fq.bz2"))

    def thread_1_runtime(rparser):
        for read in rparser:
            pass

    def thread_2_runtime(rparser):
        for readnum, read in enumerate(rparser):
            if 0 == readnum:
                pass
            # assert "895:1:1:1761:13189 2:N:0:NNNNN" == read.name, read.name

    t1 = threading.Thread(target=thread_1_runtime, args=[rparser])
    t2 = threading.Thread(target=thread_2_runtime, args=[rparser])

    t1.start()
    t2.start()

    t1.join()
    t2.join()


def test_read_truncated():

    rparser = ReadParser(utils.get_test_data("truncated.fq"))
    try:
        for read in rparser:
            pass
        assert 0, "No exception raised on a truncated file"
    except ValueError as err:
        assert "Sequence is empty" in str(err), str(err)


def test_iterator_identities():

    rparser = \
        ReadParser(utils.get_test_data("test-abund-read-paired.fa"))
    assert rparser is rparser.__iter__()
    assert rparser is rparser.iter_reads()


@attr('known_failing')
def test_read_pair_iterator_in_error_mode():
    assert 0

    rparser = \
        ReadParser(utils.get_test_data("test-abund-read-paired.fa"))

    # If walks like an iterator and quacks like an iterator...
    rpi = rparser.iter_read_pairs()
    assert "__iter__" in dir(rpi)
    assert "next" in dir(rpi)

    # Are the alleged pairs actually pairs?
    read_pairs_1 = []
    for read_1, read_2 in rpi:
        read_pairs_1.append([read_1, read_2])
        assert read_1.name[: 19] == read_2.name[: 19]

    # Reload parser.
    # Note: No 'rewind' or 'reset' capability at the time of this writing.
    rparser = \
        ReadParser(utils.get_test_data("test-abund-read-paired.fa"))

    # Ensure that error mode is the default mode.
    read_pairs_2 = []
    for read_1, read_2 \
            in rparser.iter_read_pairs(ReadParser.PAIR_MODE_ERROR_ON_UNPAIRED):
        read_pairs_2.append([read_1, read_2])
    matches = \
        list(map(
            lambda rp1, rp2: rp1[0].name == rp2[0].name,
            read_pairs_1, read_pairs_2
        ))
    assert all(matches)  # Assert ALL the matches. :-]


def test_read_pair_iterator_in_error_mode_xfail():

    rparser = \
        ReadParser(utils.get_test_data("test-abund-read-impaired.fa"))

    failed = True
    try:
        for rpair in rparser.iter_read_pairs():
            pass
        failed = False
    except ValueError as exc:
        assert "Invalid read pair detected" in str(exc), str(exc)
    assert failed


@attr('known_failing')
def test_read_pair_iterator_in_ignore_mode():
    assert 0

    rparser = \
        ReadParser(utils.get_test_data("test-abund-read-impaired.fa"))

    read_pairs = []
    for read_1, read_2 \
            in rparser.iter_read_pairs(ReadParser.PAIR_MODE_IGNORE_UNPAIRED):
        read_pairs.append([read_1, read_2])
        assert read_1.name[: 19] == read_2.name[: 19]
    assert 2 == len(read_pairs)


def test_constructor():

    # Note: Using a data file with only one read.
    try:
        rparser = ReadParser(utils.get_test_data("single-read.fq"), "a")
        assert 0, ("ReadParser's constructor shouldn't accept a character for "
                   "the number of threads")
    except TypeError as err:
        print(str(err))
    try:
        rparser = ReadParser("non-existent-file-name")
        assert 0, "ReadParser shouldn't accept a non-existant file name"
    except ValueError as err:
        print(str(err))
    except OSError as err:
        print(str(err))


def test_iternext():
    try:
        rparser = ReadParser(utils.get_test_data("fakelump.fa.stoptags.txt"))
        read_pairs = []
        for read_1, read_2 in rparser.iter_read_pairs():
            read_pairs.append(read_1, read_2)
        assert 0, "Shouldn't be able to iterate over non FASTA file"
    except OSError as err:
        print(str(err))
    except ValueError as err:
        print(str(err))
# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
