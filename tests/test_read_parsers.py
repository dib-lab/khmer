# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,invalid-name

# Tests for the ReadParser and Read classes.
from khmer import Read
from khmer import ReadParser
from screed import Record
from . import khmer_tst_utils as utils
import pytest
from functools import reduce  # pylint: disable=redefined-builtin


def test_read_type_basic():
    # test that basic properties of khmer.Read behave like screed.Record
    # Constructing without mandatory arguments should raise an exception
    with pytest.raises(TypeError):
        Read()

    name = "895:1:1:1246:14654 1:N:0:NNNNN"
    sequence = "ACGT"
    r = Read(name, sequence)
    s = Record(name, sequence)

    for x in (r, s):
        assert x.name == name
        assert x.sequence == sequence
        assert not hasattr(x, 'quality'), x
        assert not hasattr(x, 'description'), x


def test_read_quality_none():
    r = Read(name="test", sequence="ACGT", quality=None)
    assert not hasattr(r, 'quality')


def test_read_type_attributes():
    r = Read(sequence='ACGT', quality='good', name='1234', description='desc')
    assert r.sequence == 'ACGT'
    assert r.cleaned_seq == 'ACGT'
    assert r.quality == 'good'
    assert r.name == '1234'
    assert r.description == 'desc'


def test_read_type_cleaned_seq():
    r = Read(sequence='acgtnN', name='1234')
    assert r.sequence == 'acgtnN'
    assert r.cleaned_seq == 'ACGTAA'


def test_read_properties():

    # Note: Using a data file with only one read.
    rparser = ReadParser(utils.get_test_data("single-read.fq"))

    # Check the properties of all one reads in data set.
    for read in rparser:
        assert read.name == "895:1:1:1246:14654 1:N:0:NNNNN"
        assert read.sequence == "CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT"
        # if an attribute is empty it shouldn't exist
        assert not hasattr(read, 'annotations')
        assert read.quality == """][aaX__aa[`ZUZ[NONNFNNNNNO_____^RQ_"""


def test_read_properties_fa():

    # Note: Using a data file with only one read.
    rparser = ReadParser(utils.get_test_data("single-read.fa"))

    # Check the properties of all one reads in data set.
    for read in rparser:
        print(read.name)
        assert read.name == "895:1:1:1246:14654 1:N:0:NNNNN"
        assert read.sequence == "CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT"
        # if an attribute is empty it shouldn't exist
        assert not hasattr(read, 'quality')


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


@pytest.mark.multithread
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
        for _ in rparser:
            n_reads += 1
    except ValueError as err:
        assert "Sequence is empty" in str(err), str(err)
    assert rparser.num_reads == 1, "%d valid reads in file, got %d" % (
        n_reads, rparser.num_reads)


def test_gzip_decompression():
    reads_count = 0
    rparser = ReadParser(utils.get_test_data("100-reads.fq.gz"))
    for _ in rparser:
        reads_count += 1

    assert 100 == reads_count


def test_gzip_decompression_truncated():

    rparser = ReadParser(utils.get_test_data("100-reads.fq.truncated.gz"))
    try:
        for _ in rparser:
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))


def test_gzip_decompression_truncated_pairiter():

    rparser = ReadParser(utils.get_test_data("100-reads.fq.truncated.gz"))
    try:
        for _ in rparser.iter_read_pairs():
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))
    except ValueError as err:
        print(str(err))


def test_bzip2_decompression():

    reads_count = 0
    rparser = ReadParser(utils.get_test_data("100-reads.fq.bz2"))
    for _ in rparser:
        reads_count += 1

    assert 100 == reads_count


def test_bzip2_decompression_truncated():

    rparser = ReadParser(utils.get_test_data("100-reads.fq.truncated.bz2"))
    try:
        for _ in rparser:
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))
    except ValueError as err:
        print(str(err))


def test_bzip2_decompression_truncated_pairiter():

    rparser = ReadParser(utils.get_test_data("100-reads.fq.truncated.bz2"))
    try:
        for _ in rparser.iter_read_pairs():
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))
    except ValueError as err:
        print(str(err))


def test_badbzip2():
    try:
        rparser = ReadParser(utils.get_test_data("test-empty.fa.bz2"))
        for _ in rparser:
            pass
        assert 0, "this should fail"
    except OSError as err:
        print(str(err))
    except ValueError as err:
        print(str(err))


@pytest.mark.multithread
def test_with_multiple_threads(testfile="test-reads.fq.bz2"):

    import operator
    import threading

    reads_count_1thr = 0
    rparser = ReadParser(utils.get_test_data(testfile))
    for _ in rparser:
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


@pytest.mark.multithread
def test_with_multiple_threads_big():
    test_with_multiple_threads(testfile="test-large.fa")


@pytest.mark.multithread
def test_old_illumina_pair_mating():

    import threading

    rparser = ReadParser(utils.get_test_data("test-reads.fa"))

    def thread_1_runtime(rparser):
        for _ in rparser:
            pass

    def thread_2_runtime(rparser):
        for readnum, _ in enumerate(rparser):
            if 0 == readnum:
                pass

    t1 = threading.Thread(target=thread_1_runtime, args=[rparser])
    t2 = threading.Thread(target=thread_2_runtime, args=[rparser])

    t1.start()
    t2.start()

    t1.join()
    t2.join()


@pytest.mark.multithread
def test_casava_1_8_pair_mating():

    import threading

    # Note: This file, when used in conjunction with a 64 KiB per-thread
    #       prefetch buffer, tests the paired read mating logic with the
    #       Casava >= 1.8 read name format.
    rparser = ReadParser(utils.get_test_data("test-reads.fq.bz2"))

    def thread_1_runtime(rparser):
        for _ in rparser:
            pass

    def thread_2_runtime(rparser):
        for readnum, _ in enumerate(rparser):
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
        for _ in rparser:
            pass
        assert 0, "No exception raised on a truncated file"
    except ValueError as err:
        assert "Sequence is empty" in str(err), str(err)


def test_iterator_identities():

    rparser = \
        ReadParser(utils.get_test_data("test-abund-read-paired.fa"))
    assert rparser is rparser.__iter__()
    assert rparser is rparser.iter_reads()


@pytest.mark.known_failing
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
    matches = [(rp1, rp2) for rp1, rp2 in zip(read_pairs_1, read_pairs_2)
               if rp1[0].name == rp2[0].name]
    assert all(matches)  # Assert ALL the matches. :-]


@pytest.mark.linux
def test_read_pair_iterator_in_error_mode_xfail():

    rparser = \
        ReadParser(utils.get_test_data("test-abund-read-impaired.fa"))

    failed = True
    try:
        for _ in rparser.iter_read_pairs():
            pass
        failed = False
    except ValueError as exc:
        assert "Invalid read pair" in str(exc), str(exc)
    assert failed


def test_read_pair_iterator_in_error_mode_xfail_osxsafe():

    rparser = \
        ReadParser(utils.get_test_data("test-abund-read-impaired.fa"))

    failed = True
    try:
        for _ in rparser.iter_read_pairs():
            pass
        failed = False
    except ValueError:
        pass
    assert failed


@pytest.mark.known_failing
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
        ReadParser(utils.get_test_data("single-read.fq"), "a")
        assert 0, ("ReadParser's constructor shouldn't accept a character for "
                   "the number of threads")
    except TypeError as err:
        print(str(err))
    try:
        ReadParser("non-existent-file-name")
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


def test_clean_seq():
    for read in ReadParser(utils.get_test_data("test-abund-read-3.fa")):
        clean = read.sequence.upper().replace("N", "A")
        assert clean == read.cleaned_seq


def test_error_badly_formatted_file():
    fname = utils.get_temp_filename('badly-formatted.fa')
    with open(fname, 'w') as f:
        f.write("not-sequence")

    with pytest.raises(OSError) as e:
        ReadParser(fname)

    assert e.match("contains badly formatted sequence")


def test_error_file_does_not_exist():
    fname = utils.get_temp_filename('does-not-exist.fa')

    with pytest.raises(OSError) as e:
        ReadParser(fname)

    assert e.match("does not exist")

# vim: set filetype=python tabstop=4 softtabstop=4 shiftwidth=4 expandtab:
# vim: set textwidth=79:
