
import random

from khmer import QFCounttable
import khmer
from tests import khmer_tst_utils as utils
from khmer import ReadParser

import random
import pytest


MAX_COUNT = 255
MAX_BIGCOUNT = 65535

sketchSize = 1048576


DNA = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"


def teardown():
    utils.cleanup()


@pytest.fixture(params=[khmer.QFCounttable,khmer.BufferedQFCounttable])
def getSketch(request):
    return request.param


def test_count_1(getSketch):
    print("start")
    hi = getSketch(12, sketchSize,8)

    kmer = 'G' * 12
    hashval = hi.hash('G' * 12)

    assert hi.get(kmer) == 0
    assert hi.get(hashval) == 0

    hi.count(kmer)
    assert hi.get(kmer) == 1
    assert hi.get(hashval) == 1

    hi.count(kmer)
    assert hi.get(kmer) == 2
    assert hi.get(hashval) == 2

    kmer = 'G' * 11

    with pytest.raises(ValueError):
        hi.hash(kmer)


def test_count_2(getSketch):
    hi = getSketch(12, sketchSize,8)
    print("done")
    kmer = 'G' * 12
    hashval = hi.hash('G' * 12)

    assert hi.get(kmer) == 0
    assert hi.get(hashval) == 0

    hi.count(kmer)
    assert hi.get(kmer) == 1
    assert hi.get(hashval) == 1

    hi.count(hashval)                     # count hashes same as strings
    assert hi.get(kmer) == 2
    assert hi.get(hashval) == 2


# def test_read_write(getSketch):
#     print("Start")
#     fname = str.encode(utils.get_temp_filename('zzz'))
#     rng = random.Random(1)
#     ctm = getSketch(20, sketchSize,8)
#
#     kmers = ["".join(rng.choice("ACGT") for _ in range(20))
#              for n in range(400)]
#     for kmer in kmers:
#         ctm.add(kmer)
#
#     ctm.save(fname)
#
#     # print("Finish")
#     # # on purpose choose parameters that are different from sct
#     ctm2 = getSketch.load(fname)
#     ctm2.load(fname)
#     # assert ctm.ksize() == ctm2.ksize()
#     # for kmer in kmers:
#     #     assert ctm.get(kmer) == ctm2.get(kmer)
#     #
#     #
#

def test_maxcount_with_bigcount(getSketch):
    # hashtable should not saturate, if use_bigcount is set.
    kh = getSketch(4, 128,8)

    last_count = None
    for _ in range(0, 10000):
        kh.count('AAAA')
        c = kh.get('AAAA')
        print(c)
        if c == last_count:
            break

        last_count = c

    assert c == 10000, "should be able to count to 1000: %d" % c


def test_get_ksize(getSketch):
    kh = getSketch(22, 16,8)
    assert kh.ksize() == 22


#
# def test_read_write():
#     rng = random.Random(1)
#
#     qf = QFCounttable(20, 1024 * 4)
#
#     kmers = ["".join(rng.choice("ACGT") for _ in range(20))
#              for n in range(400)]
#     for kmer in kmers:
#         qf.add(kmer)
#
#     fname = utils.get_temp_filename('zzz')
#
#     qf.save(fname)
#
#     # on purpose choose parameters that are different from sct
#     qf2 = QFCounttable.load(fname)
#     assert qf.ksize() == qf2.ksize()
#     for kmer in kmers:
#         assert qf.get(kmer) == qf2.get(kmer)
