
import random

from khmer import QFCounttable

from . import khmer_tst_utils as utils


def test_read_write():
    rng = random.Random(1)

    qf = QFCounttable(20, 1024 * 4)

    kmers = ["".join(rng.choice("ACGT") for _ in range(20))
             for n in range(400)]
    for kmer in kmers:
        qf.add(kmer)

    fname = utils.get_temp_filename('zzz')

    qf.save(fname)

    # on purpose choose parameters that are different from sct
    qf2 = QFCounttable.load(fname)
    assert qf.ksize() == qf2.ksize()
    for kmer in kmers:
        assert qf.get(kmer) == qf2.get(kmer)
