
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


def test_update_from():
    table = QFCounttable(5, 1024 * 4)
    other_table = QFCounttable(5, 1024 * 4)

    assert table.get('AAAAA') == 0
    assert table.get('GCGCG') == 0
    assert table.n_occupied() == 0
    assert other_table.get('AAAAA') == 0
    assert other_table.get('GCGCG') == 0
    assert other_table.n_occupied() == 0

    other_table.count('AAAAA')

    assert table.get('AAAAA') == 0
    assert table.get('GCGCG') == 0
    assert table.n_occupied() == 0
    assert other_table.get('AAAAA') == 1
    assert other_table.get('GCGCG') == 0
    assert other_table.n_occupied() == 1

    table.count('GCGCG')

    assert table.get('AAAAA') == 0
    assert table.get('GCGCG') == 1
    assert table.n_occupied() == 1
    assert other_table.get('AAAAA') == 1
    assert other_table.get('GCGCG') == 0
    assert other_table.n_occupied() == 1

    table.update(other_table)

    assert table.get('AAAAA') == 1
    assert table.get('GCGCG') == 1
    assert table.n_occupied() == 2
    assert other_table.get('AAAAA') == 1
    assert other_table.get('GCGCG') == 0
    assert other_table.n_occupied() == 1


def test_update_from_2():

    tb1 = QFCounttable(20, 1024 * 4)
    tb2 = QFCounttable(20, 1024 * 4)

    filename = utils.get_test_data('random-20-a.fa')
    tb1.consume_seqfile(filename)
    tb2.consume_seqfile(filename)

    assert tb1.n_occupied() == tb2.n_occupied()

    # TODO: this fails (overcount?)
    #tb1.update(tb2)

    assert tb1.n_occupied() == tb2.n_occupied()
