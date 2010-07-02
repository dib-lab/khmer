import os
import tempfile
import shutil

import khmer
import screed
from screed.fasta import fasta_iter

thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

def load_fa_seq_names(filename):
    fp = open(filename)
    records = list(fasta_iter(fp))
    names = [ r['name'] for r in records ]
    return names

class Test_Filter(object):
    def setup(self):
        self.tempdir = tempfile.mkdtemp()

    def teardown(self):
        shutil.rmtree(self.tempdir)

    def test_filter(self):
        ht = khmer.new_hashtable(10, 4**10)

        filename = os.path.join(thisdir, 'test_data/simple_1.fa')
        outname = os.path.join(self.tempdir, 'test_filter.out')

        (total_reads, n_consumed) = ht.consume_fasta(filename)
        assert total_reads == 3, total_reads
        assert n_consumed == 63, n_consumed

        (total_reads, n_seq_kept) = khmer.filter_fasta_file(ht, filename,
                                                            total_reads,
                                                            outname, 2)
        assert n_seq_kept == 2, n_seq_kept

        names = load_fa_seq_names(outname)
        assert names == ['1', '2']

    def test_filter_n(self):
        ht = khmer.new_hashtable(10, 4**10)

        filename = os.path.join(thisdir, 'test_data/simple_2.fa')
        outname = os.path.join(self.tempdir, 'test_filter.out')

        (total_reads, n_consumed) = ht.consume_fasta(filename)
        assert total_reads == 4, total_reads
        assert n_consumed == 63, n_consumed

        (total_reads, n_seq_kept) = khmer.filter_fasta_file(ht, filename,
                                                            total_reads,
                                                            outname, 1)
        assert n_seq_kept == 3, n_seq_kept

        names = load_fa_seq_names(outname)
        assert names == ['1', '2', '3']

    def test_readmask_1(self):
        filename = os.path.join(thisdir, 'test_data/simple_2.fa')
        outname = os.path.join(self.tempdir, 'test_filter.out')

        readmask = khmer.new_readmask(4)
        readmask.set(1, False)
        readmask.set(2, False)
        readmask.set(3, False)

        readmask.filter_fasta_file(filename, outname)

        names = load_fa_seq_names(outname)
        assert names == ['1'], names

    def test_readmask_2(self):
        filename = os.path.join(thisdir, 'test_data/simple_2.fa')
        outname = os.path.join(self.tempdir, 'test_filter.out')

        readmask = khmer.new_readmask(4)
        readmask.set(0, False)
        readmask.set(1, True)
        readmask.set(2, False)
        readmask.set(3, False)

        readmask.filter_fasta_file(filename, outname)

        names = load_fa_seq_names(outname)
        assert names == ['2'], names
