import os
import tempfile
import shutil

import khmer
try:
   import screed
   from screed.fasta import fasta_iter
except ImportError:
   pass

import nose

thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

def load_fa_seq_names(filename):
    try:
       fasta_iter
    except NameError:
       raise nose.SkipTest

    fp = open(filename)
    records = list(fasta_iter(fp))
    names = [ r['name'] for r in records ]
    return names

class Test_Filter(object):
    def setup(self):
        self.tempdir = tempfile.mkdtemp()

    def teardown(self):
        shutil.rmtree(self.tempdir)

    def test_abund(self):
        ht = khmer.new_hashtable(10, 4**10)

        filename = os.path.join(thisdir, 'test-abund-read.fa')
        outname = os.path.join(self.tempdir, 'test_abund.out')

        ht.consume_fasta(filename)
        ht.output_fasta_kmer_pos_freq(filename, outname)
        
        fd = open(outname, "r")

        output = fd.readlines()
        assert len(output) == 1

        output = output[0]
        output = output.strip().split()

        assert ['1']*(114-10+1) == output

        fd.close()

    def test_filter_limit_n(self):
        ht = khmer.new_hashtable(4, 4**4)

        filename = os.path.join(thisdir, 'test_data/simple_3.fa')
        outname = os.path.join(self.tempdir, 'test_filter.out')

        (total_reads, n_consumed) = ht.consume_fasta(filename)
        assert total_reads == 2, total_reads

        (total_reads, n_seq_kept) = khmer.filter_fasta_file_limit_n(ht, filename,
                                                                    total_reads,
                                                                    outname, 2,
                                                                    7)


        assert total_reads == 2
        assert n_seq_kept == 1 

 
        (total_reads, n_seq_kept) = khmer.filter_fasta_file_limit_n(ht, filename,
                                                                    total_reads,
                                                                    outname, 2,
                                                                    4)

        assert total_reads == 2
        assert n_seq_kept == 2


    def test_filter(self):
        ht = khmer.new_hashtable(10, 4**10)

        filename = os.path.join(thisdir, 'test_data/simple_1.fa')
        outname = os.path.join(self.tempdir, 'test_filter.out')

        (total_reads, n_consumed) = ht.consume_fasta(filename)
        assert total_reads == 3, total_reads
        assert n_consumed == 63, n_consumed

        (total_reads, n_seq_kept) = khmer.filter_fasta_file_any(ht, filename,
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

        (total_reads, n_seq_kept) = khmer.filter_fasta_file_any(ht, filename,
                                                                total_reads,
                                                                outname, 1)
        assert n_seq_kept == 3, n_seq_kept

        names = load_fa_seq_names(outname)
        assert names == ['1', '2', '3']

    def test_consume_build_readmask(self):
        ht = khmer.new_hashtable(10, 4**10)

        filename = os.path.join(thisdir, 'test_data/simple_2.fa')
        outname = os.path.join(self.tempdir, 'test_filter.out')

        # sequence #4 (index 3) is bad; the new readmask should have that.
        x = ht.consume_fasta_build_readmask(filename)
        (total_reads, n_consumed, readmask) = x
        
        assert total_reads == 4, total_reads
        assert n_consumed == 63, n_consumed
        assert readmask.get(0)
        assert readmask.get(1)
        assert readmask.get(2)
        assert not readmask.get(3)
        
    def test_consume_update_readmask(self):
        ht = khmer.new_hashtable(10, 4**10)

        filename = os.path.join(thisdir, 'test_data/simple_2.fa')
        outname = os.path.join(self.tempdir, 'test_filter.out')

        readmask = khmer.new_readmask(4)

        # sequence #4 (index 3) is bad; the new readmask should have that.
        (total_reads, n_consumed) = ht.consume_fasta(filename, 0, 0,
                                                     readmask, True)
        assert total_reads == 4, total_reads
        assert n_consumed == 63, n_consumed
        assert readmask.get(0)
        assert readmask.get(1)
        assert readmask.get(2)
        assert not readmask.get(3)

    def test_consume_no_update_readmask(self):
        ht = khmer.new_hashtable(10, 4**10)

        filename = os.path.join(thisdir, 'test_data/simple_2.fa')
        outname = os.path.join(self.tempdir, 'test_filter.out')

        readmask = khmer.new_readmask(4)

        # sequence #4 (index 3) is bad; the new readmask should NOT have that.
        (total_reads, n_consumed) = ht.consume_fasta(filename, 0, 0,
                                                     readmask, False)
        assert total_reads == 4, total_reads
        assert n_consumed == 63, n_consumed
        assert readmask.get(0)
        assert readmask.get(1)
        assert readmask.get(2)
        assert readmask.get(3)          # NOT updated

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
