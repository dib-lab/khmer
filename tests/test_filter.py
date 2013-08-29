import khmer
import screed
from screed.fasta import fasta_iter

import khmer_tst_utils as utils

def teardown():
   utils.cleanup()

def load_fa_seq_names(filename):
    fp = open(filename)
    records = list(fasta_iter(fp))
    names = [ r['name'] for r in records ]
    return names

class Test_Filter(object):
    def test_abund(self):
        ht = khmer.new_hashtable(10, 4**10)

        filename = utils.get_test_data('test-abund-read.fa')
        outname = utils.get_temp_filename('test_abund.out')

        ht.consume_fasta(filename)
        ht.output_fasta_kmer_pos_freq(filename, outname)
        
        fd = open(outname, "r")

        output = fd.readlines()
        assert len(output) == 1

        output = output[0]
        output = output.strip().split()

        assert ['1']*(114-10+1) == output

        fd.close()

def test_filter_sodd():
   K = 32
   HASHTABLE_SIZE=int(8e7)
   N_HT = 4
   MAX_SODD=3
   
   ht = khmer.new_hashbits(K, HASHTABLE_SIZE, N_HT)
   filename = utils.get_test_data('../../data/high-sodd.fa')

   ht.consume_fasta(filename)

   seq = "CGTTAGTTGCGGTGCCGACCGGCAAACTTGGTTTTGCCAAAAATTTTTACAGTTAGAAATTATTCACAAAGTTGCACCGGAATTCGGTTACAAACGTCATTCTAACTAAT"
   trim_seq, trim_at = ht.trim_on_sodd(seq, MAX_SODD)
   assert trim_seq == "CGTTAGTTGCGGTGCCGACCGGCAAACTTGGT"

   seq = "ACAAAATTCCACATATAGTCATAATTGTGGGCAATTTTCGTCCCAAATTAGTTAGAATGACGTTTGTAACCGAATTCCGGTGCAACTTTGTGAATAATTTCTAACTGTAAAAAT"
   trim_seq, trim_at = ht.trim_on_sodd(seq, MAX_SODD)
   assert trim_seq == "ACAAAATTCCACATATAGTCATAATTGTGGGCAATT"

   seq = "GCACGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTG"
   trim_seq, trim_at = ht.trim_on_sodd(seq, MAX_SODD)
   assert trim_seq == seq
