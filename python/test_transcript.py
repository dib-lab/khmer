import khmer
import os

K = 17
HASHTABLE_SIZE=4**12+1
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

def test_get_transcript():
   est = "AACCGGTTAAACCCGGGTTTAAAACCCCGGGGTTTT"

   filename = os.path.join(thisdir, 'test-transcript.fa')
   estfile = os.path.join(thisdir, 'test-est.fa')

   ht = khmer.new_hashbits(K, HASHTABLE_SIZE)

   ht.consume_fasta(estfile)
   total_reads, n_consumed = ht.consume_fasta(filename)

   print total_reads, n_consumed

   rm = ht.filter_file_connected(est, filename, total_reads)

   assert rm.get(0)
   assert rm.get(1)
   assert not rm.get(2)
   

