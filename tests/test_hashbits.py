import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer

import khmer
try:
   import screed
   from screed.fasta import fasta_iter
except ImportError:
   import nose
   raise nose.SkipTest

def test_filter_if_present():
   ht = khmer.new_hashbits(32, 1e6, 2)

   maskfile = os.path.join(thisdir, 'test-data', 'filter-test-A.fa')
   inputfile = os.path.join(thisdir, 'test-data', 'filter-test-B.fa')
   outfile = os.path.join(thisdir, 'test-data', 'filter-test-C.fa')

   ht.consume_fasta(maskfile)
   ht.filter_if_present(inputfile, outfile)

   records = list(fasta_iter(open(outfile)))
   assert len(records) == 1
   assert records[0]['name'] == '3'

def test_combine_pe():
   inpfile = os.path.join(thisdir, 'test-data', 'combine_parts_1.fa')
   ht = khmer.new_hashbits(32, 1, 1)

   ht.consume_partitioned_fasta(inpfile)
   assert ht.count_partitions() == (2, 0)

   s1 = "CATGCAGAAGTTCCGCAACCATACCGTTCAGT"
   pid1 = ht.get_partition_id(s1)
   
   s2 = "CAAATGTACATGCACTTAAAATCATCCAGCCG"
   pid2 = ht.get_partition_id(s2)

   assert pid1 == 2
   assert pid2 == 80293

   ht.join_partitions(pid1, pid2)
   
   pid1 = ht.get_partition_id(s1)
   pid2 = ht.get_partition_id(s2)

   assert pid1 == pid2
   assert ht.count_partitions() == (1, 0)

def test_load_partitioned():
   inpfile = os.path.join(thisdir, 'test-data', 'combine_parts_1.fa')
   ht = khmer.new_hashbits(32, 1, 1)

   ht.consume_partitioned_fasta(inpfile)
   assert ht.count_partitions() == (2, 0)

   s1 = "CATGCAGAAGTTCCGCAACCATACCGTTCAGT"
   assert ht.get(s1)
   
   s2 = "CAAATGTACATGCACTTAAAATCATCCAGCCG"
   assert ht.get(s2)

   s3 = "CATGCAGAAGTTCCGCAACCATACCGTTCAGTTCCTGGTGGCTA"[-32:]
   assert ht.get(s3)
