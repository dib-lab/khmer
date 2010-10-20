import os
thisdir = os.path.dirname(__file__)
thisdir = os.path.abspath(thisdir)

import khmer

import khmer
try:
   import screed
   from screed.fasta import fasta_iter
except ImportError:
   pass

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
