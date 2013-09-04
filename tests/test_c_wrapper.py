#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
#
import os
import khmer

import khmer_tst_utils as utils

reads_filename = utils.get_test_data('test-reads.fa')

N_READS = 25000

def teardown():
    utils.cleanup()

class GoodException(Exception):
    pass

def callback_raise(info, n_reads, other):
    raise GoodException

def setup():
    khmer.set_reporting_callback(None)

def teardown():
    khmer.reset_reporting_callback()

def test_raise_in_consume_fasta():
    return ## @CTB
    kh = khmer.new_hashtable(4, 4**4)

    try:
        n, _ = kh.consume_fasta(reads_filename, 0, 0, callback_raise)
        print n
        assert 0
    except GoodException:
        pass
    except:
        raise
