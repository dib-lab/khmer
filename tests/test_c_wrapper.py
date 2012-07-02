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
    
def test_raise_in_consume_fasta_build_readmask():
    return ## @CTB
    kh = khmer.new_hashtable(4, 4**4)

    try:
        kh.consume_fasta_build_readmask(reads_filename, 0, 0, callback_raise)
        assert 0
    except GoodException:
        pass
    except:
        raise

def test_raise_in_consume_fasta():
    return ## @CTB
    kh = khmer.new_hashtable(4, 4**4)

    try:
        n, _ = kh.consume_fasta(reads_filename, 0, 0, None, False, callback_raise)
        print n
        assert 0
    except GoodException:
        pass
    except:
        raise

def test_raise_in_readmask_filter_fasta_file():
    return # @@CTB fix
    readmask = khmer.new_readmask(N_READS)

    tstfile = utils.get_temp_filename('tst')

    try:
        readmask.filter_fasta_file(reads_filename, tstfile, callback_raise)
        assert 0
    except GoodException:
        pass
    except:
        raise

def test_raise_in_fasta_file_to_minmax():
    return # @@CTB fix
    ht = khmer.new_hashtable(4, 4**4)

    try:
        ht.fasta_file_to_minmax(reads_filename, N_READS, None, callback_raise)
        assert 0
    except GoodException:
        pass
    except:
        raise
    
def test_raise_in_filter_fasta_file_max():
    return ## @CTB
    ht = khmer.new_hashtable(4, 4**4)

    khmer.reset_reporting_callback()

    mmt = ht.fasta_file_to_minmax(reads_filename, N_READS)

    try:
        ht.filter_fasta_file_any(mmt, 2, None, callback_raise)
        assert 0
    except GoodException:
        pass
    except:
        raise

def test_bad_mmt_in_filter_fasta_file_max():
    ht = khmer.new_hashtable(4, 4**4)

    try:
        ht.filter_fasta_file_any("hi", 2)
        assert 0
    except TypeError:
        pass                            # expected
    except:
        raise

def test_bad_readmask_in_filter_fasta_file_limit_n():
    ht = khmer.new_hashtable(4, 4**4)

    khmer.reset_reporting_callback()

    mmt = ht.fasta_file_to_minmax(reads_filename, N_READS)

    try:
        ht.filter_fasta_file_limit_n(mmt, 2, 2, "hi")
        assert 0
    except TypeError:
        pass
    except:
        raise
        
def test_bad_readmask_in_filter_fasta_file_max():
    ht = khmer.new_hashtable(4, 4**4)

    khmer.reset_reporting_callback()

    mmt = ht.fasta_file_to_minmax(reads_filename, N_READS)

    try:
        ht.filter_fasta_file_any(mmt, 2, "hi")
        assert 0
    except TypeError:
        pass                            # expected
    except:
        raise

def test_nonbool_in_consume_fasta():
    return ## @CTB

    kh = khmer.new_hashtable(4, 4**4)

    try:
        kh.consume_fasta(reads_filename, 0, 0, "hi", False, callback_raise)
        assert 0
    except TypeError:
        pass
    except:
        raise

