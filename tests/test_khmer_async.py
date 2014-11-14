#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
# pylint: disable=missing-docstring,protected-access
import khmer
import time

from screed.fasta import fasta_iter
import screed

import khmer_tst_utils as utils
from nose.plugins.attrib import attr

def teardown():
    utils.cleanup()

def test_n_processed():
    filename = utils.get_test_data('test-reads.fa')
    ht = khmer.new_counting_hash(20, 1e7, 4)
    
    asd = khmer.AsyncDiginorm(ht)
    asd.start(filename, 1000, 1)
    time.sleep(1)
    asd.stop()

    assert asd.n_processed() == 25000L

def test_n_kept():
    filename = utils.get_test_data('test-reads.fa')
    ht = khmer.new_counting_hash(20, 1e7, 4)
    
    asd = khmer.AsyncDiginorm(ht)
    asd.start(filename, 1000, 1)
    time.sleep(1)
    asd.stop()

    assert asd.n_kept() == 25000L

def test_n_parsed():
    filename = utils.get_test_data('test-reads.fa')
    ht = khmer.new_counting_hash(20, 1e7, 4)
    
    asd = khmer.AsyncDiginorm(ht)
    asd.start(filename, 1000, 1)
    time.sleep(1)
    asd.stop()

    assert asd.n_parsed() == 25000L

def test_async_diginorm():
    filename = utils.get_test_data('test-reads.fa')
    ht = khmer.new_counting_hash(20, 1e7, 4)
    
    asd = khmer.AsyncDiginorm(ht)
    asd.start(filename, 5, 1)
    time.sleep(1)
    asd.stop()

    assert asd.n_kept() < 24995L
