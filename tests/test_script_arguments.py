#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
"""
Tests for various argument-handling code.
"""

import khmer_tst_utils as utils

import khmer.file


def test_check_space():
    fakelump_fa = utils.get_test_data('fakelump.fa')

    try:
        khmer.file.check_space([fakelump_fa], _testhook_free_space=0)
        assert 0, "this should fail"
    except SystemExit, e:
        print str(e)
        pass


def test_check_tablespace():
    try:
        khmer.file.check_space_for_hashtable(1e9, _testhook_free_space=0)
        assert 0, "this should fail"
    except SystemExit, e:
        print str(e)
        pass
