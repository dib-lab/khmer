#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#
"""
Tests for various argument-handling code.
"""
from __future__ import print_function, unicode_literals
from __future__ import absolute_import

import sys
import io
import collections
from . import khmer_tst_utils as utils

import argparse
import khmer.kfile
from khmer import khmer_args


def test_check_space():
    fakelump_fa = utils.get_test_data('fakelump.fa')

    save_stderr, sys.stderr = sys.stderr, io.StringIO()
    try:
        khmer.kfile.check_space(
            [fakelump_fa], force=False, _testhook_free_space=0)
        assert 0, "this should fail"
    except SystemExit as e:
        print(str(e))
    finally:
        sys.stderr = save_stderr


def test_check_tablespace():
    save_stderr, sys.stderr = sys.stderr, io.StringIO()

    parser = khmer_args.build_counting_args()
    args = parser.parse_args(['-M', '1e9'])

    try:
        khmer.kfile.check_space_for_hashtable(args, 'countgraph', force=False,
                                              _testhook_free_space=0)
        assert 0, "this should fail"
    except SystemExit as e:
        print(str(e))
    finally:
        sys.stderr = save_stderr


def test_check_space_force():
    fakelump_fa = utils.get_test_data('fakelump.fa')

    save_stderr, sys.stderr = sys.stderr, io.StringIO()
    try:
        khmer.kfile.check_space(
            [fakelump_fa], force=True, _testhook_free_space=0)
        assert True, "this should pass"
    except SystemExit as e:
        print(str(e))
    finally:
        sys.stderr = save_stderr


def test_check_tablespace_force():
    save_stderr, sys.stderr = sys.stderr, io.StringIO()

    parser = khmer_args.build_counting_args()
    args = parser.parse_args(['-M', '1e9'])

    try:
        khmer.kfile.check_space_for_hashtable(args, 'countgraph', True,
                                              _testhook_free_space=0)
        assert True, "this should pass"
    except SystemExit as e:
        print(str(e))
    finally:
        sys.stderr = save_stderr


def test_invalid_file_warn():
    save_stderr, sys.stderr = sys.stderr, io.StringIO()
    try:
        khmer.kfile.check_valid_file_exists(["nonexistent", "nonexistent2"])
        assert sys.stderr.getvalue().count("\n") == 2,  \
            "Should produce two warning lines"
    except SystemExit as e:
        print(str(e))
    finally:
        sys.stderr = save_stderr


FakeArgparseObject = collections.namedtuple('FakeArgs',
                                            ['ksize', 'n_tables',
                                             'min_tablesize',
                                             'max_memory_usage'])


def test_create_countgraph_1():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    min_tablesize = khmer_args.DEFAULT_MIN_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, min_tablesize, max_mem)

    countgraph = khmer_args.create_countgraph(args)
    assert countgraph.hashsizes() == [2499997L, 2499989L, 2499983L, 2499967L]
    assert sum(countgraph.hashsizes()) < max_mem, sum(countgraph.hashsizes())
