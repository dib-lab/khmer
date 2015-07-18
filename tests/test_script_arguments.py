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
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import sys


# For map(long, [list of ints]) cross-version hackery
if sys.version_info.major > 2:
    long = int


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
    outfile = utils.get_test_data('truncated.fq')
    save_stderr, sys.stderr = sys.stderr, io.StringIO()

    parser = khmer_args.build_counting_args()
    args = parser.parse_args(['-M', '1e9'])

    try:
        mem_args = args.max_memory_usage
        tablesize = khmer_args.calculate_tablesize(args.max_tablesize,
                                                   args.n_tables, 'countgraph',
                                                   max_memory_use=mem_args)
        khmer.kfile.check_space_for_hashtable(outfile, tablesize,
                                              False, _testhook_free_space=0)
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

    outfile = utils.get_test_data('truncated')

    parser = khmer_args.build_counting_args()
    args = parser.parse_args(['-M', '1e9'])

    try:
        mem_args = args.max_memory_usage
        tablesize = khmer_args.calculate_tablesize(args.max_tablesize,
                                                   args.n_tables, 'countgraph',
                                                   max_memory_use=mem_args)
        khmer.kfile.check_space_for_hashtable(outfile, tablesize,
                                              True, _testhook_free_space=0)
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
                                             'max_tablesize',
                                             'max_memory_usage'])


def test_create_countgraph_1():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    countgraph = khmer_args.create_countgraph(args)
    expected_hashsz = utils.longify([2499997, 2499989, 2499983, 2499967])
    assert countgraph.hashsizes() == expected_hashsz, countgraph.hashsizes()
    assert sum(countgraph.hashsizes()) < max_mem, sum(countgraph.hashsizes())


def test_create_countgraph_2():
    # tests overriding ksize by passing into create_nodegraph explicitly.

    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    countgraph = khmer_args.create_countgraph(args, ksize=15)
    assert countgraph.ksize() == 15


def test_create_countgraph_3():
    # tests too-big ksize

    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    old_stderr = sys.stderr
    sys.stderr = capture = StringIO()

    try:
        countgraph = khmer_args.create_countgraph(args, ksize=35)
        assert 0, "should not reach this"
    except SystemExit:
        err = capture.getvalue()
        assert 'khmer only supports k-mer sizes <= 32.' in err, err
    finally:
        sys.stderr = old_stderr


def test_create_countgraph_4_multiplier():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    countgraph = khmer_args.create_countgraph(args, multiplier=2.0)
    assert sum(countgraph.hashsizes()) < max_mem / 2.0, \
        sum(countgraph.hashsizes())


def test_create_nodegraph_1():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    nodegraph = khmer_args.create_nodegraph(args)
    expected_hashsz = utils.longify([19999999, 19999981, 19999963, 19999927])
    assert nodegraph.hashsizes() == expected_hashsz, nodegraph.hashsizes()

    assert sum(nodegraph.hashsizes()) / \
        8.0 < max_mem, sum(nodegraph.hashsizes())


def test_create_nodegraph_2():
    # tests overriding ksize by passing into create_nodegraph explicitly.

    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    nodegraph = khmer_args.create_nodegraph(args, ksize=15)
    assert nodegraph.ksize() == 15


def test_create_nodegraph_3():
    # tests too-big ksize

    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    old_stderr = sys.stderr
    sys.stderr = capture = StringIO()

    try:
        nodegraph = khmer_args.create_nodegraph(args, ksize=35)
        assert 0, "should not reach this"
    except SystemExit:
        err = capture.getvalue()
        assert 'khmer only supports k-mer sizes <= 32.' in err, err


def test_create_nodegraph_4_multiplier():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    nodegraph = khmer_args.create_nodegraph(args, multiplier=2.0)
    assert sum(nodegraph.hashsizes()) / 8.0 < max_mem / 2.0, \
        sum(nodegraph.hashsizes())


def test_report_on_config_bad_hashtype():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    try:
        khmer_args.report_on_config(args, 'foograph')
        assert 0, "the previous statement should raise an exception"
    except AssertionError:
        raise
    except Exception as err:
        assert "unknown graph type: foograph" in str(err), str(err)


def test_fail_calculate_foograph_size():
    # tests unknown graph type

    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem)

    try:
        nodegraph = khmer_args.calculate_tablesize(max_tablesize, n_tables,
                                                   'foograph',
                                                   max_memory_use=max_mem)
        assert 0, "previous statement should fail"
    except AssertionError:
        raise
    except Exception as err:
        assert "unknown graph type: foograph" in str(err), str(err)
