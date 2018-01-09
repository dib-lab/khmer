# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2014-2015, Michigan State University.
# Copyright (C) 2015-2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
# pylint: disable=missing-docstring,invalid-name,no-member
"""
Tests for various argument-handling code.
"""

import sys
import io
import collections
from . import khmer_tst_utils as utils
import pytest

import khmer.kfile
from khmer import khmer_args
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


# For map(long, [list of ints]) cross-version hackery
if sys.version_info.major > 2:
    long = int  # pylint: disable=redefined-builtin


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


@pytest.mark.parametrize('graph_type,buckets_per_byte', [
    ('countgraph', 1),
    ('smallcountgraph', 2),
    ('nodegraph', 8),
])
def test_check_tablespace(graph_type, buckets_per_byte):
    oldstderr = sys.stderr
    sys.stderr = StringIO()

    outfile = utils.get_test_data('truncated.fq')
    parser = khmer_args.build_counting_args()
    args = parser.parse_args(['-M', '16G'])

    buckets_per_table = khmer_args.calculate_graphsize(args, graph_type)
    total_buckets = buckets_per_table * args.n_tables
    space_needed = total_buckets / buckets_per_byte

    # First, try with insufficient space
    with pytest.raises(SystemExit) as se:
        khmer.kfile.check_space_for_graph(outfile, space_needed, force=False,
                                          _testhook_free_space=10e9)
    assert 'ERROR: Not enough free space' in str(se)

    # Now, try with insufficient space, but in force mode
    khmer.kfile.check_space_for_graph(outfile, space_needed, force=True,
                                      _testhook_free_space=10e9)
    assert 'WARNING: Not enough free space' in sys.stderr.getvalue()

    # Finally, try with sufficient space
    sys.stderr = StringIO()
    khmer.kfile.check_space_for_graph(outfile, space_needed, force=False,
                                      _testhook_free_space=20e9)
    assert sys.stderr.getvalue() == ''
    sys.stderr = oldstderr


@pytest.mark.parametrize('graph_type,exp_buckets', [
    ('qfcounttable', '2.4 million buckets'),
    ('countgraph', '3.0 million buckets'),
    ('smallcountgraph', '6.0 million buckets'),
    ('nodegraph', '24.0 million buckets'),
])
def test_check_tablespace_nodegraph(graph_type, exp_buckets):
    parser = khmer_args.build_counting_args()
    args = parser.parse_args(['-M', '3G'])
    buckets_per_table = khmer_args.calculate_graphsize(args, graph_type)
    total_buckets = buckets_per_table * args.n_tables
    sizestr = '{:.1f} million buckets'.format(float(total_buckets) / 1e9)
    assert sizestr == exp_buckets


def test_normal_help(capsys):
    # check -x and -N are hidden by default with --help
    parser = khmer_args.build_graph_args()

    with pytest.raises(SystemExit):
        parser.parse_args(['-h'])

    out, err = capsys.readouterr()
    assert "--max-tablesize" not in out
    assert '--n_tables' not in out


def test_expert_help(capsys):
    # check -x and -N are hidden by default but appear with --help-expert
    old_argv = sys.argv[:]
    sys.argv.append('--help-expert')
    parser = khmer_args.build_graph_args()

    with pytest.raises(SystemExit):
        parser.parse_args(['-h', '--help-expert'])

    out, err = capsys.readouterr()
    assert "--max-tablesize" in out
    assert '--n_tables' in out

    sys.argv = old_argv


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
        tablesize = khmer_args.calculate_graphsize(args, 'countgraph')
        khmer.kfile.check_space_for_graph(outfile, tablesize,
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


def test_check_valid_stdin_nowarn():
    save_stderr, sys.stderr = sys.stderr, io.StringIO()
    try:
        khmer.kfile.check_valid_file_exists(["-"])
        err = sys.stderr.getvalue()
        assert err.count("\n") == 0, err
    except SystemExit as e:
        print(str(e))
    finally:
        sys.stderr = save_stderr


FakeArgparseObject = collections.namedtuple('FakeArgs',
                                            ['ksize', 'n_tables',
                                             'max_tablesize',
                                             'max_memory_usage',
                                             'unique_kmers',
                                             'small_count',
                                             'force'])


def test_create_countgraph_1():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

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

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    countgraph = khmer_args.create_countgraph(args, ksize=15)
    assert countgraph.ksize() == 15


def test_create_countgraph_3():
    # tests too-big ksize

    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    old_stderr = sys.stderr
    sys.stderr = capture = StringIO()

    try:
        khmer_args.create_countgraph(args, ksize=35)
        assert 0, "should not reach this"
    except SystemExit:
        err = capture.getvalue()
        assert 'khmer only supports k-mer sizes <= 32.' in err, err
    finally:
        sys.stderr = old_stderr


def test_create_countgraph_4():
    # tests too-big n_tables WITHOUT force

    ksize = khmer_args.DEFAULT_K
    n_tables = 21  # some number larger than 20
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    old_stderr = sys.stderr
    sys.stderr = capture = StringIO()

    try:
        khmer_args.create_countgraph(args, ksize=None)
        assert 0, "should not reach this"
    except SystemExit:
        err = capture.getvalue()
        assert 'khmer only supports number of tables <= 20.' in err, err
    finally:
        sys.stderr = old_stderr


def test_create_countgraph_5():
    # tests too-big n_tables WITH force

    ksize = khmer_args.DEFAULT_K
    n_tables = 21  # some number larger than 20
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 1)

    old_stderr = sys.stderr
    sys.stderr = capture = StringIO()

    try:
        khmer_args.create_countgraph(args, ksize=None)
        message = "Warning: Maximum recommended number of tables is 20, " + \
                  "discarded by force nonetheless!"
        assert message in capture.getvalue()
    except SystemExit as e:
        print(str(e))
    finally:
        sys.stderr = old_stderr


def test_create_countgraph_4_multiplier():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    countgraph = khmer_args.create_countgraph(args, multiplier=2.0)
    assert sum(countgraph.hashsizes()) < max_mem * 2.0, \
        sum(countgraph.hashsizes())


def test_create_nodegraph_1():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

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

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    nodegraph = khmer_args.create_nodegraph(args, ksize=15)
    assert nodegraph.ksize() == 15


def test_create_nodegraph_3():
    # tests too-big ksize

    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    sys.stderr = capture = StringIO()

    try:
        khmer_args.create_nodegraph(args, ksize=35)
        assert 0, "should not reach this"
    except SystemExit:
        err = capture.getvalue()
        assert 'khmer only supports k-mer sizes <= 32.' in err, err


def test_create_nodegraph_4():
    # tests too-big number of tables WITHOUT force

    ksize = khmer_args.DEFAULT_K
    n_tables = 21  # some number larger than 20
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    sys.stderr = capture = StringIO()

    try:
        khmer_args.create_nodegraph(args, ksize=None)
        assert 0, "should not reach this"
    except SystemExit:
        err = capture.getvalue()
        assert 'khmer only supports number of tables <= 20.' in err, err


def test_create_nodegraph_5():
    # tests too-big number of tables WITH force

    ksize = khmer_args.DEFAULT_K
    n_tables = 21  # some number larger than 20
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 1)

    sys.stderr = capture = StringIO()

    try:
        khmer_args.create_nodegraph(args, ksize=None)
        message = "Warning: Maximum recommended number of tables is 20, " + \
                  "discarded by force nonetheless!"
        assert message in capture.getvalue()
    except SystemExit as e:
        print(str(e))


def test_create_nodegraph_4_multiplier():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    nodegraph = khmer_args.create_nodegraph(args, multiplier=2.0)
    assert sum(nodegraph.hashsizes()) / 8.0 < max_mem * 2.0, \
        sum(nodegraph.hashsizes())


def test_report_on_config_bad_graphtype():
    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    try:
        khmer_args.report_on_config(args, 'foograph')
        assert 0, "the previous statement should raise an exception"
    except ValueError as err:
        assert "unknown graph type: foograph" in str(err), str(err)


def test_fail_calculate_foograph_size():
    # tests unknown graph type

    ksize = khmer_args.DEFAULT_K
    n_tables = khmer_args.DEFAULT_N_TABLES
    max_tablesize = khmer_args.DEFAULT_MAX_TABLESIZE
    max_mem = 1e7

    args = FakeArgparseObject(ksize, n_tables, max_tablesize, max_mem, 0,
                              False, 0)

    try:
        khmer_args.calculate_graphsize(args, 'foograph')
        assert 0, "previous statement should fail"
    except ValueError as err:
        assert "unknown graph type: foograph" in str(err), str(err)


def test_memory_setting():
    assert khmer_args.memory_setting('1') == 1.0
    assert khmer_args.memory_setting('42') == 42.0
    assert khmer_args.memory_setting('10000') == 1e4
    assert khmer_args.memory_setting('2.3e5') == 230000.0
    assert khmer_args.memory_setting('1e9') == 1e9
    assert khmer_args.memory_setting('1K') == 1e3
    assert khmer_args.memory_setting('3.14m') == 3.14e6
    assert khmer_args.memory_setting('8G') == 8e9
    assert khmer_args.memory_setting('8g') == 8e9
    assert khmer_args.memory_setting('16T') == 16e12
    try:
        _ = khmer_args.memory_setting('16Tb')
        assert False, 'previous command should have failed'
    except ValueError as err:
        assert 'cannot parse memory setting' in str(err)
    try:
        _ = khmer_args.memory_setting('16E')
        assert False, 'previous command should have failed'
    except ValueError as err:
        assert 'cannot parse memory setting' in str(err)
    try:
        _ = khmer_args.memory_setting('16Ki')
        assert False, 'previous command should have failed'
    except ValueError as err:
        assert 'cannot parse memory setting' in str(err)
    try:
        _ = khmer_args.memory_setting('b0gu$G')
        assert False, 'previous command should have failed'
    except ValueError as err:
        assert 'cannot parse memory setting' in str(err)
