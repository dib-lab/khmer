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
from . import khmer_tst_utils as utils

import khmer.kfile


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
    try:
        khmer.kfile.check_space_for_hashtable(
            1e9, force=False, _testhook_free_space=0)
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
    try:
        khmer.kfile.check_space_for_hashtable(
            1e9, force=True, _testhook_free_space=0)
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
