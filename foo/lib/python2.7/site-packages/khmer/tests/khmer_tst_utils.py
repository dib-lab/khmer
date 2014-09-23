#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2013. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#
import tempfile
import os
import shutil
from pkg_resources import Requirement, resource_filename, ResolutionError


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(
            Requirement.parse("khmer"), "khmer/tests/test-data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'test-data',
                                filename)
    return filepath

cleanup_list = []


def get_temp_filename(filename, tempdir=None):
    if tempdir is None:
        tempdir = tempfile.mkdtemp(prefix='khmertest_')
        cleanup_list.append(tempdir)

    return os.path.join(tempdir, filename)


def cleanup():
    global cleanup_list

    for path in cleanup_list:
        shutil.rmtree(path, ignore_errors=True)
    cleanup_list = []
