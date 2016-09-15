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

# pylint: disable=C0111,C0103,E1103,W0612

from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
import os
import os.path
from io import StringIO
import traceback
import glob
import imp

import pytest

from . import khmer_tst_utils as utils
import screed
from .test_scripts import _make_counting


def scriptpath(script):
    return script


def teardown():
    utils.cleanup()


def _sandbox_scripts():
    sandbox_path = os.path.join(os.path.dirname(__file__), "../sandbox")
    if not os.path.exists(sandbox_path):
        pytest.skip("sandbox scripts are only tested in a repository")

    path = os.path.join(sandbox_path, "*.py")
    return [os.path.normpath(s) for s in glob.glob(path)]


@pytest.mark.parametrize("filename", _sandbox_scripts())
def test_import_succeeds(filename):
    try:
        mod = imp.load_source('__zzz', filename)
    except:
        print(traceback.format_exc())
        raise AssertionError("%s cannot be imported" % (filename,))

    oldargs = sys.argv
    sys.argv = [filename]

    oldout, olderr = sys.stdout, sys.stderr
    sys.stdout = StringIO()
    sys.stderr = StringIO()

    try:
        try:
            global_dict = {'__name__': '__main__'}
            exec(  # pylint: disable=exec-used
                compile(open(filename).read(), filename, 'exec'),
                global_dict)
        except (ImportError, SyntaxError) as err:
            print("{0}".format(err))
            raise AssertionError("%s cannot be exec'd" % (filename),
                                 "{0}".format(traceback))
        except:  # pylint: disable=bare-except
            pass                        # other failures are expected :)
    finally:
        sys.argv = oldargs
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdout, sys.stderr = oldout, olderr


def test_collect_reads():
    outfile = utils.get_temp_filename('out.graph')
    infile = utils.get_test_data('test-reads.fa')
    script = 'collect-reads.py'
    args = ['-M', '1e7', outfile, infile]

    status, out, err = utils.runscript(script, args, sandbox=True)

    assert status == 0
    assert os.path.exists(outfile)


def test_saturate_by_median():
    infile = utils.get_test_data('test-reads.fa')
    script = 'saturate-by-median.py'
    args = ['-M', '1e7', infile]

    status, out, err = utils.runscript(script, args, sandbox=True)

    assert status == 0


def test_count_kmers_1():
    infile = utils.copy_test_data('random-20-a.fa')
    ctfile = _make_counting(infile)

    script = scriptpath('count-kmers.py')
    args = [ctfile, infile]

    status, out, err = utils.runscript(script, args, os.path.dirname(infile),
                                       sandbox=True)

    out = out.splitlines()
    assert 'TTGTAACCTGTGTGGGGTCG,1' in out


def test_count_kmers_2_single():
    infile = utils.copy_test_data('random-20-a.fa')

    script = scriptpath('count-kmers-single.py')
    args = ['-x', '1e7', '-k', '20', '-N', '2', infile]

    status, out, err = utils.runscript(script, args, os.path.dirname(infile),
                                       sandbox=True)

    out = out.splitlines()
    assert 'TTGTAACCTGTGTGGGGTCG,1' in out


def test_multirename_fasta():
    infile1 = utils.copy_test_data('test-multi.fa')
    multioutfile = utils.get_temp_filename('out.fa')
    infile2 = utils.copy_test_data('multi-output.fa')
    args = ['assembly', infile1]
    _, out, err = utils.runscript('multi-rename.py', args, sandbox=True)
    r = open(infile2).read()
    assert r in out


def test_extract_compact_dbg_1():
    infile = utils.get_test_data('simple-genome.fa')
    outfile = utils.get_temp_filename('out.gml')
    args = ['-x', '1e4', '-o', outfile, infile]
    _, out, err = utils.runscript('extract-compact-dbg.py', args, sandbox=True)

    print(out)
    print(err)
    assert os.path.exists(outfile)

    assert '174 segments, containing 2803 nodes' in out
