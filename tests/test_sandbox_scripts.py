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


# sandbox script tests are only run when the tests are loaded from
# the repository
IN_REPOSITORY = os.path.exists(os.path.join(os.path.dirname(__file__),
                                            "../sandbox"))


def _sandbox_scripts():
    if not IN_REPOSITORY:
        return []

    sandbox_path = os.path.join(os.path.dirname(__file__), "../sandbox")
    path = os.path.join(sandbox_path, "*.py")
    return [os.path.normpath(s) for s in glob.glob(path)]


@pytest.mark.parametrize("filename", _sandbox_scripts())
def test_import_succeeds(filename, tmpdir):
    try:
        mod = imp.load_source('__zzz', filename)
    except:
        print(traceback.format_exc())
        raise AssertionError("%s cannot be imported" % (filename,))

    with tmpdir.as_cwd():
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


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_sweep_reads():
    readfile = utils.copy_test_data('test-sweep-reads.fa')
    contigfile = utils.copy_test_data('test-sweep-contigs.fp')
    in_dir = os.path.dirname(contigfile)

    script = scriptpath('sweep-reads.py')
    args = ['-k', '25', '--prefix', 'test', '--label-by-pid',
            contigfile, readfile, 'junkfile.fa']

    status, out, err = utils.runscript(
        script, args, in_dir, sandbox=True)

    # check if the bad file was skipped without issue
    assert 'ERROR' in err, err
    assert 'skipping' in err, err

    out1 = os.path.join(in_dir, 'test_0.fa')
    out2 = os.path.join(in_dir, 'test_1.fa')
    mout = os.path.join(in_dir, 'test_multi.fa')
    oout = os.path.join(in_dir, 'test_orphaned.fa')

    print(os.listdir(in_dir))

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert os.path.exists(mout)
    assert os.path.exists(oout)
    seqs1 = set([r.name for r in screed.open(out1)])
    seqs2 = set([r.name for r in screed.open(out2)])
    seqsm = set([r.name for r in screed.open(mout)])
    seqso = set([r.name for r in screed.open(oout)])

    print(seqs1)
    print(seqs2)
    print(seqsm)
    print(seqso)
    assert seqs1 == set(['read1_p0\t0', 'read2_p0\t0'])
    assert seqs2 == set(['read3_p1\t1'])
    assert (seqsm == set(['read4_multi\t0\t1']) or
            seqsm == set(['read4_multi\t1\t0']))
    assert seqso == set(['read5_orphan'])


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_sweep_reads_fq():
    readfile = utils.copy_test_data('test-sweep-reads.fq')
    contigfile = utils.copy_test_data('test-sweep-contigs.fp')
    in_dir = os.path.dirname(contigfile)

    script = scriptpath('sweep-reads.py')
    args = ['-k', '25', '--prefix', 'test', '--label-by-pid',
            contigfile, readfile, 'junkfile.fa']

    status, out, err = utils.runscript(
        script, args, in_dir, sandbox=True)

    # check if the bad file was skipped without issue
    assert 'ERROR' in err, err
    assert 'skipping' in err, err

    out1 = os.path.join(in_dir, 'test_0.fq')
    out2 = os.path.join(in_dir, 'test_1.fq')
    mout = os.path.join(in_dir, 'test_multi.fq')
    oout = os.path.join(in_dir, 'test_orphaned.fq')

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert os.path.exists(mout)
    assert os.path.exists(oout)
    print(open(out1).read())

    print(os.listdir(in_dir))

    seqs1 = set([r.name for r in screed.open(out1)])
    seqs2 = set([r.name for r in screed.open(out2)])
    seqsm = set([r.name for r in screed.open(mout)])
    seqso = set([r.name for r in screed.open(oout)])

    print(seqs1)
    print(seqs2)
    print(seqsm)
    print(seqso)
    assert seqs1 == set(['read1_p0\t0', 'read2_p0\t0'])
    assert seqs2 == set(['read3_p1\t1'])
    assert (seqsm == set(['read4_multi\t0\t1']) or
            seqsm == set(['read4_multi\t1\t0']))
    assert seqso == set(['read5_orphan'])

    seqs1 = set([r.quality for r in screed.open(out1)])
    seqs2 = set([r.quality for r in screed.open(out2)])
    seqsm = set([r.quality for r in screed.open(mout)])
    seqso = set([r.quality for r in screed.open(oout)])


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_sweep_reads_2():

    infile = utils.copy_test_data('random-20-X2.fa')
    inref = utils.copy_test_data('random-20-a.fa')

    wdir = os.path.dirname(inref)
    script = scriptpath('sweep-reads.py')
    args = ['-m', '50', '-k', '20', '-l', '9', '-b', '60', '--prefix',
            'test', '--label-by-seq', inref, infile]
    status, out, err = utils.runscript(script, args, wdir, sandbox=True)

    for i in range(99):
        p = os.path.join(wdir, 'test_{i}.fa'.format(i=i))
        print(p, err, out)
        assert os.path.exists(p)
        os.remove(p)
    assert os.path.exists(os.path.join(wdir, 'test.counts.csv'))
    assert os.path.exists(os.path.join(wdir, 'test.dist.txt'))
    assert not os.path.exists(os.path.join(wdir, 'test_multi.fa'))


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_sweep_reads_3():

    infile = utils.copy_test_data('random-20-a.fa')
    wdir = os.path.dirname(infile)
    script = scriptpath('sweep-reads.py')
    args = ['-m', '75', '-k', '20', '-l', '1', '--prefix',
            'test', '--label-by-group', '10', infile, infile]
    status, out, err = utils.runscript(script, args, wdir, sandbox=True)

    for i in range(10):
        p = os.path.join(wdir, 'test_{i}.fa'.format(i=i))
        print(p, err, out)
        assert os.path.exists(p)
        os.remove(p)

    counts_fn = os.path.join(wdir, 'test.counts.csv')
    with open(counts_fn) as cfp:
        for line in cfp:
            _, _, c = line.partition(',')
            assert int(c) in [9, 10]

    assert os.path.exists(counts_fn)
    assert os.path.exists(os.path.join(wdir, 'test.dist.txt'))
    assert not os.path.exists(os.path.join(wdir, 'test_multi.fa'))


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_collect_reads():
    outfile = utils.get_temp_filename('out.graph')
    infile = utils.get_test_data('test-reads.fa')
    script = 'collect-reads.py'
    args = ['-M', '1e7', outfile, infile]

    status, out, err = utils.runscript(script, args, sandbox=True)

    assert status == 0
    assert os.path.exists(outfile)


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_saturate_by_median():
    infile = utils.get_test_data('test-reads.fa')
    script = 'saturate-by-median.py'
    args = ['-M', '1e7', infile]

    status, out, err = utils.runscript(script, args, sandbox=True)

    assert status == 0


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_count_kmers_1():
    infile = utils.copy_test_data('random-20-a.fa')
    ctfile = _make_counting(infile)

    script = scriptpath('count-kmers.py')
    args = [ctfile, infile]

    status, out, err = utils.runscript(script, args, os.path.dirname(infile),
                                       sandbox=True)

    out = out.splitlines()
    assert 'TTGTAACCTGTGTGGGGTCG,1' in out


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_count_kmers_2_single():
    infile = utils.copy_test_data('random-20-a.fa')

    script = scriptpath('count-kmers-single.py')
    args = ['-x', '1e7', '-k', '20', '-N', '2', infile]

    status, out, err = utils.runscript(script, args, os.path.dirname(infile),
                                       sandbox=True)

    out = out.splitlines()
    assert 'TTGTAACCTGTGTGGGGTCG,1' in out


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_multirename_fasta():
    infile1 = utils.copy_test_data('test-multi.fa')
    multioutfile = utils.get_temp_filename('out.fa')
    infile2 = utils.copy_test_data('multi-output.fa')
    args = ['assembly', infile1]
    _, out, err = utils.runscript('multi-rename.py', args, sandbox=True)
    r = open(infile2).read()
    assert r in out


@pytest.mark.skipif(not IN_REPOSITORY,
                    reason='executing outside of the repository')
def test_extract_compact_dbg_1():
    infile = utils.get_test_data('simple-genome.fa')
    outfile = utils.get_temp_filename('out.gml')
    args = ['-x', '1e4', '-o', outfile, infile]
    _, out, err = utils.runscript('extract-compact-dbg.py', args, sandbox=True)

    print(out)
    print(err)
    assert os.path.exists(outfile)

    assert '174 segments, containing 2803 nodes' in out
