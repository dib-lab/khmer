#
# This file is part of khmer, http://github.com/ged-lab/khmer/, and is
# Copyright (C) Michigan State University, 2014. It is licensed under
# the three-clause BSD license; see doc/LICENSE.txt.
# Contact: khmer-project@idyll.org
#

# pylint: disable=C0111,C0103,E1103,W0612

import sys
import os
import os.path
import shutil
from cStringIO import StringIO
import traceback
import nose
import glob
import imp

import khmer_tst_utils as utils
import khmer
import khmer.file
import screed


def scriptpath(script):
    return script


def teardown():
    utils.cleanup()


def test_import_all():
    sandbox_path = os.path.join(os.path.dirname(__file__), "../sandbox")
    if not os.path.exists(sandbox_path):
        raise nose.SkipTest("sandbox scripts are only tested in a repository")

    path = os.path.join(sandbox_path, "*.py")
    scripts = glob.glob(path)
    for s in scripts:
        s = os.path.normpath(s)
        yield _checkImportSucceeds('test_sandbox_scripts.py', s)


class _checkImportSucceeds(object):
    def __init__(self, tag, filename):
        self.tag = tag
        self.filename = filename
        self.description = '%s: test import %s' % (self.tag,
                                                   os.path.split(filename)[-1])

    def __call__(self):
        try:
            mod = imp.load_source('__zzz', self.filename)
        except:
            print traceback.format_exc()
            raise AssertionError("%s cannot be imported" % (self.filename,))

        ###

        oldargs = sys.argv
        sys.argv = [self.filename]

        oldout, olderr = sys.stdout, sys.stderr
        sys.stdout = StringIO()
        sys.stderr = StringIO()

        try:
            try:
                global_dict = {'__name__': '__main__'}
                execfile(self.filename, global_dict)
            except (ImportError, SyntaxError):
                print traceback.format_exc()
                raise AssertionError("%s cannot be exec'd" % (self.filename,))
            except:
                pass                        # other failures are expected :)
        finally:
            sys.argv = oldargs
            out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
            sys.stdout, sys.stderr = oldout, olderr


def test_sweep_reads():
    readfile = utils.get_temp_filename('reads.fa')
    contigfile = utils.get_temp_filename('contigs.fp')
    in_dir = os.path.dirname(contigfile)

    shutil.copyfile(utils.get_test_data('test-sweep-reads.fa'), readfile)
    shutil.copyfile(utils.get_test_data('test-sweep-contigs.fp'), contigfile)

    script = scriptpath('sweep-reads.py')
    args = ['-k', '25', '--prefix', 'test', '--label-by-pid',
            contigfile, readfile, 'junkfile.fa']

    status, out, err = utils.runscript(
        script, args, in_dir, fail_ok=True, sandbox=True)

    # check if the bad file was skipped without issue
    assert 'ERROR' in err, err
    assert 'skipping' in err, err

    out1 = os.path.join(in_dir, 'test_0.fa')
    out2 = os.path.join(in_dir, 'test_1.fa')
    mout = os.path.join(in_dir, 'test_multi.fa')
    oout = os.path.join(in_dir, 'test_orphaned.fa')

    print os.listdir(in_dir)

    seqs1 = set([r.name for r in screed.open(out1)])
    seqs2 = set([r.name for r in screed.open(out2)])
    seqsm = set([r.name for r in screed.open(mout)])
    seqso = set([r.name for r in screed.open(oout)])

    print seqs1
    print seqs2
    print seqsm
    print seqso
    assert seqs1 == set(['read1_p0\t0', 'read2_p0\t0'])
    assert seqs2 == set(['read3_p1\t1'])
    assert (seqsm == set(['read4_multi\t0\t1']) or
            seqsm == set(['read4_multi\t1\t0']))
    assert seqso == set(['read5_orphan'])


def test_sweep_reads_fq():
    readfile = utils.get_temp_filename('reads.fa')
    contigfile = utils.get_temp_filename('contigs.fp')
    in_dir = os.path.dirname(contigfile)

    shutil.copyfile(utils.get_test_data('test-sweep-reads.fq'), readfile)
    shutil.copyfile(utils.get_test_data('test-sweep-contigs.fp'), contigfile)

    script = scriptpath('sweep-reads.py')
    args = ['-k', '25', '--prefix', 'test', '--label-by-pid',
            contigfile, readfile, 'junkfile.fa']

    status, out, err = utils.runscript(
        script, args, in_dir, fail_ok=True, sandbox=True)

    # check if the bad file was skipped without issue
    assert 'ERROR' in err, err
    assert 'skipping' in err, err

    out1 = os.path.join(in_dir, 'test_0.fq')
    out2 = os.path.join(in_dir, 'test_1.fq')
    mout = os.path.join(in_dir, 'test_multi.fq')
    oout = os.path.join(in_dir, 'test_orphaned.fq')

    print open(out1).read()

    print os.listdir(in_dir)

    seqs1 = set([r.name for r in screed.open(out1)])
    seqs2 = set([r.name for r in screed.open(out2)])
    seqsm = set([r.name for r in screed.open(mout)])
    seqso = set([r.name for r in screed.open(oout)])

    print seqs1
    print seqs2
    print seqsm
    print seqso
    assert seqs1 == set(['read1_p0\t0', 'read2_p0\t0'])
    assert seqs2 == set(['read3_p1\t1'])
    assert (seqsm == set(['read4_multi\t0\t1']) or
            seqsm == set(['read4_multi\t1\t0']))
    assert seqso == set(['read5_orphan'])

    seqs1 = set([r.accuracy for r in screed.open(out1)])
    seqs2 = set([r.accuracy for r in screed.open(out2)])
    seqsm = set([r.accuracy for r in screed.open(mout)])
    seqso = set([r.accuracy for r in screed.open(oout)])


def test_sweep_reads_2():

    infile = utils.get_temp_filename('seqs.fa')
    inref = utils.get_temp_filename('ref.fa')
    shutil.copyfile(utils.get_test_data('random-20-X2.fa'), infile)
    shutil.copyfile(utils.get_test_data('random-20-a.fa'), inref)
    wdir = os.path.dirname(inref)
    script = scriptpath('sweep-reads.py')
    args = ['-m', '50', '-k', '20', '-l', '9', '-b', '60', '--prefix',
            'test', '--label-by-seq', inref, infile]
    status, out, err = utils.runscript(script, args, wdir, sandbox=True)

    for i in xrange(99):
        p = os.path.join(wdir, 'test_{i}.fa'.format(i=i))
        print p, err, out
        assert os.path.exists(p)
        os.remove(p)
    assert os.path.exists(os.path.join(wdir, 'test.counts.csv'))
    assert os.path.exists(os.path.join(wdir, 'test.dist.txt'))
    assert not os.path.exists(os.path.join(wdir, 'test_multi.fa'))


def test_sweep_reads_3():

    infile = utils.get_temp_filename('seqs.fa')
    shutil.copyfile(utils.get_test_data('random-20-a.fa'), infile)
    wdir = os.path.dirname(infile)
    script = scriptpath('sweep-reads.py')
    args = ['-m', '75', '-k', '20', '-l', '1', '--prefix',
            'test', '--label-by-group', '10', infile, infile]
    status, out, err = utils.runscript(script, args, wdir, sandbox=True)

    for i in xrange(10):
        p = os.path.join(wdir, 'test_{i}.fa'.format(i=i))
        print p, err, out
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


def test_trim_low_abund_1():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", infile]
    utils.runscript('trim-low-abund.py', args, in_dir, sandbox=True)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_trim_low_abund_1_duplicate_filename_err():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-C', '1', infile, infile]
    try:
        utils.runscript('trim-low-abund.py', args, in_dir, sandbox=True)
        raise Exception("should not reach this")
    except AssertionError:
        # an error should be raised by passing 'infile' twice.
        pass


def test_trim_low_abund_2():
    infile = utils.get_temp_filename('test.fa')
    infile2 = utils.get_temp_filename('test2.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile2)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-C', '1', infile, infile2]
    utils.runscript('trim-low-abund.py', args, in_dir, sandbox=True)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

# make sure that FASTQ records are retained.


def test_trim_low_abund_3_fq_retained():
    infile = utils.get_temp_filename('test.fq')
    infile2 = utils.get_temp_filename('test2.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile)
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile2)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-C', '1', infile, infile2]
    utils.runscript('trim-low-abund.py', args, in_dir, sandbox=True)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

    # check for 'accuracy' string.
    seqs = set([r.accuracy for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert '##################' in seqs


# test that the -V option does not trim sequences that are low abundance


def test_trim_low_abund_4_retain_low_abund():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-V', infile]
    utils.runscript('trim-low-abund.py', args, in_dir, sandbox=True)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

# test that the -V option *does* trim sequences that are low abundance


def test_trim_low_abund_5_trim_high_abund():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-3.fa'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-V', infile]
    utils.runscript('trim-low-abund.py', args, in_dir, sandbox=True)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs

    # trimmed sequence @ error
    assert 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGC' in seqs

# test that -V/-Z setting - should not trip if -Z is set high enough.


def test_trim_low_abund_6_trim_high_abund_Z():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-3.fa'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-V', '-Z', '25', infile]
    utils.runscript('trim-low-abund.py', args, in_dir, sandbox=True)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs

    # untrimmed seq.
    badseq = 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCgtgCCGCAGCTGTCGTCAGGG' \
             'GATTTCCGGGCGG'
    assert badseq in seqs       # should be there, untrimmed
