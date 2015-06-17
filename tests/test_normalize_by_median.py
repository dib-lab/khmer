from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#

import os
import shutil
import threading
import io

import screed
import khmer

from . import khmer_tst_utils as utils
from nose.plugins.attrib import attr
from .test_scripts import _make_counting


def test_normalize_by_median_indent():
    infile = utils.get_test_data('paired-mixed.fa.pe')
    hashfile = utils.get_test_data('normC20k20.ct')
    outfile = utils.get_temp_filename('paired-mixed.fa.pe.keep')
    script = 'normalize-by-median.py'
    args = ['--loadtable', hashfile, '-o', outfile, infile]
    (status, out, err) = utils.runscript(script, args)
    assert status == 0, (out, err)
    assert os.path.exists(outfile)


def test_normalize_by_median():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 98' in err, err

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [r.sequence for r in screed.open(outfile)]
    assert len(seqs) == 1, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs
    assert "IOErrors" not in err


def test_normalize_by_median_unpaired_final_read():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('single-read.fq'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', '-p',  infile]
    try:
        (status, out, err) = utils.runscript(script, args, in_dir)
        raise Exception("Shouldn't get to this")
    except AssertionError as e:
        out = str(e)
        assert "ERROR: Unpaired reads when require_paired" in out, out


def test_normalize_by_median_unforced_badfile():
    CUTOFF = '1'

    infile = utils.get_temp_filename("potatoes")
    outfile = infile + '.keep'
    in_dir = os.path.dirname(infile)
    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', infile]
    try:
        (status, out, err) = utils.runscript(script, args, in_dir)
        raise Exception("Shouldn't get to this")
    except AssertionError as e:
        out = str(e)
        assert "ERROR: [Errno 2] No such file or directory:" in out, out

    if os.path.exists(outfile):
        assert False, '.keep file should have been removed: '


def test_normalize_by_median_contradictory_args():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    outfile = utils.get_temp_filename('report.out')

    shutil.copyfile(utils.get_test_data('test-large.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', '1', '-k', '17', '--force-single', '-p', '-R',
            outfile, infile]
    try:
        (status, out, err) = utils.runscript(script, args, in_dir)
        raise Exception("Shouldn't get to this")
    except AssertionError as e:
        out = str(e)
        assert "cannot both be set" in out, out


def test_normalize_by_median_stdout_3():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', infile, '--out', '-']
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 98' in err, err
    assert 'in /dev/stdout' in err, err
    assert "IOErrors" not in err


@attr('known_failing')
def test_normalize_by_median_known_good():
    CUTOFF = '2'

    infile = utils.get_temp_filename('test.fa.gz')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data('100k-filtered.fa.gz'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '20', '-x', '4e6', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile
    iter_known = screed.open(utils.get_test_data('100k-filtered.fa.keep.gz'))
    iter_out = screed.open(outfile)
    try:
        for rknown, rout in zip(iter_known, iter_out):
            assert rknown.name == rout.name
    except Exception as e:
        print(e)
        assert False


def test_normalize_by_median_report_fp():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    outfile = utils.get_temp_filename('report.out')

    shutil.copyfile(utils.get_test_data('test-large.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', '1', '-k', '17', '-R', outfile, infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert "fp rate estimated to be 0.626" in err, err
    report = open(outfile, 'r')
    line = report.readline()
    assert "100000 25232 0.25232" in line, line


def test_normalize_by_median_unpaired_and_paired():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-paired.fa'), infile)

    unpairedfile = utils.get_temp_filename('test1.fa', tempdir=in_dir)
    shutil.copyfile(utils.get_test_data('random-20-a.fa'), unpairedfile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', '-u', unpairedfile, '-p', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 4029' in err, err

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile


def test_normalize_by_median_count_kmers_PE():
    CUTOFF = '1'
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    # The test file has one pair of identical read except for the last base
    # The 2nd read should be discarded in the unpaired mode
    # but kept in the paired end mode adding only one more unique kmer
    shutil.copyfile(utils.get_test_data('paired_one.base.dif.fa'), infile)
    script = 'normalize-by-median.py'

    args = ['-C', CUTOFF, '-k', '17', '--force-single', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert 'Total number of unique k-mers: 98' in err, err
    assert 'kept 1 of 2 or 50%' in err, err

    args = ['-C', CUTOFF, '-k', '17', '-p', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert 'Total number of unique k-mers: 99' in err, err
    assert 'kept 2 of 2 or 100%' in err, err


def test_normalize_by_median_double_file_name():
    infile = utils.get_temp_filename('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'normalize-by-median.py'
    args = [utils.get_test_data('test-abund-read-2.fa'), infile]

    try:
        (status, out, err) = utils.runscript(script, args, in_dir)
    except AssertionError as e:
        assert "Duplicate filename--Cannot handle this!" in str(e), str(e)


def test_normalize_by_median_overwrite():
    outfile = utils.get_temp_filename('test.fa.keep')
    shutil.copyfile(utils.get_test_data('test-abund-read.fa'), outfile)
    in_dir = os.path.dirname(outfile)

    CUTOFF = '1'
    infile = utils.get_temp_filename('test.fa', in_dir)
    shutil.copyfile(utils.get_test_data('test-abund-read-3.fa'), infile)
    script = 'normalize-by-median.py'

    args = ['-C', CUTOFF, '-k', '17', '-o', outfile, infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert os.path.exists(outfile), outfile
    seqs = [r.sequence for r in screed.open(outfile)]
    assert len(seqs) == 1, seqs
    assert 'GACAGCgtgCCGCA' in seqs[0], seqs


def test_normalize_by_median_version():
    script = 'normalize-by-median.py'
    args = ['--version']
    status, out, err = utils.runscript(script, args)

    errlines = err.splitlines()
    for err in errlines:
        if err.startswith('||') or \
           not err.strip():
            continue
        break

    print(errlines)
    print(err)

    assert err.startswith('khmer ')


def test_normalize_by_median_2():
    CUTOFF = '2'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [r.sequence for r in screed.open(outfile)]
    assert len(seqs) == 2, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs
    assert seqs[1] == 'GGTTGACGGGGCTCAGGG', seqs


def test_normalize_by_median_paired():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-paired.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-p', '-k', '17', infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [r.sequence for r in screed.open(outfile)]
    assert len(seqs) == 2, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs
    assert seqs[1].startswith('GGTTGACGGGGCTCAGGG'), seqs


def test_normalize_by_median_paired_fq():
    CUTOFF = '20'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-paired.fq'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-p', '-k', '17', infile]
    _, out, err = utils.runscript(script, args, in_dir)
    print(out)
    print(err)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [r.sequence for r in screed.open(outfile)]
    assert len(seqs) == 6, len(seqs)
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs
    assert seqs[1].startswith('GGTTGACGGGGCTCAGGG'), seqs

    names = [r.name for r in screed.open(outfile, parse_description=False)]
    assert len(names) == 6, names
    assert '895:1:37:17593:9954 1::FOO' in names, names
    assert '895:1:37:17593:9954 2::FOO' in names, names


def test_normalize_by_median_impaired():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-impaired.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-p', '-k', '17', infile]
    _, out, err = utils.runscript(script, args, in_dir, fail_ok=True)
    assert 'ERROR: Unpaired reads ' in err, err


def test_normalize_by_median_force():
    CUTOFF = '1'

    corrupt_infile = utils.get_temp_filename('test-corrupt.fq')
    good_infile = utils.get_temp_filename('test-good.fq',
                                          tempdir=os.path.dirname(
                                              corrupt_infile))

    in_dir = os.path.dirname(good_infile)

    shutil.copyfile(utils.get_test_data('test-error-reads.fq'), corrupt_infile)
    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), good_infile)

    script = 'normalize-by-median.py'
    args = ['-f', '-C', CUTOFF, '-k', '17', corrupt_infile, good_infile]

    (status, out, err) = utils.runscript(script, args, in_dir)

    assert '*** Skipping' in err
    assert '** IOErrors' in err


def test_normalize_by_median_no_bigcount():
    infile = utils.get_temp_filename('test.fa')
    hashfile = utils.get_temp_filename('test-out.ct')
    outfile = infile + '.keep'
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = 'normalize-by-median.py'
    args = ['-C', '1000', '-k 8', '--savetable', hashfile, infile]

    (status, out, err) = utils.runscript(script, args, in_dir)
    assert status == 0, (out, err)
    print((out, err))

    assert os.path.exists(hashfile), hashfile
    kh = khmer.load_counting_hash(hashfile)

    assert kh.get('GGTTGACG') == 255


def test_normalize_by_median_empty():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-empty.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile


def test_normalize_by_median_emptycountingtable():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-empty.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '--loadtable', infile, infile]
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert 'ValueError' in err, (status, out, err)


def test_normalize_by_median_fpr():
    MIN_TABLESIZE_PARAM = 1

    infile = utils.get_temp_filename('test-fpr.fq')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), infile)

    script = 'normalize-by-median.py'
    args = ['-f', '-k 17', '-x ' + str(MIN_TABLESIZE_PARAM), infile]

    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)

    assert os.path.exists(infile + '.keep')
    assert '** ERROR: the graph structure is too small' in err, err


def write_by_chunks(infile, outfile, CHUNKSIZE=8192):
    ifile = io.open(infile, 'rb')
    ofile = io.open(outfile, 'wb')
    chunk = ifile.read(CHUNKSIZE)
    while len(chunk) > 0:
        ofile.write(chunk)
        chunk = ifile.read(CHUNKSIZE)
    ifile.close()
    ofile.close()


def test_normalize_by_median_streaming():
    CUTOFF = '20'

    infile = utils.get_test_data('100-reads.fq.gz')
    in_dir = os.path.dirname(infile)
    fifo = utils.get_temp_filename('fifo')
    outfile = utils.get_temp_filename('outfile')

    # Use a fifo to copy stdout to a file for checking
    os.mkfifo(fifo)
    thread = threading.Thread(target=write_by_chunks, args=(fifo, outfile))
    thread.start()

    # Execute diginorm
    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', '-o', fifo, infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    # Merge the thread
    thread.join()

    assert os.path.exists(outfile), outfile
    with open(outfile) as fp:
        linecount = sum(1 for _ in fp)
    assert linecount == 400


def test_diginorm_basic_functionality_1():
    # each of these pairs has both a multicopy sequence ('ACTTCA...') and
    # a random sequence.  With 'C=1' and '-p', all should be kept.
    CUTOFF = ['-C', '1']
    PAIRING = ['-p']

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('dn-test-all-paired-all-keep.fa'),
                    infile)

    script = 'normalize-by-median.py'
    args = list(CUTOFF) + list(PAIRING) + ['-k', '15', infile]
    _, out, err = utils.runscript(script, args, in_dir)
    print(out)
    print(err)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])

    assert seqs == set(['a/1', 'a/2',
                        'b/1', 'b/2',
                        'c/1', 'c/2',
                        'd/1', 'd/2']), seqs


def test_diginorm_basic_functionality_2():
    # each of these pairs has both a multicopy sequence ('ACTTCA...')
    # and a random sequence ('G...').  With 'C=1' and '--force-
    # single', only random seqs should be kept, together with one copy
    # of the multicopy sequence.
    CUTOFF = ['-C', '1']
    PAIRING = ['--force-single']

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('dn-test-all-paired-all-keep.fa'),
                    infile)

    script = 'normalize-by-median.py'
    args = list(CUTOFF) + list(PAIRING) + ['-k', '15', infile]
    _, out, err = utils.runscript(script, args, in_dir)
    print(out)
    print(err)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])

    assert seqs == set(['a/1', 'a/2',
                        'b/2',
                        'c/1',
                        'd/2']), seqs


def test_diginorm_basic_functionality_3():
    # This data is entirely unpaired, but with one duplicate ('A...').
    # and a random sequence ('G...').  With 'C=1' only three seqs should
    # be left, with no other complaints.

    CUTOFF = ['-C', '1']
    PAIRING = []

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('dn-test-none-paired.fa'),
                    infile)

    script = 'normalize-by-median.py'
    args = list(CUTOFF) + list(PAIRING) + ['-k', '15', infile]
    _, out, err = utils.runscript(script, args, in_dir)
    print(out)
    print(err)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])

    assert seqs == set(['a/1',
                        'b/2',
                        'd/1']), seqs


def test_diginorm_basic_functionality_4():
    # This data is mixed paired/unpaired, but with one duplicate ('A...').
    # and a random sequence ('G...').  With 'C=2' all of the sequences
    # should be kept.

    CUTOFF = ['-C', '1']
    PAIRING = ['-p']

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('dn-test-some-paired-all-keep.fa'),
                    infile)

    script = 'normalize-by-median.py'
    args = list(CUTOFF) + list(PAIRING) + ['-k', '15', infile]
    _, out, err = utils.runscript(script, args, in_dir)
    print(out)
    print(err)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])

    assert seqs == set(['a/1', 'a/2',
                        'b/2',
                        'c/1', 'c/2',
                        'd/2']), seqs


def test_diginorm_basic_functionality_4():
    # each of these pairs has both a multicopy sequence ('ACTTCA...') and
    # a random sequence.  With 'C=1' and '-p', all should be
    CUTOFF = ['-C', '1']
    PAIRING = ['-p']

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('dn-test-all-paired-all-keep.fa'),
                    infile)

    script = 'normalize-by-median.py'
    args = list(CUTOFF) + list(PAIRING) + ['-k', '15', infile]
    _, out, err = utils.runscript(script, args, in_dir)
    print(out)
    print(err)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])

    assert seqs == set(['a/1', 'a/2',
                        'b/1', 'b/2',
                        'c/1', 'c/2',
                        'd/1', 'd/2']), seqs
