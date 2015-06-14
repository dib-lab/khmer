from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#

# pylint: disable=C0111,C0103,E1103,W0612

import json
import sys
import os
import stat
import shutil
from io import StringIO
import traceback
from nose.plugins.attrib import attr
import subprocess
import threading
import bz2
import io

from . import khmer_tst_utils as utils
import khmer
import khmer.kfile
import screed


def teardown():
    utils.cleanup()


def test_check_space():
    # @CTB this probably belongs in a new test file, along with other
    # tests of the file.py module.
    khmer.kfile.check_space(
        ['', utils.get_test_data('test-abund-read-2.fa')], False)


def test_load_into_counting():
    script = 'load-into-counting.py'
    args = ['-x', '1e3', '-N', '2', '-k', '20', '-t']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 89' in err, err
    assert os.path.exists(outfile)


def test_load_into_counting_nonwritable():
    script = 'load-into-counting.py'
    args = ['-x', '1e3', '-N', '2', '-k', '20', '-t']

    outfile = utils.get_temp_filename('test-nonwritable')
    with open(outfile, 'w') as fout:
        fout.write("This file is non-writable (after this)")

    os.chmod(outfile, stat.S_IWOTH | stat.S_IRUSR)
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert 'does not have write permission; exiting' in err, err
    assert status == 1, status


@attr('huge')
def test_load_into_counting_toobig():
    script = 'load-into-counting.py'
    args = ['-x', '1e12', '-N', '2', '-k', '20', '-t', '--force']

    outfile = utils.get_temp_filename('out.kh')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == -1, status
    assert "MemoryError" in err, err


def test_load_into_counting_fail():
    script = 'load-into-counting.py'
    args = ['-x', '1e2', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == 1, status
    print(err)
    assert "** ERROR: the graph structure is too small" in err


def test_load_into_counting_multifile():
    script = 'load-into-counting.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-t']

    outfile = utils.get_temp_filename('out.kh')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile, infile, infile, infile, infile,
                 infile, infile, infile, infile, infile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 95' in err, err
    assert os.path.exists(outfile)


def test_load_into_counting_tsv():
    script = 'load-into-counting.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-t', '-s', 'tsv']

    outfile = utils.get_temp_filename('out.ct')
    tabfile = outfile + '.info.tsv'
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 95' in err, err
    assert os.path.exists(outfile)
    assert os.path.exists(tabfile)
    with open(tabfile) as tabfh:
        tabfile_lines = tabfh.readlines()
    assert len(tabfile_lines) == 2
    outbase = os.path.basename(outfile)
    tsv = [outbase, '0.000', '95', '1001', infile]
    expected_tsv_line = '\t'.join(tsv) + '\n'
    assert tabfile_lines[1] == expected_tsv_line, tabfile_lines


def test_load_into_counting_json():
    script = 'load-into-counting.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-t', '-s', 'json']

    outfile = utils.get_temp_filename('out.ct')
    jsonfile = outfile + '.info.json'
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 95' in err, err
    assert os.path.exists(outfile)
    assert os.path.exists(jsonfile)

    with open(jsonfile) as jsonfh:
        got_json = json.load(jsonfh)
    outbase = os.path.basename(outfile)
    expected_json = {
        "files": [infile],
        "ht_name": outbase,
        "num_kmers": 95,
        "num_reads": 1001,
        "fpr": 9.024965705097741e-11,
        "mrinfo_version": "0.2.0",
    }

    assert got_json == expected_json, got_json


def test_load_into_counting_bad_summary_fmt():
    script = 'load-into-counting.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-s', 'badfmt']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status != 0, status
    assert "invalid choice: 'badfmt'" in err, err


def _make_counting(infilename, SIZE=1e7, N=2, K=20, BIGCOUNT=True):
    script = 'load-into-counting.py'
    args = ['-x', str(SIZE), '-N', str(N), '-k', str(K)]

    if not BIGCOUNT:
        args.append('-b')

    outfile = utils.get_temp_filename('out.ct')

    args.extend([outfile, infilename])

    utils.runscript(script, args)
    assert os.path.exists(outfile)

    return outfile


def test_filter_abund_1():
    script = 'filter-abund.py'

    infile = utils.get_temp_filename('test.fa')
    n_infile = utils.get_temp_filename('test-fastq-n-reads.fq')

    in_dir = os.path.dirname(infile)
    n_in_dir = os.path.dirname(n_infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    shutil.copyfile(utils.get_test_data('test-fastq-n-reads.fq'), n_infile)

    counting_ht = _make_counting(infile, K=17)
    n_counting_ht = _make_counting(n_infile, K=17)

    args = [counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    n_outfile = n_infile + '.abundfilt'
    n_outfile2 = n_infile + '2.abundfilt'

    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])

    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

    args = [n_counting_ht, n_infile]
    utils.runscript(script, args, n_in_dir)

    seqs = set([r.sequence for r in screed.open(n_infile)])
    assert os.path.exists(n_outfile), n_outfile

    args = [n_counting_ht, n_infile, '-o', n_outfile2]
    utils.runscript(script, args, in_dir)
    assert os.path.exists(n_outfile2), n_outfile2


def test_filter_abund_2():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = ['-C', '1', counting_ht, infile, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

# make sure that FASTQ records are retained.


def test_filter_abund_3_fq_retained():
    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = ['-C', '1', counting_ht, infile, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

    # check for 'quality' string.
    quals = set([r.quality for r in screed.open(outfile)])
    assert len(quals) == 2, quals
    assert '##################' in quals


# make sure that FASTQ names are properly parsed, both formats.


def test_filter_abund_4_fq_casava_18():
    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.paired2.fq'),
                    infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = [counting_ht, infile, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile, parse_description=False)])
    assert 'pair:foo 1::N' in seqs, seqs


def test_filter_abund_1_singlefile():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'filter-abund-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-t', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 98' in err, err

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_filter_abund_2_singlefile():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    tabfile = utils.get_temp_filename('test-savetable.ct')

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'filter-abund-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-t', '--savetable',
            tabfile, infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 98' in err, err

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_filter_abund_2_singlefile_fq_casava_18():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.paired2.fq'),
                    infile)

    script = 'filter-abund-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile, parse_description=False)])
    assert 'pair:foo 1::N' in seqs, seqs


def test_filter_abund_4_retain_low_abund():
    # test that the -V option does not trim sequences that are low abundance
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = ['-V', counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_filter_abund_5_trim_high_abund():
    # test that the -V option *does* trim sequences that are high abundance
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-3.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = ['-V', counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs

    # trimmed sequence @ error
    assert 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGC' in seqs


def test_filter_abund_6_trim_high_abund_Z():
    # test that -V/-Z settings interact properly -
    # trimming should not happen if -Z is set high enough.

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-3.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = ['-V', '-Z', '25', counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs

    # untrimmed seq.
    badseq = 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCgtgCCGCAGCTGTCGTCAGGG' \
             'GATTTCCGGGCGG'
    assert badseq in seqs       # should be there, untrimmed


def test_filter_abund_7_retain_Ns():
    # check that filter-abund retains sequences with Ns, and treats them as As.

    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    # copy test file over to test.fq & load into counting table
    shutil.copyfile(utils.get_test_data('test-filter-abund-Ns.fq'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = ['-C', '3', counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    # test for a sequence with an 'N' in it --
    names = set([r.name for r in screed.open(outfile, parse_description=0)])
    assert '895:1:37:17593:9954 1::FOO_withN' in names, names

    # check to see if that 'N' was properly changed to an 'A'
    seqs = set([r.sequence for r in screed.open(outfile)])
    assert 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAG' not in seqs, seqs

    # ...and that an 'N' remains in the output sequences
    found_N = False
    for s in seqs:
        if 'N' in s:
            found_N = True
    assert found_N, seqs


def test_filter_abund_single_8_retain_Ns():
    # check that filter-abund-single retains
    # sequences with Ns, and treats them as As.

    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    # copy test file over to test.fq & load into counting table
    shutil.copyfile(utils.get_test_data('test-filter-abund-Ns.fq'), infile)

    script = 'filter-abund-single.py'
    args = ['-k', '17', '-x', '1e7', '-N', '2', '-C', '3', infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    # test for a sequence with an 'N' in it --
    names = set([r.name for r in screed.open(outfile, parse_description=0)])
    assert '895:1:37:17593:9954 1::FOO_withN' in names, names

    # check to see if that 'N' was properly changed to an 'A'
    seqs = set([r.sequence for r in screed.open(outfile)])
    assert 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAG' not in seqs, seqs

    # ...and that an 'N' remains in the output sequences
    found_N = False
    for s in seqs:
        if 'N' in s:
            found_N = True
    assert found_N, seqs


def test_filter_stoptags():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    stopfile = utils.get_temp_filename('stoptags', in_dir)

    # first, copy test-abund-read-2.fa to 'test.fa' in the temp dir.
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    # now, create a file with some stop tags in it --
    K = 18
    kh = khmer.new_hashbits(K, 1, 1)
    kh.add_stop_tag('GTTGACGGGGCTCAGGGG')
    kh.save_stop_tags(stopfile)
    del kh

    # finally, run filter-stoptags.
    script = 'filter-stoptags.py'
    args = ['-k', str(K), stopfile, infile, infile]
    utils.runscript(script, args, in_dir)

    # verify that the basic output file exists
    outfile = infile + '.stopfilt'
    assert os.path.exists(outfile), outfile

    # it should contain only one unique sequence, because we've trimmed
    # off everything after the beginning of the only long sequence in there.
    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs, seqs


def test_filter_stoptags_fq():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    stopfile = utils.get_temp_filename('stoptags', in_dir)

    # first, copy test-abund-read-2.fa to 'test.fa' in the temp dir.
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile)

    # now, create a file with some stop tags in it --
    K = 18
    kh = khmer.new_hashbits(K, 1, 1)
    kh.add_stop_tag('GTTGACGGGGCTCAGGGG')
    kh.save_stop_tags(stopfile)
    del kh

    # finally, run filter-stoptags.
    script = 'filter-stoptags.py'
    args = ['-k', str(K), stopfile, infile, infile]
    utils.runscript(script, args, in_dir)

    # verify that the basic output file exists
    outfile = infile + '.stopfilt'
    assert os.path.exists(outfile), outfile

    # it should contain only one unique sequence, because we've trimmed
    # off everything after the beginning of the only long sequence in there.
    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs, seqs

    # make sure that record names are carried through unparsed
    names = [r.name for r in screed.open(outfile, parse_description=False)]
    names = set(names)
    assert 'seq 1::BAR' in names


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
    args = ['-C', CUTOFF, '-k', '17', '-p', infile]
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


def test_count_median():
    infile = utils.get_temp_filename('test.fa')
    outfile = infile + '.counts'

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = 'count-median.py'
    args = [counting_ht, infile, outfile]
    utils.runscript(script, args)

    assert os.path.exists(outfile), outfile

    data = [x.strip() for x in open(outfile)]
    data = set(data)
    assert len(data) == 2, data
    assert 'seq 1001 1001.0 0.0 18' in data
    assert '895:1:37:17593:9954/1 1 103.803741455 303.702941895 114' in data


def test_count_median_fq():
    infile = utils.get_temp_filename('test.fa')
    outfile = infile + '.counts'

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = 'count-median.py'
    args = [counting_ht, infile, outfile]
    utils.runscript(script, args)

    assert os.path.exists(outfile), outfile

    data = [x.strip() for x in open(outfile)]
    data = set(data)
    assert len(data) == 2, data
    assert 'seq 1001 1001.0 0.0 18' in data
    assert '895:1:37:17593:9954 1 103.803741455 303.702941895 114' in data


def test_count_median_fq_csv():
    infile = utils.get_temp_filename('test.fa')
    outfile = infile + '.counts'

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = 'count-median.py'
    args = ['--csv', counting_ht, infile, outfile]
    utils.runscript(script, args)

    assert os.path.exists(outfile), outfile

    data = [x.strip() for x in open(outfile)]
    data = set(data)
    assert len(data) == 4, data
    assert 'name,median,average,stddev,seqlen' in data
    assert 'seq,1001,1001.0,0.0,18' in data

    # verify that sequence names remain unparsed with '--csv'
    names = set([line.split(',')[0] for line in data])
    assert '895:1:37:17593:9954 1::FOO' in names, names


def test_load_graph():
    script = 'load-graph.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    ht = khmer.load_hashbits(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph worked :)
    subset = ht.do_subset_partition(0, 0)
    x = ht.subset_count_partitions(subset)
    assert x == (1, 0), x


def test_oxli_build_graph():
    script = 'oxli'
    args = ['build-graph', '-x', '1e7', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    ht = khmer.load_hashbits(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph worked :)
    subset = ht.do_subset_partition(0, 0)
    x = ht.subset_count_partitions(subset)
    assert x == (1, 0), x


def test_load_graph_no_tags():
    script = 'load-graph.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert not os.path.exists(tagset_file), tagset_file

    assert khmer.load_hashbits(ht_file)

    # can't think of a good way to make sure this worked, beyond just
    # loading the ht file...


def test_oxli_build_graph_no_tags():
    script = 'oxli'
    args = ['build-graph', '-x', '1e7', '-N', '2', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert not os.path.exists(tagset_file), tagset_file

    assert khmer.load_hashbits(ht_file)

    # can't think of a good way to make sure this worked, beyond just
    # loading the ht file...


def test_load_graph_fail():
    script = 'load-graph.py'
    args = ['-x', '1e3', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == 1, status
    assert "** ERROR: the graph structure is too small" in err


def test_oxli_build_graph_fail():
    script = 'oxli'
    args = ['build-graph', '-x', '1e3', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status == 1, status
    assert "** ERROR: the graph structure is too small" in err


def test_load_graph_write_fp():
    script = 'load-graph.py'
    args = ['-x', '1e5', '-N', '2', '-k', '20']  # use small HT

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    info_file = outfile + '.info'
    assert os.path.exists(info_file), info_file
    data = [x.strip() for x in open(info_file)]
    data = set(data)
    assert '3959 unique k-mers' in data, data
    assert 'false positive rate estimated to be 0.002' in data


def test_oxli_build_graph_write_fp():
    script = 'oxli'
    # use small HT
    args = ['build-graph', '-x', '1e5', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    info_file = outfile + '.info'
    assert os.path.exists(info_file), info_file
    data = [x.strip() for x in open(info_file)]
    data = set(data)
    assert '3959 unique k-mers' in data
    assert 'false positive rate estimated to be 0.002' in data


def test_load_graph_multithread():
    script = 'load-graph.py'

    outfile = utils.get_temp_filename('test')
    infile = utils.get_test_data('test-reads.fa')

    args = ['-N', '4', '-x', '1e7', '-T', '8', outfile, infile]

    (status, out, err) = utils.runscript(script, args)


def test_oxli_build_graph_multithread():
    script = 'oxli'

    outfile = utils.get_temp_filename('test')
    infile = utils.get_test_data('test-reads.fa')

    args = ['build-graph', '-N', '4', '-x', '1e7', '-T', '8', outfile, infile]

    (status, out, err) = utils.runscript(script, args)


def _make_graph(infilename, min_hashsize=1e7, n_hashes=2, ksize=20,
                do_partition=False,
                annotate_partitions=False,
                stop_big_traverse=False):
    script = 'load-graph.py'
    args = ['-x', str(min_hashsize), '-N', str(n_hashes), '-k', str(ksize)]

    outfile = utils.get_temp_filename('out')
    infile = infilename

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile + '.pt'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    if do_partition:
        script = 'partition-graph.py'
        args = [outfile]
        if stop_big_traverse:
            args.insert(0, '--no-big-traverse')
        utils.runscript(script, args)

        script = 'merge-partitions.py'
        args = [outfile, '-k', str(ksize)]
        utils.runscript(script, args)

        final_pmap_file = outfile + '.pmap.merged'
        assert os.path.exists(final_pmap_file)

        if annotate_partitions:
            script = 'annotate-partitions.py'
            args = ["-k", str(ksize), outfile, infilename]

            in_dir = os.path.dirname(outfile)
            utils.runscript(script, args, in_dir)

            baseinfile = os.path.basename(infilename)
            assert os.path.exists(os.path.join(in_dir, baseinfile + '.part'))

    return outfile


def _DEBUG_make_graph(infilename, min_hashsize=1e7, n_hashes=2, ksize=20,
                      do_partition=False,
                      annotate_partitions=False,
                      stop_big_traverse=False):
    script = 'load-graph.py'
    args = ['-x', str(min_hashsize), '-N', str(n_hashes), '-k', str(ksize)]

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data(infilename)

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile + '.ct'
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    if do_partition:
        print(">>>> DEBUG: Partitioning <<<")
        script = 'partition-graph.py'
        args = [outfile]
        if stop_big_traverse:
            args.insert(0, '--no-big-traverse')
        utils.runscript(script, args)

        print(">>>> DEBUG: Merging Partitions <<<")
        script = 'merge-partitions.py'
        args = [outfile, '-k', str(ksize)]
        utils.runscript(script, args)

        final_pmap_file = outfile + '.pmap.merged'
        assert os.path.exists(final_pmap_file)

        if annotate_partitions:
            print(">>>> DEBUG: Annotating Partitions <<<")
            script = 'annotate-partitions.py'
            args = ["-k", str(ksize), outfile, infilename]

            in_dir = os.path.dirname(outfile)
            utils.runscript(script, args, in_dir)

            baseinfile = os.path.basename(infilename)
            assert os.path.exists(os.path.join(in_dir, baseinfile + '.part'))

    return outfile


def test_partition_graph_1():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    script = 'partition-graph.py'
    args = [graphbase]

    utils.runscript(script, args)

    script = 'merge-partitions.py'
    args = [graphbase, '-k', str(20)]
    utils.runscript(script, args)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0), x          # should be exactly one partition.


def test_partition_graph_nojoin_k21():
    # test with K=21
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'), ksize=21)

    script = 'partition-graph.py'
    args = [graphbase]

    utils.runscript(script, args)

    script = 'merge-partitions.py'
    args = [graphbase, '-k', str(21)]
    utils.runscript(script, args)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (99, 0), x          # should be 99 partitions at K=21


def test_partition_graph_nojoin_stoptags():
    # test with stoptags
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    # add in some stop tags
    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.add_stop_tag('TTGCATACGTTGAGCCAGCG')
    stoptags_file = graphbase + '.stoptags'
    ht.save_stop_tags(stoptags_file)
    del ht

    # run script with stoptags option
    script = 'partition-graph.py'
    args = ['--stoptags', stoptags_file, graphbase]

    utils.runscript(script, args)

    script = 'merge-partitions.py'
    args = [graphbase, '-k', str(20)]
    utils.runscript(script, args)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (2, 0), x          # should be 2 partitions


def test_partition_graph_big_traverse():
    graphbase = _make_graph(utils.get_test_data('biglump-random-20-a.fa'),
                            do_partition=True, stop_big_traverse=False)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (1, 0), x          # should be exactly one partition.


def test_partition_graph_no_big_traverse():
    # do NOT exhaustively traverse
    graphbase = _make_graph(utils.get_test_data('biglump-random-20-a.fa'),
                            do_partition=True, stop_big_traverse=True)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_hashbits(graphbase + '.pt')
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x[0] == 4, x       # should be four partitions, broken at knot.


def test_partition_find_knots_execute():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    script = 'partition-graph.py'
    args = [graphbase]

    utils.runscript(script, args)

    script = 'find-knots.py'
    args = [graphbase]
    utils.runscript(script, args)

    stoptags_file = graphbase + '.stoptags'
    assert os.path.exists(stoptags_file)


def test_annotate_partitions():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(seqfile, do_partition=True)
    in_dir = os.path.dirname(graphbase)

    # get the final pmap file
    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    script = 'annotate-partitions.py'
    args = ["-k", "20", graphbase, seqfile]
    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    assert '2' in parts
    assert len(parts) == 1


def test_annotate_partitions_2():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(seqfile, do_partition=True,
                            ksize=21)
    in_dir = os.path.dirname(graphbase)

    # get the final pmap file
    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    script = 'annotate-partitions.py'
    args = ["-k", "21", graphbase, seqfile]
    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    print(parts)
    assert len(parts) == 99, len(parts)


def test_extract_partitions():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)


def test_extract_partitions_header_whitespace():
    seqfile = utils.get_test_data('test-overlap2.fa')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'test-overlap2.fa.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)

    dist = open(distfile).readline()
    assert dist.strip() == '1 11957 11957 11957'

    parts = [r.name.split('\t')[1]
             for r in screed.open(partfile, parse_description=False)]
    assert len(parts) == 13538, len(parts)
    parts = set(parts)
    assert len(parts) == 12601, len(parts)


def test_extract_partitions_fq():
    seqfile = utils.get_test_data('random-20-a.fq')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fq.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fq')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    screed_iter = screed.open(partfile, parse_description=False)
    names = [r.name.split('\t')[0] for r in screed_iter]
    assert '35 1::FOO' in names
    assert '46 1::FIZ' in names

    screed_iter = screed.open(partfile, parse_description=False)
    parts = [r.name.split('\t')[1] for r in screed_iter]

    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)

    quals = set([r.quality for r in screed.open(partfile)])
    quals = list(quals)
    assert quals[0], quals


def test_extract_partitions_output_unassigned():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['-U', 'extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    unassigned_file = os.path.join(in_dir, 'extracted.unassigned.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)
    assert os.path.exists(unassigned_file)

    dist = open(distfile).readline()
    assert dist.strip() == '99 1 1 99'

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    assert len(parts) == 99, len(parts)
    parts = set(parts)
    assert len(parts) == 1, len(parts)


def test_extract_partitions_no_output_groups():
    seqfile = utils.get_test_data('random-20-a.fq')
    graphbase = _make_graph(
        seqfile, do_partition=True, annotate_partitions=True)
    in_dir = os.path.dirname(graphbase)

    # get the final part file
    partfile = os.path.join(in_dir, 'random-20-a.fq.part')

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['-n', 'extracted', partfile]

    # We expect a sys.exit -> we need the test to be tolerant
    _, out, err = utils.runscript(script, args, in_dir, fail_ok=True)
    assert "NOT outputting groups! Beware!" in err
    # Group files are created after output_groups is
    # checked. They should not exist in this scenario
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    assert not os.path.exists(groupfile)


def test_extract_partitions_pid_0():
    basefile = utils.get_test_data('random-20-a.fa.part')
    partfile = utils.get_temp_filename('random-20-a.fa.part')
    shutil.copyfile(basefile, partfile)

    in_dir = os.path.dirname(partfile)
    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['-U', 'extracted', partfile]

    utils.runscript(script, args, in_dir)

    distfile = os.path.join(in_dir, 'extracted.dist')
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')
    unassigned_file = os.path.join(in_dir, 'extracted.unassigned.fa')
    assert os.path.exists(distfile)
    assert os.path.exists(groupfile)
    assert os.path.exists(unassigned_file)

    # Assert unassigned file not empty
    unassigned_content = open(unassigned_file).readline()
    assert unassigned_content.strip().split('\t')[0] != ''


def test_extract_partitions_multi_groups():
    basefile = utils.get_test_data('random-20-a.fa.part')
    partfile = utils.get_temp_filename('random-20-a.fa.part')
    shutil.copyfile(basefile, partfile)

    in_dir = os.path.dirname(partfile)

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['-m', '1', '-X', '1', 'extracted', partfile]

    utils.runscript(script, args, in_dir)

    # Multiple group files are created after should be created
    groupfile1 = os.path.join(in_dir, 'extracted.group0000.fa')
    groupfile2 = os.path.join(in_dir, 'extracted.group0001.fa')
    groupfile3 = os.path.join(in_dir, 'extracted.group0002.fa')
    assert os.path.exists(groupfile1)
    assert os.path.exists(groupfile2)
    assert os.path.exists(groupfile3)


def test_extract_partitions_no_groups():
    empty_file = utils.get_temp_filename('empty-file')
    basefile = utils.get_test_data('empty-file')

    shutil.copyfile(basefile, empty_file)
    in_dir = os.path.dirname(empty_file)

    # ok, now run extract-partitions.
    script = 'extract-partitions.py'
    args = ['extracted', empty_file]

    _, _, err = utils.runscript(script, args, in_dir, fail_ok=True)
    assert "ERROR: Input file", "is empty; Exiting." in err
    # No group files should be created
    groupfile = os.path.join(in_dir, 'extracted.group0000.fa')

    assert not os.path.exists(groupfile)


def test_abundance_dist():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    htfile = _make_counting(infile, K=17)

    script = 'abundance-dist.py'
    args = ['-z', htfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert line == '1 96 96 0.98', line
        line = fp.readline().strip()
        assert line == '1001 2 98 1.0', line

    os.remove(outfile)
    args = ['-z', '--csv', htfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert (line == 'abundance,count,cumulative,cumulative_fraction'), line
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_nobigcount():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    htfile = _make_counting(infile, K=17, BIGCOUNT=False)

    script = 'abundance-dist.py'
    args = ['-z', htfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert line == '1 96 96 0.98', line
        line = fp.readline().strip()
        assert line == '255 2 98 1.0', line


def test_abundance_dist_single():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '-t', infile,
            outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 98' in err, err

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert line == '1 96 96 0.98', line
        line = fp.readline().strip()
        assert line == '1001 2 98 1.0', line


def test_abundance_dist_threaded():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '-t', '--threads', '18',
            infile, outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 98' in err, err

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert line == '1 96 96 0.98', line
        line = fp.readline().strip()
        assert line == '1001 2 98 1.0', line


def test_abundance_dist_single_csv():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '--csv', infile,
            outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert (line == 'abundance,count,cumulative,cumulative_fraction'), line
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_single_nobigcount():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '-b', infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert line == '1 96 96 0.98', line
        line = fp.readline().strip()
        assert line == '255 2 98 1.0', line


def test_abundance_dist_single_nosquash():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '-t', infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert line == '1 96 96 0.98', line
        line = fp.readline().strip()
        assert line == '1001 2 98 1.0', line


def test_abundance_dist_single_savetable():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    tabfile = utils.get_temp_filename('test-savetable.ct')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '-t', '--savetable',
            tabfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()
        assert line == '1 96 96 0.98', line
        line = fp.readline().strip()
        assert line == '1001 2 98 1.0', line


def test_do_partition():
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = 'do-partition.py'
    args = ["-k", "20", graphbase, seqfile]

    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)
    assert '2' in parts
    assert len(parts) == 1


def test_do_partition_2():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fa')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = 'do-partition.py'
    args = ["-k", "21", graphbase, seqfile]

    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fa.part')

    parts = [r.name.split('\t')[1] for r in screed.open(partfile)]
    parts = set(parts)

    assert len(parts) == 99, len(parts)


def test_do_partition_2_fq():
    # test with K=21 (no joining of sequences)
    seqfile = utils.get_test_data('random-20-a.fq')
    graphbase = utils.get_temp_filename('out')
    in_dir = os.path.dirname(graphbase)

    script = 'do-partition.py'
    args = ["-k", "21", graphbase, seqfile]

    utils.runscript(script, args, in_dir)

    partfile = os.path.join(in_dir, 'random-20-a.fq.part')

    screed_iter = screed.open(partfile, parse_description=False)
    names = [r.name.split('\t')[0] for r in screed_iter]
    assert '35 1::FOO' in names
    assert '46 1::FIZ' in names


def test_interleave_reads_1_fq():
    # test input files
    infile1 = utils.get_test_data('paired.fq.1')
    infile2 = utils.get_test_data('paired.fq.2')

    # correct output
    ex_outfile = utils.get_test_data('paired.fq')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, infile2, '-o', outfile]

    utils.runscript(script, args)

    r = open(ex_outfile).read()
    q = open(outfile).read()

    assert r == q, (r, q)


def test_interleave_reads_broken_fq():
    # test input files
    infile1 = utils.get_test_data('paired-broken.fq.1')
    infile2 = utils.get_test_data('paired-broken.fq.2')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, infile2, '-o', outfile]

    status, out, err = utils.runscript(script, args, fail_ok=True)
    assert status == 1
    assert 'ERROR: Input files contain different number of records.' in err


def test_interleave_reads_broken_fq_2():
    # test input files
    infile1 = utils.get_test_data('paired-broken2.fq.1')
    infile2 = utils.get_test_data('paired-broken2.fq.2')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, infile2, '-o', outfile]

    status, out, err = utils.runscript(script, args, fail_ok=True)
    assert status == 1
    assert "ERROR: This doesn't look like paired data!" in err


def test_interleave_reads_broken_fq_3():
    # test input files
    infile1 = utils.get_test_data('paired-broken3.fq.1')
    infile2 = utils.get_test_data('paired-broken3.fq.2')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, infile2, '-o', outfile]

    status, out, err = utils.runscript(script, args, fail_ok=True)
    assert status == 1
    assert "ERROR: This doesn't look like paired data!" in err


def test_interleave_reads_broken_fq_4():
    # test input files
    infile1 = utils.get_test_data('paired-mixed-broken.fq')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, '-o', outfile]

    status, out, err = utils.runscript(script, args, fail_ok=True)
    assert status == 1
    assert "ERROR: given only one filename, that doesn't contain _R1_" in err


def test_interleave_reads_2_fa():
    # test input files
    infile1 = utils.get_test_data('paired.fa.1')
    infile2 = utils.get_test_data('paired.fa.2')

    # correct output
    ex_outfile = utils.get_test_data('paired.fa')

    # actual output file
    outfile = utils.get_temp_filename('out.fa')

    script = 'interleave-reads.py'
    args = [infile1, infile2, '-o', outfile]

    utils.runscript(script, args)

    n = 0
    for r, q in zip(screed.open(ex_outfile), screed.open(outfile)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0


def test_make_initial_stoptags():
    # gen input files using load-graph.py -t
    # should keep test_data directory size down
    # or something like that
    # this assumes (obv.) load-graph works properly
    bzinfile = utils.get_temp_filename('test-reads.fq.bz2')
    shutil.copyfile(utils.get_test_data('test-reads.fq.bz2'), bzinfile)
    in_dir = os.path.dirname(bzinfile)

    genscript = 'load-graph.py'
    genscriptargs = ['test-reads', 'test-reads.fq.bz2']
    utils.runscript(genscript, genscriptargs, in_dir)

    # test input file gen'd by load-graphs
    infile = utils.get_temp_filename('test-reads.pt')
    infile2 = utils.get_temp_filename('test-reads.tagset', in_dir)

    # get file to compare against
    ex_outfile = utils.get_test_data('test-reads.stoptags')

    # actual output file
    outfile1 = utils.get_temp_filename('test-reads.stoptags', in_dir)

    script = 'make-initial-stoptags.py'
    # make-initial-stoptags has weird file argument syntax
    # read the code before modifying
    args = ['test-reads']

    utils.runscript(script, args, in_dir)
    assert os.path.exists(outfile1), outfile1


def test_extract_paired_reads_1_fa():
    # test input file
    infile = utils.get_test_data('paired-mixed.fa')

    ex_outfile1 = utils.get_test_data('paired-mixed.fa.pe')
    ex_outfile2 = utils.get_test_data('paired-mixed.fa.se')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired-mixed.fa.pe')
    in_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired-mixed.fa.se', in_dir)

    script = 'extract-paired-reads.py'
    args = [infile]

    utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0


def test_extract_paired_reads_2_fq():
    # test input file
    infile = utils.get_test_data('paired-mixed.fq')

    ex_outfile1 = utils.get_test_data('paired-mixed.fq.pe')
    ex_outfile2 = utils.get_test_data('paired-mixed.fq.se')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired-mixed.fq.pe')
    in_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired-mixed.fq.se', in_dir)

    script = 'extract-paired-reads.py'
    args = [infile]

    utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1, parse_description=False),
                    screed.open(outfile1, parse_description=False)):
        n += 1
        assert r.name == q.name, (r.name, q.name, n)
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2, parse_description=False),
                    screed.open(outfile2, parse_description=False)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0


def execute_split_paired_streaming(ifilename):
    fifo = utils.get_temp_filename('fifo')
    in_dir = os.path.dirname(fifo)
    outfile1 = utils.get_temp_filename('paired-1.fa')
    outfile2 = utils.get_temp_filename('paired-2.fa')
    script = 'split-paired-reads.py'
    args = [fifo, '-1', outfile1, '-2', outfile2]

    # make a fifo to simulate streaming
    os.mkfifo(fifo)

    thread = threading.Thread(target=utils.runscript,
                              args=(script, args, in_dir))
    thread.start()
    ifile = open(ifilename, 'r')
    fifofile = open(fifo, 'w')
    chunk = ifile.read(4)
    while len(chunk) > 0:
        fifofile.write(chunk)
        chunk = ifile.read(4)
    fifofile.close()
    thread.join()
    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2


def test_split_paired_streaming():
    o = execute_split_paired_streaming(utils.get_test_data('paired.fa'))


def test_split_paired_reads_1_fa():
    # test input file
    infile = utils.get_test_data('paired.fa')

    ex_outfile1 = utils.get_test_data('paired.fa.1')
    ex_outfile2 = utils.get_test_data('paired.fa.2')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired.fa.1')
    in_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired.fa.2', in_dir)

    script = 'split-paired-reads.py'
    args = [infile]

    utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0


def test_split_paired_reads_2_fq():
    # test input file
    infile = utils.get_test_data('paired.fq')

    ex_outfile1 = utils.get_test_data('paired.fq.1')
    ex_outfile2 = utils.get_test_data('paired.fq.2')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired.fq.1')
    in_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired.fq.2', in_dir)

    script = 'split-paired-reads.py'
    args = [infile]

    utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0


def test_split_paired_reads_2_mixed_fq_require_pair():
    # test input file
    infile = utils.get_temp_filename('test.fq')
    shutil.copyfile(utils.get_test_data('paired-mixed.fq'), infile)
    in_dir = os.path.dirname(infile)

    script = 'split-paired-reads.py'
    args = ['-p', infile]

    status, out, err = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status == 1
    assert "is not part of a pair" in err


def test_split_paired_reads_2_mixed_fq():
    # test input file
    infile = utils.get_temp_filename('test.fq')
    shutil.copyfile(utils.get_test_data('paired-mixed-2.fq'), infile)
    in_dir = os.path.dirname(infile)

    script = 'split-paired-reads.py'
    args = [infile]

    status, out, err = utils.runscript(script, args, in_dir)
    assert status == 0
    assert "split 11 sequences (7 left, 4 right)" in err, err


def test_split_paired_reads_2_mixed_fq_broken_pairing_format():
    # test input file
    infile = utils.get_temp_filename('test.fq')
    shutil.copyfile(utils.get_test_data('paired-mixed-broken.fq'), infile)
    in_dir = os.path.dirname(infile)

    script = 'split-paired-reads.py'
    args = [infile]

    status, out, err = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status == 1
    assert "Unrecognized format" in err


def test_split_paired_reads_3_output_dir():
    # test input file
    infile = utils.get_test_data('paired.fq')

    ex_outfile1 = utils.get_test_data('paired.fq.1')
    ex_outfile2 = utils.get_test_data('paired.fq.2')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired.fq.1')
    output_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired.fq.2', output_dir)

    script = 'split-paired-reads.py'
    args = ['--output-dir', output_dir, infile]

    utils.runscript(script, args)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0


def test_split_paired_reads_3_output_files():
    # test input file
    infile = utils.get_test_data('paired.fq')

    ex_outfile1 = utils.get_test_data('paired.fq.1')
    ex_outfile2 = utils.get_test_data('paired.fq.2')

    # actual output files...
    outfile1 = utils.get_temp_filename('xxx')
    output_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('yyy', output_dir)

    script = 'split-paired-reads.py'
    args = ['-1', outfile1, '-2', outfile2, infile]

    utils.runscript(script, args)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0


def test_split_paired_reads_3_output_files_left():
    # test input file
    infile = utils.get_test_data('paired.fq')

    ex_outfile1 = utils.get_test_data('paired.fq.1')
    ex_outfile2 = utils.get_test_data('paired.fq.2')

    # actual output files...
    outfile1 = utils.get_temp_filename('xxx')
    output_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('paired.fq.2', output_dir)

    script = 'split-paired-reads.py'
    args = ['-o', output_dir, '-1', outfile1, infile]

    utils.runscript(script, args)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0


def test_split_paired_reads_3_output_files_right():
    # test input file
    infile = utils.get_test_data('paired.fq')

    ex_outfile1 = utils.get_test_data('paired.fq.1')
    ex_outfile2 = utils.get_test_data('paired.fq.2')

    # actual output files...
    outfile1 = utils.get_temp_filename('paired.fq.1')
    output_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('yyy', output_dir)

    script = 'split-paired-reads.py'
    args = ['-2', outfile2, '-o', output_dir, infile]

    utils.runscript(script, args)

    assert os.path.exists(outfile1), outfile1
    assert os.path.exists(outfile2), outfile2

    n = 0
    for r, q in zip(screed.open(ex_outfile1), screed.open(outfile1)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2), screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0


def test_sample_reads_randomly():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-reads.fa'), infile)

    script = 'sample-reads-randomly.py'
    # fix random number seed for reproducibility
    args = ['-N', '10', '-M', '12000', '-R', '1']
    args.append(infile)
    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])
    print(list(sorted(seqs)))

    if sys.version_info.major == 2:
        answer = {'850:2:1:1859:11742/1', '850:2:1:1859:11742/2',
                  '850:2:1:2131:17360/1', '850:2:1:2131:17360/2',
                  '850:2:1:2416:7565/1', '850:2:1:2416:7565/2',
                  '850:2:1:2490:13491/1', '850:2:1:2490:13491/2',
                  '850:2:1:2962:3999/1', '850:2:1:2962:3999/2',
                  '850:2:1:3096:20321/1', '850:2:1:3096:20321/2',
                  '850:2:1:3164:6414/1', '850:2:1:3164:6414/2',
                  '850:2:1:3206:13876/1', '850:2:1:3206:13876/2',
                  '850:2:1:3631:20919/1', '850:2:1:3631:20919/2',
                  '850:2:1:3655:15581/1', '850:2:1:3655:15581/2'}
    else:
        answer = {'850:2:1:1257:3404/1', '850:2:1:1257:3404/2',
                  '850:2:1:1362:19357/1', '850:2:1:1362:19357/2',
                  '850:2:1:1396:5659/1', '850:2:1:1396:5659/2',
                  '850:2:1:2063:11124/1', '850:2:1:2063:11124/2',
                  '850:2:1:2121:12070/1', '850:2:1:2121:12070/2',
                  '850:2:1:2528:15779/1', '850:2:1:2528:15779/2',
                  '850:2:1:2581:12886/1', '850:2:1:2581:12886/2',
                  '850:2:1:2864:8505/1', '850:2:1:2864:8505/2',
                  '850:2:1:3000:2015/1', '850:2:1:3000:2015/2',
                  '850:2:1:3302:5025/1', '850:2:1:3302:5025/2'}

    assert seqs == answer


def test_sample_reads_randomly_force_single():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-reads.fa'), infile)

    script = 'sample-reads-randomly.py'
    # fix random number seed for reproducibility
    args = ['-N', '10', '-M', '12000', '-R', '1', '--force_single']
    args.append(infile)
    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile)])
    print(list(sorted(seqs)))

    if sys.version_info.major == 2:
        answer = {'850:2:1:2399:20086/2',
                  '850:2:1:2273:13309/1',
                  '850:2:1:2065:16816/1',
                  '850:2:1:1984:7162/2',
                  '850:2:1:2691:14602/1',
                  '850:2:1:1762:5439/1',
                  '850:2:1:2503:4494/2',
                  '850:2:1:2263:11143/2',
                  '850:2:1:1792:15774/2',
                  '850:2:1:2084:17145/1'}
    else:
        answer = {'850:2:1:1199:4197/1',
                  '850:2:1:1251:16575/2',
                  '850:2:1:1267:6790/2',
                  '850:2:1:1601:4443/1',
                  '850:2:1:1625:19325/1',
                  '850:2:1:1832:14607/2',
                  '850:2:1:1946:20852/2',
                  '850:2:1:2401:4896/2',
                  '850:2:1:2562:1308/1',
                  '850:2:1:3123:15968/2'}

    assert seqs == answer


def test_sample_reads_randomly_fq():
    infile = utils.get_temp_filename('test.fq.gz')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-reads.fq.gz'), infile)

    script = 'sample-reads-randomly.py'
    # fix random number seed for reproducibility
    args = ['-N', '10', '-M', '12000', '-R', '1']
    args.append(infile)
    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset'
    assert os.path.exists(outfile), outfile

    if sys.version_info.major == 2:
        answer = {'850:2:1:2399:20086/2',
                  '850:2:1:1762:5439 1::FOO',
                  '850:2:1:2065:16816/1',
                  '850:2:1:2263:11143/2',
                  '850:2:1:1792:15774/2',
                  '850:2:1:2691:14602/1',
                  '850:2:1:2503:4494 1::FOO',
                  '850:2:1:2084:17145/1',
                  '850:2:1:1984:7162 1::FOO',
                  '850:2:1:2273:13309 1::FOO'}
    else:
        answer = {'850:2:1:1199:4197 1::FOO',
                  '850:2:1:1251:16575/2',
                  '850:2:1:1267:6790/2',
                  '850:2:1:1601:4443 1::FOO',
                  '850:2:1:1625:1932 1::FOO1',
                  '850:2:1:1832:14607 1::FOO',
                  '850:2:1:1946:20852 1::FOO',
                  '850:2:1:2401:4896/2',
                  '850:2:1:2562:1308/1',
                  '850:2:1:3123:15968/2'}

    seqs = set([r.name for r in screed.open(outfile,
                                            parse_description=False)])
    print(list(sorted(seqs)))
    assert seqs == answer


def test_fastq_to_fasta():

    script = 'fastq-to-fasta.py'
    clean_infile = utils.get_temp_filename('test-clean.fq')
    n_infile = utils.get_temp_filename('test-n.fq')

    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), clean_infile)
    shutil.copyfile(utils.get_test_data('test-fastq-n-reads.fq'), n_infile)

    clean_outfile = clean_infile + '.keep.fa'
    n_outfile = n_infile + '.keep.fa'

    in_dir = os.path.dirname(clean_infile)
    in_dir_n = os.path.dirname(n_infile)

    args = [clean_infile, '-n', '-o', clean_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 2, len(out.splitlines())
    assert "No lines dropped" in err

    names = [r.name for r in screed.open(clean_outfile,
                                         parse_description=False)]
    assert '895:1:1:1246:14654 1:N:0:NNNNN' in names, names

    args = [n_infile, '-n', '-o', n_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) == 2
    assert "No lines dropped" in err

    args = [clean_infile, '-o', clean_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 2
    assert "0 lines dropped" in err

    args = [n_infile, '-o', n_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) == 2, out
    assert "4 lines dropped" in err, err

    args = [clean_infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) > 2
    assert "0 lines dropped" in err

    args = [n_infile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) > 2
    assert "4 lines dropped" in err


def test_extract_long_sequences_fa():

    script = 'extract-long-sequences.py'
    fa_infile = utils.get_temp_filename('test.fa')

    shutil.copyfile(utils.get_test_data('paired-mixed.fa'), fa_infile)

    fa_outfile = fa_infile + '.keep.fa'

    in_dir_fa = os.path.dirname(fa_infile)

    args = [fa_infile, '-l', '10', '-o', fa_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_fa)

    countlines = sum(1 for line in open(fa_outfile))
    assert countlines == 22, countlines

    names = [r.name for r in screed.open(fa_outfile, parse_description=False)]
    assert "895:1:37:17593:9954/1" in names
    assert "895:1:37:17593:9954/2" in names


def test_extract_long_sequences_fq():

    script = 'extract-long-sequences.py'
    fq_infile = utils.get_temp_filename('test.fq')

    shutil.copyfile(utils.get_test_data('paired-mixed.fq'), fq_infile)

    fq_outfile = fq_infile + '.keep.fq'

    in_dir_fq = os.path.dirname(fq_infile)

    args = [fq_infile, '-l', '10', '-o', fq_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_fq)

    countlines = sum(1 for line in open(fq_outfile))
    assert countlines == 44, countlines

    names = [r.name for r in screed.open(fq_outfile, parse_description=False)]
    assert "895:1:37:17593:9954 1::foo" in names
    assert "895:1:37:17593:9954 2::foo" in names


def test_sample_reads_randomly_S():
    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), infile)

    script = 'sample-reads-randomly.py'

    # fix random number seed for reproducibility
    args = ['-N', '10', '-R', '1', '-S', '3']

    badargs = list(args)
    badargs.extend(['-o', 'test', infile, infile])
    (status, out, err) = utils.runscript(script, badargs, in_dir, fail_ok=True)
    assert status == 1, (status, out, err)
    assert "Error: cannot specify -o with more than one sample" in err

    args.append(infile)

    utils.runscript(script, args, in_dir)

    outfile = infile + '.subset.0'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile, parse_description=True)])
    print(list(sorted(seqs)))

    print(seqs)
    if sys.version_info.major == 2:
        answer = {'895:1:1:1303:14389', '895:1:1:1347:3237',
                  '895:1:1:1295:6189', '895:1:1:1308:20421',
                  '895:1:1:1320:11648', '895:1:1:1352:5369',
                  '895:1:1:1318:10532', '895:1:1:1363:11839',
                  '895:1:1:1355:13535', '895:1:1:1349:15165'}
    else:
        answer = {'895:1:1:1290:11501', '895:1:1:1303:14389',
                  '895:1:1:1307:4308', '895:1:1:1308:2539',
                  '895:1:1:1331:1766', '895:1:1:1333:2512',
                  '895:1:1:1347:3237', '895:1:1:1363:11839',
                  '895:1:1:1378:18986', '895:1:1:1383:3089'}

    assert seqs == answer

    outfile = infile + '.subset.1'
    assert os.path.exists(outfile), outfile

    seqs = set([r.name for r in screed.open(outfile, parse_description=True)])
    print(list(sorted(seqs)))

    if sys.version_info.major == 2:
        answer = set(['895:1:1:1303:14389', '895:1:1:1373:4848',
                      '895:1:1:1357:19736', '895:1:1:1347:3237',
                      '895:1:1:1338:7557', '895:1:1:1388:11093',
                      '895:1:1:1296:1784', '895:1:1:1290:11501',
                      '895:1:1:1355:13535', '895:1:1:1303:6251'])
    else:
        answer = {'895:1:1:1255:18861', '895:1:1:1276:16426',
                  '895:1:1:1303:6251', '895:1:1:1308:20421',
                  '895:1:1:1314:10430', '895:1:1:1351:14718',
                  '895:1:1:1355:13535', '895:1:1:1358:4953',
                  '895:1:1:1362:3983', '895:1:1:1363:9988'}
    assert seqs == answer

    seqs = set([r.name for r in screed.open(outfile, parse_description=True)])
    print(list(sorted(seqs)))

    if sys.version_info.major == 2:
        answer = {'895:1:1:1303:14389', '895:1:1:1373:4848',
                  '895:1:1:1357:19736', '895:1:1:1347:3237',
                  '895:1:1:1338:7557', '895:1:1:1388:11093',
                  '895:1:1:1296:1784', '895:1:1:1290:11501',
                  '895:1:1:1355:13535', '895:1:1:1303:6251'}

    else:
        answer = {'895:1:1:1362:3983', '895:1:1:1363:9988',
                  '895:1:1:1314:10430', '895:1:1:1255:18861',
                  '895:1:1:1308:20421', '895:1:1:1358:4953',
                  '895:1:1:1351:14718', '895:1:1:1303:6251',
                  '895:1:1:1276:16426', '895:1:1:1355:13535'}

    assert seqs == answer


def test_count_overlap_invalid_datafile():
    seqfile1 = utils.get_temp_filename('test-overlap1.fa')
    in_dir = os.path.dirname(seqfile1)
    shutil.copy(utils.get_test_data('test-overlap1.fa'), seqfile1)
    htfile = _make_graph(seqfile1, ksize=20)
    outfile = utils.get_temp_filename('overlap.out', in_dir)
    script = 'count-overlap.py'
    args = ['--ksize', '20', '--n_tables', '2', '--min-tablesize', '10000000',
            htfile + '.pt', htfile + '.pt', outfile]
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    if sys.version_info.major == 2:
        assert "IOError" in err
    else:
        assert "OSError" in err


def test_count_overlap():
    seqfile1 = utils.get_temp_filename('test-overlap1.fa')
    in_dir = os.path.dirname(seqfile1)
    seqfile2 = utils.get_temp_filename('test-overlap2.fa', in_dir)
    outfile = utils.get_temp_filename('overlap.out', in_dir)
    curvefile = utils.get_temp_filename('overlap.out.curve', in_dir)
    shutil.copy(utils.get_test_data('test-overlap1.fa'), seqfile1)
    shutil.copy(utils.get_test_data('test-overlap2.fa'), seqfile2)
    htfile = _make_graph(seqfile1, ksize=20)
    script = 'count-overlap.py'
    args = ['--ksize', '20', '--n_tables', '2', '--min-tablesize', '10000000',
            htfile + '.pt', seqfile2, outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert status == 0
    assert os.path.exists(outfile), outfile
    data = [x.strip() for x in open(outfile)]
    data = set(data)
    assert '# of unique k-mers in dataset2: 759047' in data
    assert '# of overlap unique k-mers: 245621' in data
    assert os.path.exists(curvefile), curvefile
    data = [x.strip() for x in open(curvefile)]
    data = set(data)
    assert '178633 1155' in data
    assert '496285 2970' in data
    assert '752053 238627' in data


def test_count_overlap_csv():
    seqfile1 = utils.get_temp_filename('test-overlap1.fa')
    in_dir = os.path.dirname(seqfile1)
    seqfile2 = utils.get_temp_filename('test-overlap2.fa', in_dir)
    outfile = utils.get_temp_filename('overlap.out', in_dir)
    curvefile = utils.get_temp_filename('overlap.out.curve', in_dir)
    shutil.copy(utils.get_test_data('test-overlap1.fa'), seqfile1)
    shutil.copy(utils.get_test_data('test-overlap2.fa'), seqfile2)
    htfile = _make_graph(seqfile1, ksize=20)
    script = 'count-overlap.py'
    args = ['--ksize', '20', '--n_tables', '2', '--min-tablesize',
            '10000000', '--csv', htfile + '.pt', seqfile2, outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert status == 0
    assert os.path.exists(outfile), outfile
    data = [x.strip() for x in open(outfile)]
    data = set(data)
    assert '# of unique k-mers in dataset2: 759047' in data
    assert '# of overlap unique k-mers: 245621' in data
    assert os.path.exists(curvefile), curvefile
    data = [x.strip() for x in open(curvefile)]
    data = set(data)
    assert '178633,1155' in data
    assert '496285,2970' in data
    assert '752053,238627' in data


def execute_streaming_diginorm(ifilename):
    '''Helper function for the matrix of streaming tests for read_parser
    using diginorm, i.e. uncompressed fasta, gzip fasta, bz2 fasta,
    uncompressed fastq, etc.
    This is not directly executed but is run by the tests themselves
    '''
    # Get temp filenames, etc.
    fifo = utils.get_temp_filename('fifo')
    in_dir = os.path.dirname(fifo)
    script = 'normalize-by-median.py'
    args = ['-C', '1', '-k', '17', '-o', 'outfile', fifo]

    # make a fifo to simulate streaming
    os.mkfifo(fifo)

    # FIFOs MUST BE OPENED FOR READING BEFORE THEY ARE WRITTEN TO
    # If this isn't done, they will BLOCK and things will hang.
    thread = threading.Thread(target=utils.runscript,
                              args=(script, args, in_dir))
    thread.start()
    ifile = io.open(ifilename, 'rb')
    fifofile = io.open(fifo, 'wb')
    # read binary to handle compressed files
    chunk = ifile.read(8192)
    while len(chunk) > 0:
        fifofile.write(chunk)
        chunk = ifile.read(8192)

    fifofile.close()

    thread.join()

    return in_dir + '/outfile'


def execute_load_graph_streaming(filename):
    '''Helper function for the matrix of streaming tests using screed via
    filter-abund-single, i.e. uncompressed fasta, gzip fasta, bz2 fasta,
    uncompressed fastq, etc.
    This is not directly executed but is run by the tests themselves
    '''

    script = 'load-graph.py'
    args = '-x 1e7 -N 2 -k 20 out -'

    infile = utils.get_temp_filename('temp')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data(filename), infile)
    (status, out, err) = utils.runscriptredirect(script, args, infile, in_dir)

    if status != 0:
        for line in out:
            print(out)
        for line in err:
            print(err)
        assert status == 0, status
    err.seek(0)
    err = err.read()
    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = os.path.join(in_dir, 'out.pt')
    assert os.path.exists(ht_file), ht_file

    tagset_file = os.path.join(in_dir, 'out.tagset')
    assert os.path.exists(tagset_file), tagset_file

    ht = khmer.load_hashbits(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph worked :)
    subset = ht.do_subset_partition(0, 0)
    x = ht.subset_count_partitions(subset)
    assert x == (1, 0), x


def test_screed_streaming_ufa():
    # uncompressed fa
    o = execute_streaming_diginorm(utils.get_test_data('test-abund-read-2.fa'))

    pathstat = os.stat(o)
    seqs = [r.sequence for r in screed.open(o)]
    assert len(seqs) == 1, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG')


def test_screed_streaming_ufq():
    # uncompressed fq
    o = execute_streaming_diginorm(utils.get_test_data('test-fastq-reads.fq'))

    seqs = [r.sequence for r in screed.open(o)]
    assert seqs[0].startswith('CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT')


def test_screed_streaming_bzipfq():
    # bzip compressed fq
    o = execute_streaming_diginorm(utils.get_test_data('100-reads.fq.bz2'))
    seqs = [r.sequence for r in screed.open(o)]
    assert len(seqs) == 100, seqs
    assert seqs[0].startswith('CAGGCGCCCACCACCGTGCCCTCCAACCTGATGGT'), seqs


def test_screed_streaming_bzipfa():
    # bzip compressed fa
    o = execute_streaming_diginorm(
        utils.get_test_data('test-abund-read-2.fa.bz2'))

    seqs = [r.sequence for r in screed.open(o)]
    assert len(seqs) == 1, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG')


@attr('known_failing')
def test_screed_streaming_gzipfq():
    # gzip compressed fq
    o = execute_streaming_diginorm(utils.get_test_data('100-reads.fq.gz'))
    assert os.path.exists(o)
    seqs = [r.sequence for r in screed.open(o)]
    assert seqs[0].startswith('CAGGCGCCCACCACCGTGCCCTCCAACCTG')


@attr('known_failing')
def test_screed_streaming_gzipfa():
    o = execute_streaming_diginorm(
        utils.get_test_data('test-abund-read-2.fa.gz'))
    assert os.path.exists(o)
    seqs = [r.sequence for r in screed.open(o)]
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGG')


def test_read_parser_streaming_ufa():
    # uncompressed FASTA
    execute_load_graph_streaming(utils.get_test_data('random-20-a.fa'))


def test_read_parser_streaming_ufq():
    # uncompressed FASTQ
    execute_load_graph_streaming(utils.get_test_data('random-20-a.fq'))


@attr('known_failing')
def test_read_parser_streaming_bzfq():
    # bzip compressed FASTQ
    execute_load_graph_streaming(utils.get_test_data('random-20-a.fq.bz2'))


def test_read_parser_streaming_gzfq():
    # gzip compressed FASTQ
    execute_load_graph_streaming(utils.get_test_data('random-20-a.fq.gz'))


@attr('known_failing')
def test_read_parser_streaming_bzfa():
    # bzip compressed FASTA
    execute_load_graph_streaming(utils.get_test_data('random-20-a.fa.bz2'))


def test_read_parser_streaming_gzfa():
    # gzip compressed FASTA
    execute_load_graph_streaming(utils.get_test_data('random-20-a.fa.gz'))


def test_readstats():
    readstats_output = ("358 bp / 5 seqs; 71.6 average length",
                        "916 bp / 11 seqs; 83.3 average length")

    args = [utils.get_test_data("test-sweep-reads.fq"),
            utils.get_test_data("paired-mixed.fq")]
    status, out, err = utils.runscript('readstats.py', args)
    assert status == 0

    for k in readstats_output:
        assert k in out, (k, out)


def test_readstats_csv():
    readstats_output = ("358,5,71.6," +
                        utils.get_test_data("test-sweep-reads.fq"),
                        "916,11,83.3," +
                        utils.get_test_data("paired-mixed.fq"))

    args = [utils.get_test_data("test-sweep-reads.fq"),
            utils.get_test_data("paired-mixed.fq"),
            '--csv']
    status, out, err = utils.runscript('readstats.py', args)
    assert status == 0

    for k in readstats_output:
        assert k in out, (k, out)


def test_readstats_output():
    readstats_output = ("358 bp / 5 seqs; 71.6 average length",
                        "916 bp / 11 seqs; 83.3 average length")

    outfile = utils.get_temp_filename('output.txt')
    args = ["-o", outfile,
            utils.get_test_data("test-sweep-reads.fq"),
            utils.get_test_data("paired-mixed.fq")]

    status, _, _ = utils.runscript('readstats.py', args)
    assert status == 0

    out = open(outfile).read()

    for k in readstats_output:
        assert k in out, (k, out)


def test_readstats_empty():
    expected_output = "No sequences found in 2 files"

    args = [utils.get_test_data("test-empty.fa"),
            utils.get_test_data("test-empty.fa.bz2")]

    status, out, err = utils.runscript('readstats.py', args)
    assert status == 0

    assert expected_output in out


def test_trim_low_abund_1():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

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
        utils.runscript('trim-low-abund.py', args, in_dir)
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
    utils.runscript('trim-low-abund.py', args, in_dir)

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
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs

    # check for 'quality' string.
    seqs = set([r.quality for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs
    assert '##################' in seqs


# test that the -V option does not trim sequences that are low abundance


def test_trim_low_abund_4_retain_low_abund():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", '-V', infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

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
    utils.runscript('trim-low-abund.py', args, in_dir)

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
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = set([r.sequence for r in screed.open(outfile)])
    assert len(seqs) == 2, seqs

    # untrimmed seq.
    badseq = 'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCgtgCCGCAGCTGTCGTCAGGG' \
             'GATTTCCGGGCGG'
    assert badseq in seqs       # should be there, untrimmed


def test_trim_low_abund_keep_paired():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.paired.fq'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", "-V", infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = [r.name for r in screed.open(outfile)]
    assert seqs[-2:] == ['pair/1', 'pair/2'], seqs


def test_trim_low_abund_keep_paired_casava18():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.paired2.fq'),
                    infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", "-V", infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    seqs = [r.name for r in screed.open(outfile, parse_description=False)]
    assert seqs[-2:] == ['pair:foo 1::N', 'pair:foo 2::N'], seqs


def test_trim_low_abund_highfpr():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.paired.fq'), infile)

    args = ["-k", "17", "-x", "1", "-N", "1", "-V", infile]
    code, out, err = utils.runscript('trim-low-abund.py', args, in_dir,
                                     fail_ok=True)

    assert code == 1
    assert '** ERROR: the graph structure is too small' in err, err


def test_trim_low_abund_trimtest():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.paired.fq'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", "-Z", "2", "-C", "1",
            "-V", infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    for record in screed.open(outfile):
        if record.name == 'seqtrim/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCC'
        elif record.name == 'seqtrim/2':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGC'
        elif record.name == 'seqtrim2/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCA'


def test_trim_low_abund_trimtest_after_load():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    saved_table = utils.get_temp_filename('save.ct')

    shutil.copyfile(utils.get_test_data('test-abund-read-2.paired.fq'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", saved_table, infile]
    utils.runscript('load-into-counting.py', args, in_dir)

    args = ["-Z", "2", "-C", "2", "-V", '--loadtable', saved_table, infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile

    for record in screed.open(outfile):
        if record.name == 'seqtrim/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCC'
        elif record.name == 'seqtrim/2':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGC'
        elif record.name == 'seqtrim2/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCA'


def test_trim_low_abund_trimtest_savetable():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    saved_table = utils.get_temp_filename('save.ct')

    shutil.copyfile(utils.get_test_data('test-abund-read-2.paired.fq'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2",
            "-Z", "2", "-C", "2", "-V", '--savetable', saved_table, infile]
    utils.runscript('trim-low-abund.py', args, in_dir)

    outfile = infile + '.abundtrim'
    assert os.path.exists(outfile), outfile
    assert os.path.exists(saved_table)

    for record in screed.open(outfile):
        if record.name == 'seqtrim/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCC'
        elif record.name == 'seqtrim/2':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGC'
        elif record.name == 'seqtrim2/1':
            print(record.name, record.sequence)
            assert record.sequence == \
                'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCA'

# test that -o/--out option outputs to STDOUT


def test_trim_low_abund_stdout():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2", infile, "-o", "-"]
    _, out, err = utils.runscript('trim-low-abund.py', args, in_dir)

    assert 'GGTTGACGGGGCTCAGGG' in out


def test_roundtrip_casava_format_1():
    # check to make sure that extract-paired-reads produces a file identical
    # to the input file when only paired data is given.

    infile = utils.get_temp_filename('test.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('casava_18-pe.fq'), infile)

    _, out, err = utils.runscript('extract-paired-reads.py', [infile], in_dir)

    r = open(infile).read()

    outfile = infile + '.pe'
    r2 = open(outfile).read()
    assert r == r2, (r, r2)


def test_roundtrip_casava_format_2():
    # check that split-paired-reads -> interleave-reads produces a file
    # identical to input, when only paired reads are given.

    infile = utils.get_temp_filename('test.fq')
    outfile = utils.get_temp_filename('test2.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('casava_18-pe.fq'), infile)

    _, out, err = utils.runscript('split-paired-reads.py', [infile], in_dir)

    utils.runscript('interleave-reads.py', [infile + '.1',
                                            infile + '.2',
                                            '-o', outfile], in_dir)

    r = open(infile).read()
    r2 = open(outfile).read()
    assert r == r2, (r, r2)


def test_existance_failure():
    expected_output = 'ERROR: Input file'

    args = [utils.get_temp_filename('thisfiledoesnotexistatall')]

    status, out, err = utils.runscript(
        'extract-paired-reads.py', args, fail_ok=True)
    assert status == 1

    assert expected_output in err


def test_roundtrip_commented_format():
    """Split/interleave roundtrip for old style format with comments (#873).

    This should produce a file identical to the input when only paired
    reads are given.
    """
    infile = utils.get_temp_filename('test.fq')
    outfile = utils.get_temp_filename('test2.fq')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('old-style-format-w-comments.fq'),
                    infile)

    _, out, err = utils.runscript('split-paired-reads.py', [infile], in_dir)

    utils.runscript('interleave-reads.py', [infile + '.1',
                                            infile + '.2',
                                            '-o', outfile], in_dir)

    r = open(infile).read()
    r2 = open(outfile).read()
    assert r == r2, (r, r2)
