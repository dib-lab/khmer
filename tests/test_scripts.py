#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE.
# Contact: khmer-project@idyll.org
#

from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

# pylint: disable=C0111,C0103,E1103,W0612

import json
import sys
import os
import stat
import shutil
from io import StringIO
import traceback
from nose.plugins.attrib import attr
import threading
import bz2
import gzip
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
    args = ['-x', '1e3', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 94' in err, err
    assert os.path.exists(outfile)


def test_load_into_counting_autoargs_0():
    script = 'load-into-counting.py'

    outfile = utils.get_temp_filename('table')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args = ['-U', '1e7', '--fp-rate', '0.08', outfile, infile]
    (status, out, err) = utils.runscript(script, args)

    assert os.path.exists(outfile)
    assert 'INFO: Overriding default fp 0.1 with new fp: 0.08' in err, err
    assert ' tablesize is too small!' in err, err
    assert 'Estimated FP rate with current config is: 0.9999546' in err, err
    assert 'Recommended tablesize is: 1.77407e+07 bytes' in err, err


def test_load_into_counting_autoargs_1():
    script = 'load-into-counting.py'

    outfile = utils.get_temp_filename('table')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args = ['-U', '1e7', '--max-tablesize', '3e7', outfile, infile]
    (status, out, err) = utils.runscript(script, args)

    assert os.path.exists(outfile)
    assert "Ceiling is: 4.80833e+07 bytes" in err, err
    assert "set memory ceiling automatically." in err, err


def test_load_into_count_graphsize_warning():
    script = 'load-into-counting.py'
    args = ['-k', '20']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert os.path.exists(outfile)
    assert "WARNING: tablesize is default!" in err


def test_load_into_counting_max_memory_usage_parameter():
    script = 'load-into-counting.py'
    args = ['-M', '2e3', '-k', '20']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert os.path.exists(outfile)
    assert "WARNING: tablesize is default!" not in err

    kh = khmer.load_countgraph(outfile)
    assert sum(kh.hashsizes()) < 3e8


def test_load_into_counting_abundance_dist_nobig():
    script = 'load-into-counting.py'
    args = ['-x', '1e3', '-N', '2', '-k', '20', '-b']

    outfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 94' in err, err
    assert os.path.exists(outfile)

    htfile = outfile
    outfile = utils.get_temp_filename('out')
    script2 = 'abundance-dist.py'
    args = ['-z', htfile, infile, outfile]
    (status, out, err) = utils.runscript(script2, args)
    assert 'WARNING: The loaded graph has bigcount' in err, err
    assert 'bigcount' in err, err


def test_load_into_counting_abundance_dist_squashing():
    graphfile = utils.get_temp_filename('out.ct')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args = [graphfile, infile]
    script = 'load-into-counting.py'
    utils.runscript(script, args)

    histogram = utils.get_temp_filename('histogram')
    infile = utils.get_test_data('test-abund-read-2.fa')
    args = [graphfile, infile, histogram]

    script = 'abundance-dist.py'
    # make histogram
    (status, out, err) = utils.runscript(script, args)
    assert os.path.exists(histogram)
    # attempt to overwrite histogram; fail
    failed = True
    try:
        (status, out, err) = utils.runscript(script, args)
        failed = False
    except AssertionError as error:
        assert "exists; not squashing" in str(error), str(error)

    assert failed, "Expected to fail"
    # attempt to overwrite with squashing; should work
    args = ['-s', graphfile, infile, histogram]
    (status, out, err) = utils.runscript(script, args)
    assert "squashing existing file" in err, err

    histfile = open(histogram, 'r')
    lines = histfile.readlines()
    # stripping because boo whitespace
    assert lines[1].strip() == "0,0,0,0.0", lines[1]
    assert lines[2].strip() == "1,83,83,1.0", lines[2]


def test_load_into_counting_nonwritable():
    script = 'load-into-counting.py'
    args = ['-x', '1e3', '-N', '2', '-k', '20']

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
    args = ['-x', '1e12', '-N', '2', '-k', '20', '--force']

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
    args = ['-x', '1e7', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out.kh')
    infile = utils.get_test_data('test-abund-read-2.fa')

    args.extend([outfile, infile, infile, infile, infile, infile,
                 infile, infile, infile, infile, infile, infile])

    (status, out, err) = utils.runscript(script, args)
    assert 'Total number of unique k-mers: 95' in err, err
    assert os.path.exists(outfile)


def test_load_into_counting_tsv():
    script = 'load-into-counting.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-s', 'tsv']

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
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-s', 'json']

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
        u"files": [infile],
        u"ht_name": outbase,
        u"num_kmers": 95,
        u"num_reads": 1001,
        u"fpr": 9.025048735197377e-11,
        u"mrinfo_version": "0.2.0",
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


def test_filter_abund_2_stdin():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = ['-C', '1', counting_ht, '-']
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status == 1
    assert "Accepting input from stdin; output filename must be provided" \
           in str(err)

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

    seqs = set([r.name for r in screed.open(outfile)])
    assert 'pair:foo 1::N' in seqs, seqs


def test_filter_abund_1_singlefile():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'filter-abund-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', infile]
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
    tabfile = utils.get_temp_filename('test-savegraph.ct')

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'filter-abund-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '--savegraph',
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

    seqs = set([r.name for r in screed.open(outfile)])
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

    # copy test file over to test.fq & load into countgraph
    shutil.copyfile(utils.get_test_data('test-filter-abund-Ns.fq'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = ['-C', '3', counting_ht, infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    # test for a sequence with an 'N' in it --
    names = set([r.name for r in screed.open(outfile)])
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

    # copy test file over to test.fq & load into countgraph
    shutil.copyfile(utils.get_test_data('test-filter-abund-Ns.fq'), infile)

    script = 'filter-abund-single.py'
    args = ['-k', '17', '-x', '1e7', '-N', '2', '-C', '3', infile]
    utils.runscript(script, args, in_dir)

    outfile = infile + '.abundfilt'
    assert os.path.exists(outfile), outfile

    # test for a sequence with an 'N' in it --
    names = set([r.name for r in screed.open(outfile)])
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
    kh = khmer._Nodegraph(K, [1])
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
    kh = khmer._Nodegraph(K, [1])
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
    names = [r.name for r in screed.open(outfile)]
    names = set(names)
    assert 'seq 1::BAR' in names


def test_count_median():
    infile = utils.get_temp_filename('test.fa')
    outfile = infile + '.counts'

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = 'count-median.py'
    args = [counting_ht, infile, outfile]
    utils.runscript(script, args)

    assert os.path.exists(outfile), outfile

    data = [x.strip() for x in open(outfile).readlines()[1:]]
    data = set(data)
    assert len(data) == 2, data
    assert 'seq,1001,1001.0,0.0,18' in data, data
    assert '895:1:37:17593:9954/1,1,103.803741455,303.702941895,114' in data


def test_count_median_fq_csv():
    infile = utils.get_temp_filename('test.fq')
    outfile = infile + '.counts'

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = 'count-median.py'
    args = [counting_ht, infile, outfile]
    utils.runscript(script, args)

    assert os.path.exists(outfile), outfile

    data = [x.strip() for x in open(outfile)]
    data = set(data)
    assert len(data) == 4, data
    assert 'name,median,average,stddev,seqlen' in data
    assert 'seq,1001,1001.0,0.0,18' in data

    # verify that sequence names remain unparsed
    names = set([line.split(',')[0] for line in data])
    assert '895:1:37:17593:9954 1::FOO' in names, names


def test_count_median_fq_csv_stdout():
    infile = utils.get_temp_filename('test.fq')
    outfile = '-'

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fq'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = 'count-median.py'
    args = [counting_ht, infile, outfile]
    (status, out, err) = utils.runscript(script, args)

    assert 'name,median,average,stddev,seqlen' in out
    assert 'seq,1001,1001.0,0.0,18' in out


def test_load_graph():
    script = 'load-graph.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    try:
        ht = khmer.load_nodegraph(ht_file)
    except OSError as err:
        assert 0, str(err)
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

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    ht = khmer.load_nodegraph(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph worked :)
    subset = ht.do_subset_partition(0, 0)
    x = ht.subset_count_partitions(subset)
    assert x == (1, 0), x


def test_oxli_build_graph_unique_kmers_arg():
    script = 'oxli'
    args = ['build-graph', '-x', '1e7', '-N', '2', '-k', '20', '-U', '3960']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of unique k-mers: 3960' in err, err
    assert 'INFO: set memory ceiling automatically' in err, err
    assert 'Ceiling is: 1e+06 bytes' in err, err

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert os.path.exists(tagset_file), tagset_file

    ht = khmer.load_nodegraph(ht_file)
    ht.load_tagset(tagset_file)

    # check to make sure we get the expected result for this data set
    # upon partitioning (all in one partition).  This is kind of a
    # roundabout way of checking that load-graph worked :)
    subset = ht.do_subset_partition(0, 0)
    x = ht.subset_count_partitions(subset)
    assert x == (1, 0), x


def test_oxli_nocommand():
    script = 'oxli'

    (status, out, err) = utils.runscript(script, [])
    assert status == 0


def test_load_graph_no_tags():
    script = 'load-graph.py'
    args = ['-x', '1e7', '-N', '2', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert not os.path.exists(tagset_file), tagset_file

    assert khmer.load_nodegraph(ht_file)

    # can't think of a good way to make sure this worked, beyond just
    # loading the ht file...


def test_oxli_build_graph_no_tags():
    script = 'oxli'
    args = ['build-graph', '-x', '1e7', '-N', '2', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    utils.runscript(script, args)

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    tagset_file = outfile + '.tagset'
    assert not os.path.exists(tagset_file), tagset_file

    assert khmer.load_nodegraph(ht_file)

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

    ht_file = outfile
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

    ht_file = outfile
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


def test_load_graph_max_memory_usage_parameter():
    script = 'load-graph.py'
    args = ['-M', '2e7', '-k', '20', '-n']

    outfile = utils.get_temp_filename('out')
    infile = utils.get_test_data('random-20-a.fa')

    args.extend([outfile, infile])

    (status, out, err) = utils.runscript(script, args)

    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = outfile
    assert os.path.exists(ht_file), ht_file

    try:
        ht = khmer.load_nodegraph(ht_file)
    except OSError as err:
        assert 0, str(err)

    assert (sum(ht.hashsizes()) / 8.) < 2e7, ht.hashsizes()


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

    ht_file = outfile
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

    ht = khmer.load_nodegraph(graphbase)
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

    ht = khmer.load_nodegraph(graphbase)
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (99, 0), x          # should be 99 partitions at K=21


def test_partition_graph_nojoin_stoptags():
    # test with stoptags
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    # add in some stop tags
    ht = khmer.load_nodegraph(graphbase)
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

    ht = khmer.load_nodegraph(graphbase)
    ht.load_tagset(graphbase + '.tagset')
    ht.load_partitionmap(final_pmap_file)

    x = ht.count_partitions()
    assert x == (2, 0), x          # should be 2 partitions


def test_partition_graph_big_traverse():
    graphbase = _make_graph(utils.get_test_data('biglump-random-20-a.fa'),
                            do_partition=True, stop_big_traverse=False)

    final_pmap_file = graphbase + '.pmap.merged'
    assert os.path.exists(final_pmap_file)

    ht = khmer.load_nodegraph(graphbase)
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

    ht = khmer.load_nodegraph(graphbase)
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


def test_partition_find_knots_existing_stoptags():
    graphbase = _make_graph(utils.get_test_data('random-20-a.fa'))

    script = 'partition-graph.py'
    args = [graphbase]
    utils.runscript(script, args)

    script = 'make-initial-stoptags.py'
    args = [graphbase]
    utils.runscript(script, args)

    script = 'find-knots.py'
    args = [graphbase]
    (status, out, err) = utils.runscript(script, args)

    stoptags_file = graphbase + '.stoptags'
    assert os.path.exists(stoptags_file)
    assert "loading stoptags" in err, err
    assert "these output stoptags will include the already" in err, err


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
    assert dist.strip() == '1 11960 11960 11960', dist.strip()

    parts = [r.name.split('\t')[1]
             for r in screed.open(partfile)]
    assert len(parts) == 13538, len(parts)
    parts = set(parts)
    assert len(parts) == 12602, len(parts)


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

    screed_iter = screed.open(partfile)
    names = [r.name.split('\t')[0] for r in screed_iter]
    assert '35 1::FOO' in names
    assert '46 1::FIZ' in names

    screed_iter = screed.open(partfile)
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
    status, out, err = utils.runscript(script, args, in_dir)
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

    status, _, err = utils.runscript(script, args, in_dir, fail_ok=True)
    assert "ERROR: Input file", "is empty; Exiting." in err
    assert status != 0
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
        assert (line == 'abundance,count,cumulative,cumulative_fraction'), line
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_stdout():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    htfile = _make_counting(infile, K=17)

    script = 'abundance-dist.py'
    args = ['-z', htfile, infile, "-"]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert '1,96,96,0.98' in out, out
    assert '1001,2,98,1.0' in out, out


def test_abundance_dist_nobigcount():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    htfile = _make_counting(infile, K=17)

    script = 'abundance-dist.py'
    args = ['-b', '-z', htfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '255,2,98,1.0', line


def test_abundance_dist_threaded():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '--threads', '18',
            infile, outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 98' in err, err

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_single_csv():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', infile,
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
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '255,2,98,1.0', line


def test_abundance_dist_single_nosquash():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_abundance_dist_single_savegraph():
    infile = utils.get_temp_filename('test.fa')
    outfile = utils.get_temp_filename('test.dist')
    tabfile = utils.get_temp_filename('test-savegraph.ct')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'abundance-dist-single.py'
    args = ['-x', '1e7', '-N', '2', '-k', '17', '-z', '--savegraph',
            tabfile, infile, outfile]
    utils.runscript(script, args, in_dir)

    with open(outfile) as fp:
        line = fp.readline().strip()    # skip header
        line = fp.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fp.readline().strip()
        assert line == '1001,2,98,1.0', line


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

    screed_iter = screed.open(partfile)
    names = [r.name.split('\t')[0] for r in screed_iter]
    assert '35 1::FOO' in names
    assert '46 1::FIZ' in names


def test_interleave_read_stdout():
    # create input files
    infile1 = utils.get_test_data('paired-slash1.fq.1')
    infile2 = utils.get_test_data('paired-slash1.fq.2')

    # correct output
    ex_outfile = utils.get_test_data('paired-slash1.fq')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, infile2]

    (stats, out, err) = utils.runscript(script, args)

    with open(outfile, 'w') as ofile:
        ofile.write(out)

    n = 0
    for r, q in zip(screed.open(ex_outfile), screed.open(outfile)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0


def test_interleave_read_seq1_fq():
    # create input files
    infile1 = utils.get_test_data('paired-slash1.fq.1')
    infile2 = utils.get_test_data('paired-slash1.fq.2')

    # correct output
    ex_outfile = utils.get_test_data('paired-slash1.fq')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, infile2, '-o', outfile]

    utils.runscript(script, args)

    n = 0
    for r, q in zip(screed.open(ex_outfile), screed.open(outfile)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0


def test_interleave_read_badleft_badright():
    # create input files
    infile1 = utils.get_test_data('paired-broken.fq.badleft')
    infile2 = utils.get_test_data('paired-broken.fq.badright')

    # correct output
    ex_outfile = utils.get_test_data('paired-broken.fq.paired_bad')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, infile2, '-o', outfile]

    utils.runscript(script, args)

    n = 0
    for r, q in zip(screed.open(ex_outfile), screed.open(outfile)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
    assert n > 0


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


def test_interleave_reads_broken_fq_5():
    # test input files
    infile1 = utils.get_test_data('paired-broken4.fq.1')
    infile2 = utils.get_test_data('paired-broken4.fq.2')

    # actual output file
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, infile2, '-o', outfile]

    status, out, err = utils.runscript(script, args, fail_ok=True)
    assert status == 1
    assert "ERROR: This doesn't look like paired data!" in err


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


def test_make_initial_stoptags_load_stoptags():
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
    args = ['test-reads', '--stoptags', 'test-reads.stoptags']
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
    for r, q in zip(screed.open(ex_outfile1),
                    screed.open(outfile1)):
        n += 1
        assert r.name == q.name, (r.name, q.name, n)
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0

    n = 0
    for r, q in zip(screed.open(ex_outfile2),
                    screed.open(outfile2)):
        n += 1
        assert r.name == q.name
        assert r.sequence == q.sequence
        assert r.quality == q.quality
    assert n > 0


def test_extract_paired_reads_3_output_dir():
    # test input file
    infile = utils.get_test_data('paired-mixed.fa')

    ex_outfile1 = utils.get_test_data('paired-mixed.fa.pe')
    ex_outfile2 = utils.get_test_data('paired-mixed.fa.se')

    # output directory
    out_dir = utils.get_temp_filename('output')

    script = 'extract-paired-reads.py'
    args = [infile, '-d', out_dir]

    utils.runscript(script, args)

    outfile1 = os.path.join(out_dir, 'paired-mixed.fa.pe')
    outfile2 = os.path.join(out_dir, 'paired-mixed.fa.se')
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


def test_extract_paired_reads_4_output_files():
    # test input file
    infile = utils.get_test_data('paired-mixed.fa')

    ex_outfile1 = utils.get_test_data('paired-mixed.fa.pe')
    ex_outfile2 = utils.get_test_data('paired-mixed.fa.se')

    # actual output files...
    outfile1 = utils.get_temp_filename('out_pe')
    outfile2 = utils.get_temp_filename('out_se')

    script = 'extract-paired-reads.py'
    args = [infile, '-p', outfile1, '-s', outfile2]

    utils.runscript(script, args)

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


def test_extract_paired_reads_5_stdin_error():
    script = 'extract-paired-reads.py'
    args = ['-f', '/dev/stdin']

    status, out, err = utils.runscript(script, args, fail_ok=True)
    assert status == 1
    assert "output filenames must be provided." in err


def execute_extract_paired_streaming(ifilename):
    fifo = utils.get_temp_filename('fifo')
    in_dir = os.path.dirname(fifo)
    outfile1 = utils.get_temp_filename('paired.pe')
    outfile2 = utils.get_temp_filename('paired.se')
    script = 'extract-paired-reads.py'
    args = [fifo, '-p', outfile1, '-s', outfile2]

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


def test_extract_paired_streaming():
    testinput = utils.get_test_data('paired-mixed.fa')
    o = execute_extract_paired_streaming(testinput)


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
    args = [infile]

    status, out, err = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status == 1, status
    assert "Unpaired reads found" in err


def test_split_paired_reads_2_stdin_no_out():
    script = 'split-paired-reads.py'
    args = ['-']

    status, out, err = utils.runscript(script, args, fail_ok=True)
    assert status == 1
    assert "Accepting input from stdin; output filenames must " in err


def test_split_paired_reads_2_mixed_fq():
    # test input file
    infile = utils.get_temp_filename('test.fq')
    shutil.copyfile(utils.get_test_data('paired-mixed-2.fq'), infile)
    in_dir = os.path.dirname(infile)

    script = 'split-paired-reads.py'
    args = ['-0', '/dev/null', infile]

    status, out, err = utils.runscript(script, args, in_dir)
    assert status == 0
    assert "split 6 sequences (3 left, 3 right, 5 orphans)" in err, err


def test_split_paired_reads_2_mixed_fq_orphans_to_file():
    # test input file
    infile = utils.get_temp_filename('test.fq')
    shutil.copyfile(utils.get_test_data('paired-mixed-2.fq'), infile)
    in_dir = os.path.dirname(infile)
    outfile = utils.get_temp_filename('out.fq')

    script = 'split-paired-reads.py'
    args = ['-0', outfile, infile]

    status, out, err = utils.runscript(script, args, in_dir)
    assert status == 0
    assert "split 6 sequences (3 left, 3 right, 5 orphans)" in err, err

    n_orphans = len([1 for record in screed.open(outfile)])
    assert n_orphans == 5
    n_left = len([1 for record in screed.open(infile + '.1')])
    assert n_left == 3
    n_right = len([1 for record in screed.open(infile + '.2')])
    assert n_right == 3
    for filename in [outfile, infile + '.1', infile + '.2']:
        fp = gzip.open(filename)
        try:
            fp.read()
        except IOError as e:
            assert "Not a gzipped file" in str(e), str(e)
        fp.close()


def test_split_paired_reads_2_mixed_fq_gzfile():
    # test input file
    infile = utils.get_temp_filename('test.fq')
    shutil.copyfile(utils.get_test_data('paired-mixed-2.fq'), infile)
    in_dir = os.path.dirname(infile)
    outfile = utils.get_temp_filename('out.fq')

    script = 'split-paired-reads.py'
    args = ['-0', outfile, '--gzip', infile]

    status, out, err = utils.runscript(script, args, in_dir)
    assert status == 0
    assert "split 6 sequences (3 left, 3 right, 5 orphans)" in err, err

    n_orphans = len([1 for record in screed.open(outfile)])
    assert n_orphans == 5
    n_left = len([1 for record in screed.open(infile + '.1')])
    assert n_left == 3
    n_right = len([1 for record in screed.open(infile + '.2')])
    assert n_right == 3

    for filename in [outfile, infile + '.1', infile + '.2']:
        fp = gzip.open(filename)
        fp.read()                       # this will fail if not gzip file.
        fp.close()


def test_split_paired_reads_2_mixed_fq_broken_pairing_format():
    # test input file
    infile = utils.get_temp_filename('test.fq')
    shutil.copyfile(utils.get_test_data('paired-mixed-broken.fq'), infile)
    in_dir = os.path.dirname(infile)

    script = 'split-paired-reads.py'
    args = [infile]

    status, out, err = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status == 1
    assert "Unpaired reads found starting at 895:1:37:17593:9954" in err, err


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
    args = ['-d', output_dir, '-1', outfile1, infile]

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
    args = ['-2', outfile2, '-d', output_dir, infile]

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


def test_sample_reads_randomly_force_single_outfile():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-reads.fa'), infile)

    script = 'sample-reads-randomly.py'
    # fix random number seed for reproducibility
    args = ['-N', '10', '-M', '12000', '-R', '1', '--force_single', '-o',
            in_dir + '/randreads.out']

    args.append(infile)
    utils.runscript(script, args, in_dir)

    outfile = in_dir + '/randreads.out'
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

    seqs = set([r.name for r in screed.open(outfile)])
    print(list(sorted(seqs)))
    assert seqs == answer


def test_sample_reads_randomly_stdin_no_out():
    script = 'sample-reads-randomly.py'
    args = ['-']

    (status, out, err) = utils.runscript(script, args, fail_ok=True)
    assert status != 0
    assert "Accepting input from stdin; output filename" in err, err


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
    assert len(out.splitlines()) == 0, len(out.splitlines())
    assert "No lines dropped" in err

    names = [r.name for r in screed.open(clean_outfile)]
    assert '895:1:1:1246:14654 1:N:0:NNNNN' in names, names

    args = [n_infile, '-n', '-o', n_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) == 0
    assert "No lines dropped" in err

    args = [clean_infile, '-o', clean_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 0
    assert "0 lines dropped" in err

    args = [n_infile, '-o', n_outfile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) == 0, out
    assert "4 lines dropped" in err, err

    args = [clean_infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) > 0
    assert "0 lines dropped" in err

    args = [n_infile]
    (status, out, err) = utils.runscript(script, args, in_dir_n)
    assert len(out.splitlines()) > 0
    assert "4 lines dropped" in err

    args = [clean_infile, '-o', clean_outfile, '--gzip']
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 0
    assert "0 lines dropped" in err

    args = [clean_infile, '-o', clean_outfile, '--bzip']
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert len(out.splitlines()) == 0
    assert "0 lines dropped" in err


def test_fastq_to_fasta_streaming_compressed_gzip():

    script = 'fastq-to-fasta.py'
    infile = utils.get_temp_filename('test-clean.fq')
    in_dir = os.path.dirname(infile)
    fifo = utils.get_temp_filename('fifo')
    copyfilepath = utils.get_temp_filename('copied.fa.gz', in_dir)
    shutil.copyfile(utils.get_test_data('test-reads.fq.gz'), infile)

    # make a fifo to simulate streaming
    os.mkfifo(fifo)
    args = ['--gzip', '-o', fifo, infile]
    # FIFOs MUST BE OPENED FOR READING BEFORE THEY ARE WRITTEN TO
    # If this isn't done, they will BLOCK and things will hang.
    thread = threading.Thread(target=utils.runscript,
                              args=(script, args, in_dir))
    thread.start()
    copyfile = io.open(copyfilepath, 'wb')
    fifofile = io.open(fifo, 'rb')

    # read binary to handle compressed files
    chunk = fifofile.read(8192)
    while len(chunk) > 0:
        copyfile.write(chunk)
        chunk = fifofile.read(8192)

    fifofile.close()
    thread.join()
    copyfile.close()

    # verify that the seqs are there and not broken
    f = screed.open(copyfilepath)
    count = 0
    for _ in f:
        count += 1

    assert count == 25000, count
    f.close()

    # verify we're looking at a gzipped file
    gzfile = io.open(file=copyfilepath, mode='rb', buffering=8192)
    magic = b"\x1f\x8b\x08"  # gzip magic signature
    file_start = gzfile.peek(len(magic))
    assert file_start[:3] == magic, file_start[:3]


def test_fastq_to_fasta_streaming_compressed_bzip():

    script = 'fastq-to-fasta.py'
    infile = utils.get_temp_filename('test-clean.fq')
    in_dir = os.path.dirname(infile)
    fifo = utils.get_temp_filename('fifo')
    copyfilepath = utils.get_temp_filename('copied.fa.bz', in_dir)
    shutil.copyfile(utils.get_test_data('test-reads.fq.gz'), infile)

    # make a fifo to simulate streaming
    os.mkfifo(fifo)
    args = ['--bzip', '-o', fifo, infile]
    # FIFOs MUST BE OPENED FOR READING BEFORE THEY ARE WRITTEN TO
    # If this isn't done, they will BLOCK and things will hang.
    thread = threading.Thread(target=utils.runscript,
                              args=(script, args, in_dir))
    thread.start()
    copyfile = io.open(copyfilepath, 'wb')
    fifofile = io.open(fifo, 'rb')

    # read binary to handle compressed files
    chunk = fifofile.read(8192)
    while len(chunk) > 0:
        copyfile.write(chunk)
        chunk = fifofile.read(8192)

    fifofile.close()
    thread.join()
    copyfile.close()

    # verify that the seqs are there and not broken
    f = screed.open(copyfilepath)
    count = 0
    for _ in f:
        count += 1

    assert count == 25000, count
    f.close()

    # verify we're looking at a gzipped file
    bzfile = io.open(file=copyfilepath, mode='rb', buffering=8192)
    magic = b"\x42\x5a\x68"  # bzip magic signature
    file_start = bzfile.peek(len(magic))
    assert file_start[:3] == magic, file_start[:3]


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

    names = [r.name for r in screed.open(fa_outfile)]
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

    names = [r.name for r in screed.open(fq_outfile)]
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


def _execute_load_graph_streaming(filename):
    '''Helper function for the matrix of streaming tests using screed via
    filter-abund-single, i.e. uncompressed fasta, gzip fasta, bz2 fasta,
    uncompressed fastq, etc.
    This is not directly executed but is run by the tests themselves
    '''

    scripts = utils.scriptpath()
    infile = utils.get_temp_filename('temp')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data(filename), infile)

    args = '-x 1e7 -N 2 -k 20 out -'

    cmd = 'cat {infile} | {scripts}/load-graph.py {args}'.format(
        infile=infile, scripts=scripts, args=args)

    (status, out, err) = utils.run_shell_cmd(cmd, in_directory=in_dir)

    if status != 0:
        print(out)
        print(err)
        assert status == 0, status

    assert 'Total number of unique k-mers: 3960' in err, err

    ht_file = os.path.join(in_dir, 'out')
    assert os.path.exists(ht_file), ht_file

    tagset_file = os.path.join(in_dir, 'out.tagset')
    assert os.path.exists(tagset_file), tagset_file

    ht = khmer.load_nodegraph(ht_file)
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
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fa'))


def test_read_parser_streaming_ufq():
    # uncompressed FASTQ
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fq'))


@attr('known_failing')
def test_read_parser_streaming_bzfq():
    # bzip compressed FASTQ
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fq.bz2'))


def test_read_parser_streaming_gzfq():
    # gzip compressed FASTQ
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fq.gz'))


@attr('known_failing')
def test_read_parser_streaming_bzfa():
    # bzip compressed FASTA
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fa.bz2'))


def test_read_parser_streaming_gzfa():
    # gzip compressed FASTA
    _execute_load_graph_streaming(utils.get_test_data('random-20-a.fa.gz'))


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
    (status, out, err) = utils.runscript('trim-low-abund.py', args, in_dir,
                                         fail_ok=True)
    assert status == 1
    assert "Error: Cannot input the same filename multiple times." in str(err)


def test_trim_low_abund_1_stdin_err():
    args = ["-"]

    (status, out, err) = utils.runscript('trim-low-abund.py', args,
                                         fail_ok=True)
    assert status == 1
    assert "Accepting input from stdin; output filename must be provided" \
           in str(err)


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

    seqs = [r.name for r in screed.open(outfile)]
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

    args = ["-Z", "2", "-C", "2", "-V", '--loadgraph', saved_table, infile]
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


def test_trim_low_abund_trimtest_savegraph():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    saved_table = utils.get_temp_filename('save.ct')

    shutil.copyfile(utils.get_test_data('test-abund-read-2.paired.fq'), infile)

    args = ["-k", "17", "-x", "1e7", "-N", "2",
            "-Z", "2", "-C", "2", "-V", '--savegraph', saved_table, infile]
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


def test_unique_kmers_defaults():
    infile = utils.get_temp_filename('random-20-a.fa')
    shutil.copyfile(utils.get_test_data('random-20-a.fa'), infile)

    args = ['-k', '20', '-e', '0.01', infile]

    _, out, err = utils.runscript('unique-kmers.py', args,
                                  os.path.dirname(infile))

    err = err.splitlines()
    assert ('Estimated number of unique 20-mers in {0}: 3950'.format(infile)
            in err)
    assert 'Total estimated number of unique 20-mers: 3950' in err


def test_unique_kmers_streaming():
    infile = utils.get_temp_filename('random-20-a.fa')
    shutil.copyfile(utils.get_test_data('random-20-a.fa'), infile)

    args = ['-k', '20', '-e', '0.01', '--stream-out', infile]

    _, out, err = utils.runscript('unique-kmers.py', args,
                                  os.path.dirname(infile))

    err = err.splitlines()
    assert ('Estimated number of unique 20-mers in {0}: 3950'.format(infile)
            in err)
    assert 'Total estimated number of unique 20-mers: 3950' in err
    assert "CCAACCATGGTAGGTTAGGAAAGCCGCCAAATAAGTTCTTATACG" in out, out


def test_unique_kmers_report_fp():
    infile = utils.get_temp_filename('random-20-a.fa')
    shutil.copyfile(utils.get_test_data('random-20-a.fa'), infile)
    outfile = utils.get_temp_filename('report.unique')

    args = ['-k', '20', '-e', '0.01', '-R', outfile, infile]

    _, out, err = utils.runscript('unique-kmers.py', args,
                                  os.path.dirname(infile))

    err = err.splitlines()
    assert ('Estimated number of unique 20-mers in {0}: 3950'.format(infile)
            in err)
    assert 'Total estimated number of unique 20-mers: 3950' in err

    with open(outfile, 'r') as report_fp:
        outf = report_fp.read().splitlines()
        assert '3950 20 (total)' in outf
        assert '3950 20 total' in outf


def test_unique_kmers_diagnostics():
    infile = utils.get_temp_filename('random-20-a.fa')
    shutil.copyfile(utils.get_test_data('random-20-a.fa'), infile)

    args = ['-k', '20', '-e', '0.01', '--diagnostics', infile]

    _, out, err = utils.runscript('unique-kmers.py', args,
                                  os.path.dirname(infile))

    out = out.splitlines()
    assert ('expected_fp\tnumber_hashtable(Z)\t'
            'size_hashtable(H)\texpected_memory_usage' in err)


def test_unique_kmers_multiple_inputs():
    infiles = []
    for fname in ('random-20-a.fa', 'paired-mixed.fa'):
        infile = utils.get_temp_filename(fname)
        shutil.copyfile(utils.get_test_data(fname), infile)
        infiles.append(infile)

    args = ['-k', '20', '-e', '0.01']
    args += infiles

    _, out, err = utils.runscript('unique-kmers.py', args,
                                  os.path.dirname(infile))

    err = err.splitlines()
    assert ('Estimated number of unique 20-mers in {0}: 3950'
            .format(infiles[0]) in err)
    assert ('Estimated number of unique 20-mers in {0}: 232'.format(infiles[1])
            in err)
    assert 'Total estimated number of unique 20-mers: 4170' in err
