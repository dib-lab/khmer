# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2015, The Regents of the University of California.
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
# pylint: disable=missing-docstring,invalid-name,unused-variable
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
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
    args = ['--loadgraph', hashfile, '-o', outfile, infile]
    (status, out, err) = utils.runscript(script, args)
    assert status == 0, (out, err)
    assert os.path.exists(outfile)


def test_normalize_by_median_loadgraph_with_args():
    infile = utils.get_test_data("test-abund-read-2.fa")
    tablefile = utils.get_temp_filename("table")
    in_dir = os.path.dirname(tablefile)

    script = "load-into-counting.py"
    args = [tablefile, infile]
    (status, out, err) = utils.runscript(script, args)

    script = "normalize-by-median.py"
    args = ["--ksize", "7", "--loadgraph", tablefile, infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert 'WARNING: You are loading a saved k-mer countgraph from' in err, err


def test_normalize_by_median_empty_file():
    infile = utils.get_temp_filename('empty')
    shutil.copyfile(utils.get_test_data('empty-file'), infile)
    script = 'normalize-by-median.py'
    in_dir = os.path.dirname(infile)

    args = [infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'WARNING:' in err, err
    assert 'is empty' in err, err
    assert 'SKIPPED' in err, err


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
    assert "I/O Errors" not in err


def test_normalize_by_median_quiet():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', '--quiet', '-M', '2e6', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert len(out) == 0, out
    assert len(err) < 460, len(err)

    outfile = infile + '.keep'
    assert os.path.exists(outfile), outfile

    seqs = [r.sequence for r in screed.open(outfile)]
    assert len(seqs) == 1, seqs
    assert seqs[0].startswith('GGTTGACGGGGCTCAGGGGG'), seqs
    assert "I/O Errors" not in err


def test_normalize_by_median_unpaired_final_read():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('single-read.fq'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', '-p', infile]
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status != 0
    assert "ERROR: Unpaired reads when require_paired" in err, err


def test_normalize_by_median_sanity_check_0():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('single-read.fq'), infile)

    script = 'normalize-by-median.py'
    args = ['-U', '1024', '--max-mem', '60', infile]
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status != 0, status
    assert "recommended false positive ceiling of 0.1!" in err, err


def test_normalize_by_median_sanity_check_1():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-filter-abund-Ns.fq'), infile)

    script = 'normalize-by-median.py'
    args = ['-U', '83', '--max-tablesize', '17', infile]
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status != 0
    assert "Warning: The given tablesize is too small!" in err, err


def test_normalize_by_median_sanity_check_2():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-filter-abund-Ns.fq'), infile)

    script = 'normalize-by-median.py'
    args = ['-U', '83', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert "*** INFO: set memory ceiling automatically." in err, err
    assert "*** Ceiling is: 1e+06 bytes" in err, err


def test_normalize_by_median_sanity_check_3():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    tablefile = utils.get_temp_filename('table', in_dir)

    shutil.copyfile(utils.get_test_data('test-filter-abund-Ns.fq'), infile)

    script = 'normalize-by-median.py'
    args = ['-s', tablefile, '-U', '83', '--fp-rate', '0.7', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert "Overriding default fp 0.1 with new fp: 0.7" in err, err

    args = ['--loadgraph', tablefile, '-U', '83', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert "WARNING: You have asked that the graph size be auto" in err, err
    assert "NOT be set automatically" in err, err
    assert "loading an existing graph" in err, err


def test_normalize_by_median_unforced_badfile():
    CUTOFF = '1'

    infile = utils.get_temp_filename("potatoes")
    outfile = infile + '.keep'
    in_dir = os.path.dirname(infile)
    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', infile]
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status != 0
    assert "ERROR: [Errno 2] No such file or directory:" in err, err

    if os.path.exists(outfile):
        assert False, '.keep file should have been removed: '


def test_normalize_by_median_contradictory_args():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    outfile = utils.get_temp_filename('report.out')

    shutil.copyfile(utils.get_test_data('test-large.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', '1', '-k', '17', '--force_single', '-p', '-R',
            outfile, infile]
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status != 0
    assert "cannot both be set" in err, err


def test_normalize_by_median_stdout_3():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', infile, '--out', '-']
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert 'Total number of unique k-mers: 98' in err, err
    assert 'in block device' in err, err
    assert "I/O Errors" not in err


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
    # this tests basic reporting of diginorm stats => report.out, including
    # a test of aggregate stats for two input files.

    infile = utils.get_temp_filename('test.fa')
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    infile2 = utils.get_temp_filename('test2.fa')
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile2)

    in_dir = os.path.dirname(infile)
    outfile = utils.get_temp_filename('report.out')

    script = 'normalize-by-median.py'
    args = ['-C', '1', '-k', '17', '-R', outfile, infile, infile2]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile)
    report = open(outfile, 'r')
    line = report.readline().strip()
    assert line == 'total,kept,f_kept', line
    line = report.readline().strip()
    assert line == '1001,1,0.000999', line
    line = report.readline().strip()
    assert line == '2002,1,0.0004995', line


def test_normalize_by_median_report_fp_hifreq():
    # this tests high-frequency reporting of diginorm stats for a single
    # file => report.out.

    infile = utils.get_temp_filename('test.fa')
    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    in_dir = os.path.dirname(infile)
    outfile = utils.get_temp_filename('report.out')

    script = 'normalize-by-median.py'
    args = ['-C', '1', '-k', '17', '-R', outfile, infile,
            '--report-frequency', '100']
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert os.path.exists(outfile)
    report = open(outfile, 'r')
    line = report.readline().strip()
    assert line == 'total,kept,f_kept', line
    line = report.readline().strip()
    assert line == '100,1,0.01', line
    line = report.readline().strip()
    assert line == '200,1,0.005', line


@attr('huge')
def test_normalize_by_median_report_fp_huge():
    # this tests reporting of diginorm stats => report.out for a large
    # file, with the default reporting interval of once every 100k.

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    outfile = utils.get_temp_filename('report.out')

    shutil.copyfile(utils.get_test_data('test-large.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', '1', '-k', '17', '-R', outfile, infile]
    (status, out, err) = utils.runscript(script, args, in_dir)

    assert "fp rate estimated to be 0.623" in err, err
    report = open(outfile, 'r')
    line = report.readline()            # skip header
    line = report.readline()
    assert "100000,25261,0.2526" in line, line


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

    assert 'Total number of unique k-mers: 4061' in err, err

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

    args = ['-C', CUTOFF, '-k', '17', '--force_single', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert 'Total number of unique k-mers: 98' in err, err
    assert 'kept 1 of 2 or 50.0%' in err, err

    args = ['-C', CUTOFF, '-k', '17', '-p', infile]
    (status, out, err) = utils.runscript(script, args, in_dir)
    assert 'Total number of unique k-mers: 99' in err, err
    assert 'kept 2 of 2 or 100.0%' in err, err


def test_normalize_by_median_double_file_name():
    infile = utils.get_temp_filename('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)

    script = 'normalize-by-median.py'
    args = [utils.get_test_data('test-abund-read-2.fa'), infile]

    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status != 0
    assert "Duplicate filename--Cannot handle this!" in err, err


def test_normalize_by_median_stdin_no_out():
    infile = utils.get_temp_filename('test-abund-read-2.fa')
    in_dir = os.path.dirname(infile)

    script = 'normalize-by-median.py'
    args = ["-"]

    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status != 0
    assert "Accepting input from stdin; output filename" in err, err


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

    names = [r.name for r in screed.open(outfile)]
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
    status, out, err = utils.runscript(script, args, in_dir, fail_ok=True)
    status != 0
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
    assert '** I/O Errors' in err


def test_normalize_by_median_no_bigcount():
    infile = utils.get_temp_filename('test.fa')
    hashfile = utils.get_temp_filename('test-out.ct')
    outfile = infile + '.keep'
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=8)

    script = 'normalize-by-median.py'
    args = ['-C', '1000', '-k 8', '--savegraph', hashfile, infile]

    (status, out, err) = utils.runscript(script, args, in_dir)
    assert status == 0, (out, err)
    print((out, err))

    assert os.path.exists(hashfile), hashfile
    kh = khmer.load_countgraph(hashfile)

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


def test_normalize_by_median_emptycountgraph():
    CUTOFF = '1'

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('test-empty.fa'), infile)

    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '--loadgraph', infile, infile]
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status != 0
    assert 'ValueError' in err, (status, out, err)


def test_normalize_by_median_fpr():
    MAX_TABLESIZE_PARAM = 12

    infile = utils.get_temp_filename('test-fpr.fq')
    in_dir = os.path.dirname(infile)
    shutil.copyfile(utils.get_test_data('test-fastq-reads.fq'), infile)

    script = 'normalize-by-median.py'
    args = ['-f', '-k 17', '-x ' + str(MAX_TABLESIZE_PARAM), infile]

    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert status != 0
    assert os.path.exists(infile + '.keep'), infile
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


def test_normalize_by_median_streaming_0():
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


def test_normalize_by_median_streaming_1():
    CUTOFF = '20'

    infile = utils.get_test_data('test-filter-abund-Ns.fq')
    in_dir = os.path.dirname(infile)
    fifo = utils.get_temp_filename('fifo')
    outfile = utils.get_temp_filename('outfile')

    # Use a fifo to copy stdout to a file for checking
    os.mkfifo(fifo)
    thread = threading.Thread(target=write_by_chunks, args=(infile, fifo))
    thread.start()

    # Execute diginorm
    script = 'normalize-by-median.py'
    args = ['-C', CUTOFF, '-k', '17', '-o', outfile, fifo]
    (status, out, err) = utils.runscript(script, args, in_dir)

    # Merge the thread
    thread.join()

    assert os.path.exists(outfile), outfile
    assert 'Total number of unique k-mers: 98' in err, err
    assert 'fifo is empty' not in err, err


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
    PAIRING = ['--force_single']

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

    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)

    shutil.copyfile(utils.get_test_data('dn-test-some-paired-all-keep.fa'),
                    infile)

    script = 'normalize-by-median.py'
    args = list(CUTOFF) + ['-k', '15', infile]

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


def test_diginorm_basic_functionality_5():
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


def test_normalize_by_median_outfile_closed_err():
    infile1 = utils.get_test_data('paired-mixed.fa.pe')
    infile2 = utils.get_test_data("test-abund-read-2.fa")
    outfile = utils.get_temp_filename('outfile_xxx')
    script = 'normalize-by-median.py'
    args = ['-o', outfile, infile1, infile2]
    (status, out, err) = utils.runscript(script, args)
    assert status == 0, (out, err)
    assert os.path.exists(outfile)
