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


import gzip
import os
import pytest

from . import khmer_tst_utils as utils
import khmer
import khmer.kfile
import khmer.utils
import screed


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


def test_interleave_reads_no_reformat():
    infile1 = utils.get_test_data('paired.fq.1')
    infile2 = utils.get_test_data('paired.malformat.fq.2')

    ex_outfile = utils.get_test_data('paired.malformat.fq')
    outfile = utils.get_temp_filename('out.fq')

    script = 'interleave-reads.py'
    args = [infile1, infile2, '--no-reformat', '-o', outfile]

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
    infile = utils.copy_test_data('paired-mixed.fq')
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
    infile = utils.copy_test_data('paired-mixed-2.fq')
    in_dir = os.path.dirname(infile)

    script = 'split-paired-reads.py'
    args = ['-0', '/dev/null', infile]

    status, out, err = utils.runscript(script, args, in_dir)
    assert status == 0
    assert "split 6 sequences (3 left, 3 right, 5 orphans)" in err, err


def test_split_paired_reads_2_mixed_fq_orphans_to_file():
    # test input file
    infile = utils.copy_test_data('paired-mixed-2.fq')
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
    infile = utils.copy_test_data('paired-mixed-2.fq')
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
    infile = utils.copy_test_data('paired-mixed-broken.fq')
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
    testdir = utils.get_temp_filename('test')
    output_dir = os.path.join(os.path.dirname(testdir), "out")
    outfile1 = utils.get_temp_filename('paired.fq.1', output_dir)
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


def test_extract_paired_reads_unpaired():
    # test input file
    infile = utils.get_test_data('random-20-a.fa')

    # actual output files...
    outfile1 = utils.get_temp_filename('unpaired.pe.fa')
    in_dir = os.path.dirname(outfile1)
    outfile2 = utils.get_temp_filename('unpaired.se.fa', in_dir)

    script = 'extract-paired-reads.py'
    args = [infile]

    (_, _, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    assert 'no paired reads!? check file formats...' in err, err


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


def test_read_bundler():
    infile = utils.get_test_data('unclean-reads.fastq')
    records = [r for r in khmer.ReadParser(infile)]
    bundle = khmer.utils.ReadBundle(*records)

    raw_seqs = (
        'GGTTGACGGGGNNNAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGCAGCTGTCGTCAGGGGATTTCCG'
        'GGGCGGAGGCCGCAGACGCGAGTGGTGGAGG',
        'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGCAGCTGTCGTCAGGGGANNNCCG'
        'GGGCGGAGGCCGCAGACGCGAGTGGTGGAGG',
    )

    cleaned_seqs = (
        'GGTTGACGGGGAAAAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGCAGCTGTCGTCAGGGGATTTCCG'
        'GGGCGGAGGCCGCAGACGCGAGTGGTGGAGG',
        'GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGAGACAGCAGCCGCAGCTGTCGTCAGGGGAAAACCG'
        'GGGCGGAGGCCGCAGACGCGAGTGGTGGAGG',
    )

    assert bundle.num_reads == 2
    assert bundle.total_length == 200

    for read, raw_seq, clean_seq in zip(bundle.reads, raw_seqs, cleaned_seqs):
        assert read.sequence == raw_seq
        assert read.cleaned_seq == clean_seq


def test_read_bundler_single_read():
    infile = utils.get_test_data('single-read.fq')
    records = [r for r in khmer.ReadParser(infile)]
    bundle = khmer.utils.ReadBundle(*records)
    assert bundle.num_reads == 1
    assert bundle.reads[0].sequence == bundle.reads[0].cleaned_seq


def test_read_bundler_empty_file():
    infile = utils.get_test_data('empty-file')
    with pytest.raises(OSError):
        records = [r for r in khmer.ReadParser(infile)]


def test_read_bundler_empty_list():
    bundle = khmer.utils.ReadBundle(*[])
    assert bundle.num_reads == 0
