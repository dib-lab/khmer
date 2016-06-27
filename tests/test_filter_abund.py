# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2016, The Regents of the University of California.
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
# pylint: disable=missing-docstring
import os
import shutil
import khmer
import screed
from . import khmer_tst_utils as utils
from .test_scripts import _make_counting


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


def test_filter_abund_2_stdin_gzip_out():
    infile = utils.get_temp_filename('test.fa')
    in_dir = os.path.dirname(infile)
    outfile = utils.get_temp_filename('out.fa.gz')

    shutil.copyfile(utils.get_test_data('test-abund-read-2.fa'), infile)
    counting_ht = _make_counting(infile, K=17)

    script = 'filter-abund.py'
    args = ['-C', '1', counting_ht, infile, '-o', outfile, '--gzip']
    (status, out, err) = utils.runscript(script, args, in_dir, fail_ok=True)
    print(out)
    print(err)
    assert status == 0

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


def test_outfile():
    infile = utils.get_test_data('paired-mixed-witherror.fa.pe')
    outfile = utils.get_temp_filename('paired-mixed-witherror.fa.pe.abundfilt')
    script = 'filter-abund-single.py'
    args = ['-o', outfile, infile]
    (status, out, err) = utils.runscript(script, args)
    md5hash = utils._calc_md5(open(outfile, 'rb'))
    assert md5hash == 'f17122f4c0c3dc0bcc4eeb375de93040', md5hash
