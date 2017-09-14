# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
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
# pylint: disable=missing-docstring,invalid-name

# important note -- these tests do not contribute to code coverage, because
# of the use of subprocess to execute.  Most script tests should go into
# test_scripts.py for this reason.


import khmer
from khmer import Nodegraph, Countgraph
import screed
from . import khmer_tst_utils as utils
from .khmer_tst_utils import scriptpath, run_shell_cmd
from .test_scripts import _make_counting
import os.path
import difflib


def files_are_equal(a, b):
    al = open(a).readlines()
    bl = open(b).readlines()

    return al == bl


def diff_files(a, b):
    al = open(a).readlines()
    bl = open(b).readlines()

    results = "\n".join(difflib.context_diff(al, bl, fromfile=a, tofile=b))
    return results


def test_interleave_split_1():
    in1 = utils.get_test_data('paired.fq.1')
    in2 = utils.get_test_data('paired.fq.2')

    out1 = utils.get_temp_filename('a.fa')
    out2 = utils.get_temp_filename('b.fa')

    cmd = """
       {scripts}/interleave-reads.py {in1} {in2} -o -             |
       {scripts}/split-paired-reads.py -1 {out1} -2 {out2} -
    """

    cmd = cmd.format(scripts=scriptpath(),
                     in1=in1, in2=in2,
                     out1=out1, out2=out2)

    run_shell_cmd(cmd)

    assert files_are_equal(in1, out1), diff_files(in1, out1)
    assert files_are_equal(in2, out2), diff_files(in2, out2)


def test_interleave_split_2_fail():
    in1 = utils.get_test_data('paired.fq.1')
    in2 = utils.get_test_data('paired.fq.2')

    out1 = utils.get_temp_filename('a.fa')
    out2 = utils.get_temp_filename('b.fa')

    cmd = """
       {scripts}/interleave-reads.py {in1} {in2} -o -             |
       {scripts}/split-paired-reads.py -
    """

    cmd = cmd.format(scripts=scriptpath(),
                     in1=in1, in2=in2,
                     out1=out1, out2=out2)

    (status, _, err) = run_shell_cmd(cmd, fail_ok=True)
    assert status != 0
    assert "Accepting input from stdin; output filenames must be provided." \
           in err, err


def test_interleave_split_3_out1():
    in1 = utils.get_test_data('paired.fq.1')
    in2 = utils.get_test_data('paired.fq.2')

    out1 = utils.get_temp_filename('a.fa')
    out2 = utils.get_temp_filename('b.fa')

    cmd = """
       {scripts}/interleave-reads.py {in1} {in2} -o -             |
       {scripts}/split-paired-reads.py -1 - -2 {out2} - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(),
                     in1=in1, in2=in2,
                     out1=out1, out2=out2)

    run_shell_cmd(cmd)

    assert files_are_equal(in1, out1), diff_files(in1, out1)
    assert files_are_equal(in2, out2), diff_files(in2, out2)


def test_interleave_split_3_out2():
    in1 = utils.get_test_data('paired.fq.1')
    in2 = utils.get_test_data('paired.fq.2')

    out1 = utils.get_temp_filename('a.fa')
    out2 = utils.get_temp_filename('b.fa')

    cmd = """
       {scripts}/interleave-reads.py {in1} {in2} -o -             |
       {scripts}/split-paired-reads.py -1 {out1} -2 - - > {out2}
    """

    cmd = cmd.format(scripts=scriptpath(),
                     in1=in1, in2=in2,
                     out1=out1, out2=out2)

    run_shell_cmd(cmd)

    assert files_are_equal(in1, out1), diff_files(in1, out1)
    assert files_are_equal(in2, out2), diff_files(in2, out2)


def test_interleave_split_3_out0():
    in1 = utils.get_test_data('paired-mixed-broken.fq')

    out1 = utils.get_temp_filename('a.fa')
    out2 = utils.get_temp_filename('b.fa')
    out3 = utils.get_temp_filename('c.fa')

    cmd = """
       cat {in1} |
       {scripts}/split-paired-reads.py -1 {out1} -2 {out2} -0 - - > {out3}
    """

    cmd = cmd.format(scripts=scriptpath(),
                     in1=in1,
                     out1=out1, out2=out2, out3=out3)

    run_shell_cmd(cmd)

    assert files_are_equal(in1, out3), diff_files(in1, out3)
    assert len(open(out1, 'rb').read()) == 0
    assert len(open(out2, 'rb').read()) == 0


def test_extract_paired_pe():
    in1 = utils.get_test_data('paired-mixed.fq')
    out_test = utils.get_test_data('paired-mixed.fq.pe')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/extract-paired-reads.py - -p - -s /dev/null > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert files_are_equal(out1, out_test), diff_files(out1, out_test)


def test_extract_paired_se():
    in1 = utils.get_test_data('paired-mixed.fq')
    out_test = utils.get_test_data('paired-mixed.fq.se')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/extract-paired-reads.py - -p /dev/null -s - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert files_are_equal(out1, out_test), diff_files(out1, out_test)


def test_extract_paired_stdin_equivalence():
    # Use '/dev/stdin' instead of '-' to check it is treated the same way
    in1 = utils.get_test_data('paired-mixed.fq')
    out_test = utils.get_test_data('paired-mixed.fq.se')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/extract-paired-reads.py /dev/stdin -p /dev/null -s - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert files_are_equal(out1, out_test), diff_files(out1, out_test)


def test_extract_paired_se_fail():
    in1 = utils.get_test_data('paired-mixed.fq')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/extract-paired-reads.py -p /dev/null - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    (status, _, err) = run_shell_cmd(cmd, fail_ok=True)
    assert status != 0
    assert "Accepting input from stdin; output filenames must be provided." \
           in err, err


def test_norm_by_median_1():
    in1 = utils.get_test_data('paired-mixed.fq')
    out_test = utils.get_test_data('paired-mixed.fq.pe')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/extract-paired-reads.py - -p - -s /dev/null |
       {scripts}/normalize-by-median.py - -o - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert files_are_equal(out1, out_test), diff_files(out1, out_test)


def test_norm_by_median_2_fail():
    in1 = utils.get_test_data('paired-mixed.fq')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/extract-paired-reads.py - -p - -s /dev/null |
       {scripts}/normalize-by-median.py -p - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    (status, _, err) = run_shell_cmd(cmd, fail_ok=True)
    assert status != 0
    assert "Accepting input from stdin; output filename must be provided with"\
           in err, err


def test_sample_reads_randomly_1():
    in1 = utils.get_test_data('paired-mixed.fq')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/sample-reads-randomly.py - -o - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert files_are_equal(in1, out1), diff_files(in1, out1)


def test_sample_reads_randomly_2_fail():
    in1 = utils.get_test_data('paired-mixed.fq')
    out1 = utils.get_temp_filename('a.fq')

    cmd = """
       cat {in1} |
       {scripts}/sample-reads-randomly.py - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    (status, _, err) = run_shell_cmd(cmd, fail_ok=True)
    assert status != 0
    assert "Accepting input from stdin; output filename must be provided with"\
           in err, err


def test_extract_long_sequences_1():
    in1 = utils.get_test_data('paired-mixed.fa')
    out1 = utils.get_temp_filename('a.fa')

    cmd = """
       cat {in1} |
       {scripts}/extract-long-sequences.py - -l 10 > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    countlines = sum(1 for line in open(out1))
    assert countlines == 22, countlines


def test_fastq_to_fasta_1():
    in1 = utils.get_test_data('test-fastq-reads.fq')
    out1 = utils.get_temp_filename('clean.fa')
    out_test = utils.get_test_data('test-fastq-reads.fa')

    cmd = """
       cat {in1} |
       {scripts}/fastq-to-fasta.py - -o - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)
    assert files_are_equal(out1, out_test), diff_files(out1, out_test)


def test_load_into_counting_1():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.ct')

    cmd = """
       cat {in1} |
       {scripts}/load-into-counting.py -x 1e3 -N 2 -k 20 {out1} - \
       2> /dev/null
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)
    print(cmd)

    run_shell_cmd(cmd)
    assert os.path.exists(out1)
    Countgraph.load(out1)


def test_load_graph_1():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.ct')

    cmd = """
       cat {in1} |
       {scripts}/load-graph.py -x 1e3 -N 2 -k 20 {out1} - \
       2> /dev/null
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)
    print(cmd)

    run_shell_cmd(cmd)
    assert os.path.exists(out1)
    Nodegraph.load(out1)


def test_filter_abund_1():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.abundfilt')

    countgraph = _make_counting(in1, K=17)

    cmd = """
       cat {in1} |
       {scripts}/filter-abund.py {countgraph} - -o - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1,
                     countgraph=countgraph)

    run_shell_cmd(cmd)

    assert os.path.exists(out1)
    seqs = set([r.sequence for r in screed.open(out1)])

    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_filter_abund_2_fail():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.abundfilt')

    countgraph = _make_counting(in1, K=17)

    cmd = """
       cat {in1} |
       {scripts}/filter-abund.py {countgraph} - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1,
                     countgraph=countgraph)

    status, _, err = run_shell_cmd(cmd, fail_ok=True)
    assert status != 0
    assert "Accepting input from stdin; output filename must be provided with"\
           in err, err


def test_abundance_dist_1():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.dist')

    countgraph = _make_counting(in1, K=17)
    assert os.path.exists(countgraph)

    cmd = """
       cat {in1} |
       {scripts}/abundance-dist.py -z {countgraph} - - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1,
                     countgraph=countgraph)

    run_shell_cmd(cmd)

    assert os.path.exists(out1)
    with open(out1) as fpout1:
        line = fpout1.readline().strip()
        line = fpout1.readline().strip()
        assert line == '1,96,96,0.98', line
        line = fpout1.readline().strip()
        assert line == '1001,2,98,1.0', line


def test_trim_low_abund_1():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.abundtrim')

    cmd = """
       cat {in1} |
       {scripts}/trim-low-abund.py -k 17 -x 1e7 -N 2 - -o - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert os.path.exists(out1)
    seqs = set([r.sequence for r in screed.open(out1)])

    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_trim_low_abund_smallcount():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.abundtrim')

    cmd = """
       cat {in1} |
       {scripts}/trim-low-abund.py --small-count \
         -k 17 -x 1e7 -N 2 - -o - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert os.path.exists(out1)
    seqs = set([r.sequence for r in screed.open(out1)])

    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_trim_low_abund_1_gzip_o():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.abundtrim.gz')

    cmd = """
       cat {in1} |
       {scripts}/trim-low-abund.py -k 17 -x 1e7 -N 2 - -o - --gzip > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert os.path.exists(out1)
    seqs = set([r.sequence for r in screed.open(out1)])

    assert len(seqs) == 1, seqs
    assert 'GGTTGACGGGGCTCAGGG' in seqs


def test_trim_low_abund_2_fail():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.abundtrim')

    cmd = """
       cat {in1} |
       {scripts}/trim-low-abund.py -k 17 -x 1e7 -N 2 - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    (status, _, err) = run_shell_cmd(cmd, fail_ok=True)
    assert status != 0
    assert "Accepting input from stdin; output filename must be provided with"\
           in err, err


def test_count_median_1():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.counts')

    countgraph = _make_counting(in1, K=8)
    cmd = """
       cat {in1} |
       {scripts}/count-median.py {countgraph} - - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), countgraph=countgraph,
                     in1=in1, out1=out1)

    run_shell_cmd(cmd)

    assert os.path.exists(out1), out1
    data = [x.strip() for x in open(out1)]
    data = set(data)
    assert len(data) == 3, data
    assert 'seq,1001,1001.0,0.0,18' in data
    assert '895:1:37:17593:9954/1,1,103.803741455,303.702941895,114' in data


def test_readstats_1():
    in1 = utils.get_test_data('test-abund-read-2.fa')
    out1 = utils.get_temp_filename('out.stats')

    cmd = """
       cat {in1} |
       {scripts}/readstats.py --csv - > {out1}
    """

    cmd = cmd.format(scripts=scriptpath(), in1=in1, out1=out1)

    run_shell_cmd(cmd)
    assert '18114,1001,18.1,-' in open(out1).read(), open(out1).read()


def test_unique_kmers_stream_out_fasta():
    infile = utils.get_test_data('random-20-a.fa')

    cmd = "{scripts}/unique-kmers.py -k 20 -e 0.01 --stream-records {infile}"
    cmd = cmd.format(scripts=scriptpath(), infile=infile)

    (_, out, err) = run_shell_cmd(cmd)

    expected = ('Estimated number of unique 20-mers in {infile}: 3950'
                .format(infile=infile))
    assert expected in err
    assert 'Total estimated number of unique 20-mers: 3950' in err

    assert '>45' in out
    assert "ATACGCCACTCGACTTGGCTCGCCCTCGATCTAAAATAGCGGTCGTGTTGGGTTAACAA" in out


def test_unique_kmers_stream_out_fastq_with_N():
    infile = utils.get_test_data('test-filter-abund-Ns.fq')

    cmd = "{scripts}/unique-kmers.py -k 20 -e 0.01 --stream-records {infile}"
    cmd = cmd.format(scripts=scriptpath(), infile=infile)

    (_, out, err) = run_shell_cmd(cmd)

    expected = ('Estimated number of unique 20-mers in {infile}: 94'
                .format(infile=infile))
    assert expected in err
    assert 'Total estimated number of unique 20-mers: 94' in err

    assert '@895:1:37:17593:9954 1::FOO_withN' in out
    assert "GGTTGACGGGGCTCAGGGGGCGGCTGACTCCGAGNGACAGCAGCCGCAGCTGTCGTCA" in out
    assert "##########################################################" in out
