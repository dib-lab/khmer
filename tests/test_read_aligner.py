# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2013-2015, Michigan State University.
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
# pylint: disable=missing-docstring,no-member,invalid-name,unused-variable

import pytest
import khmer
from . import khmer_tst_utils as utils


def pretty_compare(a, b):
    print(len(a), len(b))

    line1 = []
    line2 = []
    line3 = []
    for (x, y) in zip(a, b):
        line1.append(x)
        line2.append(y)
        if x == y:
            line3.append('|')
        else:
            line3.append('x')

    for i in range(0, len(line1), 60):
        print("".join(line1[i:i + 60]))
        print("".join(line3[i:i + 60]))
        print("".join(line2[i:i + 60]))


def eq_(v1, v2):
    assert len(v1)
    if v1 != v2:
        pretty_compare(v1, v2)
    assert v1 == v2, (v1, v2)


def neq_(v1, v2):
    assert len(v1)
    if v1 == v2:
        pretty_compare(v1, v2)
    assert v1 != v2, (v1, v2)


def test_graph_attribute():
    ch = khmer.Countgraph(10, 1048576, 1)
    aligner = khmer.ReadAligner(ch, 0, 0)
    assert aligner.graph is ch


def test_scoring_matrix():
    ch = khmer.Countgraph(10, 1048576, 1)
    aligner = khmer.ReadAligner(ch, 0, 0)
    assert aligner.scoring_matrix == aligner.defaultScoringMatrix


def test_transition_probabilities():
    ch = khmer.Countgraph(10, 1048576, 1)
    aligner = khmer.ReadAligner(ch)
    assert aligner.transition_probabilities == \
        aligner.defaultTransitionProbabilities


def test_align_nothing():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "ACCAAGGCTCGAGATTTACC"

    aligner = khmer.ReadAligner(ch, 0, 0)
    for _ in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    score, graphAlign, readAlign, trunc = aligner.align(read)

    print(score, graphAlign, readAlign)

    assert trunc
    assert len(graphAlign) == 0
    assert len(readAlign) == 0


def test_alignnocov():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "ACCTAGGTTCGACATGTACC"
    aligner = khmer.ReadAligner(ch, trusted_cov_cutoff=0,
                                bits_theta=0)
    for _ in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    ch.consume("ACCTAGGTTCGACATGTACC")
    _, graphAlign, readAlign, trunc = aligner.align(read)

    # should be the same
    eq_(readAlign, 'ACCTAGGTTCGACATGTACC')
    eq_(graphAlign, 'ACCTAGGTTCGACATGTACC')
    assert not trunc


def test_align_middle():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "TCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, trusted_cov_cutoff=0,
                                bits_theta=0)
    for _ in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    ch.consume(read)
    _, graphAlign, readAlign, trunc = aligner.align(read)

    # should be the same
    eq_(readAlign, read)
    eq_(graphAlign, read)
    assert not trunc


@pytest.mark.known_failing
def test_align_middle_trunc():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "TCGACAAGTCCTTGACAGATGGGGGG"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for _ in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")

    # omit suffix from graph
    ch.consume(read[:-5])
    _, graphAlign, readAlign, trunc = aligner.align(read)

    # should not be the same...
    neq_(readAlign, read)
    neq_(graphAlign, read)

    eq_(readAlign, read[:-5])
    eq_(graphAlign, read[:-5])

    # ...but truncated
    assert trunc


@pytest.mark.known_failing
def test_align_middle_trunc_2():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "GGGGGGGGGGGGTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for _ in range(20):
        ch.consume("AAAAAAAAAAAATCGACAAGTCCTTGACAGAT")

    # omit prefix from graph
    ch.consume(read[12:])
    _, graphAlign, readAlign, trunc = aligner.align(read)

    # here, the alignment must start not at the beginning
    print(readAlign)
    print(graphAlign)

    eq_(readAlign, read[12:])
    eq_(graphAlign, read[12:])

    # ...but truncated
    assert trunc


def test_align_fwd_nothing():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "ACCAAGGCTCGAGATTTACC"

    aligner = khmer.ReadAligner(ch, 0, 0)
    for _ in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    score, graphAlign, readAlign, trunc, _ = aligner.align_forward(read)

    print(score, graphAlign, readAlign)

    assert trunc
    assert len(graphAlign) == 0
    assert len(readAlign) == 0


def test_align_fwd_nocov():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "ACCTAGGTTCGACATGTACC"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for _ in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    ch.consume("ACCTAGGTTCGACATGTACC")
    _, graphAlign, readAlign, trunc, _ = aligner.align_forward(read)

    # should be the same
    eq_(readAlign, 'ACCTAGGTTCGACATGTACC')
    eq_(graphAlign, 'ACCTAGGTTCGACATGTACC')
    assert not trunc


def test_align_fwd_middle():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "TCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for _ in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    ch.consume(read)
    score, graphAlign, readAlign, trunc, _ = aligner.align_forward(read)

    # should be the same
    eq_(readAlign, read)
    eq_(graphAlign, read)
    assert not trunc


@pytest.mark.known_failing
def test_align_fwd_middle_trunc():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "TCGACAAGTCCTTGACAGATGGGGGG"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")

    # omit suffix from graph
    ch.consume(read[:-5])
    score, graphAlign, readAlign, trunc, _ = aligner.align_forward(read)

    # should not be the same...
    neq_(readAlign, read)
    neq_(graphAlign, read)

    eq_(readAlign, read[:-5])
    eq_(graphAlign, read[:-5])

    # ...but truncated
    assert trunc


def test_align_fwd_middle_trunc_2():
    ch = khmer.Countgraph(10, 1048576, 1)
    read = "GGGGGGGGGGGGTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for _ in range(20):
        ch.consume("AAAAAAAAAAAATCGACAAGTCCTTGACAGAT")

    # omit prefix from graph
    ch.consume(read[12:])
    _, graphAlign, readAlign, trunc, _ = aligner.align_forward(read)

    # this will fail, because align_forward chooses the first kmer as the
    # seed.
    assert not readAlign
    assert not graphAlign
    assert trunc


def test_align_fwd_covs_1():
    K = 10

    ch = khmer.Countgraph(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("CTCGACAAGTCCTTGACAGAT")
    #           ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    for start in range(0, len(read) - K + 1):
        print(ch.get(read[start:start + K]), end=' ')
    print('')

    assert len(covs) == len(read)
    assert covs[0] == 19
    assert min(covs[1:-K]) == 20, covs
    assert max(covs) == 20, covs


def test_align_fwd_covs_2():
    K = 10

    ch = khmer.Countgraph(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("GACGACAAGTCCTTGACAGAT")
    #            ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    print(covs, g)
    for start in range(0, len(read) - K + 1):
        print(ch.get(read[start:start + K]), end=' ')
    print('')

    assert len(covs) == len(read)
    assert covs[0] == 19
    assert covs[1] == 19
    assert min(covs[2:-K]) == 20, covs
    assert max(covs) == 20, covs


def test_align_fwd_covs_3():
    K = 10

    ch = khmer.Countgraph(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("GTAGACAAGTCCTTGACAGAT")
    #             ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    print(covs, g)
    for start in range(0, len(read) - K + 1):
        print(ch.get(read[start:start + K]), end=' ')
    print('')

    assert len(covs) == len(read)
    assert covs[0] == 19
    assert covs[1] == 19
    assert covs[2] == 19
    assert min(covs[3:-K]) == 20, covs
    assert max(covs) == 20, covs


def test_align_fwd_covs_4():
    K = 10

    ch = khmer.Countgraph(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("GTCGACAAGTCCTTGACAGAG")
    #                               ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    print(covs, g)
    for start in range(0, len(read) - K + 1):
        print(ch.get(read[start:start + K]), end=' ')
    print('')

    assert len(covs) == len(read)
    assert covs[-K] == 19
    assert min(covs[:-K]) == 20, covs
    assert max(covs) == 20, covs


def test_align_fwd_covs_5():
    K = 10

    ch = khmer.Countgraph(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("GTCGACAAGTCCTTGACAGCT")
    #                              ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    print(covs, g)
    for start in range(0, len(read) - K + 1):
        print(ch.get(read[start:start + K]), end=' ')
    print('')

    assert len(covs) == len(read)
    assert covs[-K] == 19
    assert covs[-K - 1] == 19
    assert min(covs[:-K - 1]) == 20, covs
    assert max(covs) == 20, covs


@pytest.mark.known_failing
def test_simple_readalign():
    ch = khmer.Countgraph(10, 1048576, 1)
    aligner = khmer.ReadAligner(ch, 2, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACATGTCCTTGACAGAT")
    read = "ACCTAGGTTCGACAAGTACC"
    #                      ^^            ^  ^
    ch.consume("GCTTTTAAAAAGGTTCGACAAAGGCCCGGG")
    # CCCGGGCCTTTGTCGAACCTTTTTAAAAGC

    score, graphAlign, readAlign, trunc = aligner.align(read)

#                        AGCTAGGTTCGACAAGT CCT
#                        ACCTAGGTTCGACAAGTaCC
#                        --CTAGGTTCGACATGT-CC
    eq_(graphAlign, 'AGCTAGGTTCGACATGTCCT')
    eq_(readAlign, 'ACCTAGGTTCGACAAGTACC')


@pytest.mark.known_failing
def test_readalign():
    ch = khmer.Countgraph(10, 1048576, 1)
    aligner = khmer.ReadAligner(ch, 1, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    read = "ACCTAGGTTCGACATGTACC"
    #                      ^^            ^  ^

    ch.consume("GCTTTTAAAAAGGTTCGACAAAGGCCCGGG")

    score, graphAlign, readAlign, _ = aligner.align(read)

    eq_(readAlign, 'ACCTAGGTTCGACATGTACC')
    eq_(graphAlign, 'AGCTAGGTTCGACAAGTCCT')


ht_seqs = ["TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAACTGG"
           "GTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTTAACAAC"
           "CTCTTTAC",
           "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAACTGG"
           "GTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTTAACAAC"
           "CTCTTTAC",
           "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAACTGG"
           "GTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATCTTAACAAC"
           "CTCTTTAC",
           "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAACTGG"
           "GTCTGTTTCTACTGCAAACTTTCCACCAACAAGAAAAATGTCATCCTGTATTGCAATCTTAACAAC"
           "CTCTTTAC"]

queries = [
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAA"
               "CTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATC"
               "TTAACAACCTCTTTAC",
        "score": 274.76338282696173,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCG"
                     "CTTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCC"
                     "TGTGTTGCAATCTTAACAACCTCTTTAC",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
                    "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTG"
                    "TGTTGCAATCTTAACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAA"
               "CTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATC"
               "TAACAACCTCTTTAC",
        "score": 274.76338282696173,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCG"
                     "CTTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCC"
                     "TGTATTGCAATCTTAACAACCTCTTTAC",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
                    "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTG"
                    "TATTGCAATCTTAACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAAC"
               "TGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCT"
               "TAACAACCTCTTTAC",
        "score": 272.841515695261,
        "graph_aln": "TAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
                     "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCT"
                     "GTGTTGCAATCTTAACAACCTCTTTAC",
        "read_aln": "TAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCT"
                    "TTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGT"
                    "GTTGCAATCTTAACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TAAATGCGCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAAC"
               "TGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCT"
               "TAACAACCTCTTTAC",
        "score": 268.2640868672253,
        "graph_aln": "TAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
                     "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCT"
                     "GTGTTGCAATCTTAACAACCTCTTTAC",
        "read_aln": "TAAATGCGCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCT"
                    "TTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGT"
                    "GTTGCAATCTTAACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAA",
        "score": 97.37145206396536,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTAGATGTTTGATTATCAA",
        "score": 92.79402323592961,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTAGATGTTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTATTGATTATCAA",
        "score": 84.74620322710143,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGT-TTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTaTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATTGTTTGATTATCAA",
        "score": 82.2182409986759,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATaTGTTTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTAT-TGTTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTATTGATTATCAA",
        "score": 84.74620322710143,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGT-TTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTaTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTATAGATTATCAA",
        "score": 80.1687743990657,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGT-TTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTaTAGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATAATTTTGCCGCTTTAAC"
               "TGGGTCTAGTTTCTACTGCAAACTTTCCACCAACTAGTTTTTCTGCATCCTTTGTTGCAATC"
               "TTAACAACCTCTTTAC",
        "score": 237.81111469018322,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATcAATTTTGCCG"
                     "CTTTAACTGGGTCT-GTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATC"
                     "CTGTGTTGCAATCTTAACAACCTCTTTAC",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTAT-AATTTTGCCGC"
                    "TTTAACTGGGTCTaGTTTCTACTGCAAACTTTCCACCAACTAGTTTTTCTGCATCCT"
                    "TTGTTGCAATCTTAACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGAAAATAATTAAAAAAAAAAAAA"
               "AAAAAAAAAAAAAAAAAAAAAAAAAA",
        "score": 5.331560863368736,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCG"
                     "CTTTAACTGGGTCTGTTTCTACTGCAAACTTT",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGAAAATAATTAAAAAAAA"
                    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAA"
               "CTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGAAAAATGTCATCCTGTATTGCAATC"
               "TTAACAACCTCTTTAC",
        "score": 274.76338282696173,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCG"
                     "CTTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGAAAAATGTCATCC"
                     "TGTATTGCAATCTTAACAACCTCTTTAC",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
                    "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGAAAAATGTCATCCTG"
                    "TATTGCAATCTTAACAACCTCTTTAC",
        "truncated": False
    },
    {  # the motif of 32 bases are identical match to HT seqs, the rest are
        # random. TTAAATGCCCAATTTTTCCCTCTTTTCTTCTAT" is the from HT seqs
        "seq":
        "ACAAGGCCATTTGTTCGCATTCTGAAGCCGGCTTCCACCATGGTACTGGGAAACTGTCGGAATATTAAA"
        "TGCCCAATTTTTCCCTCTTTTCTTCTATCCGCAGTATGGACACTGTTTTCCTGAATTTCATTGACAGTT"
        "TAATTTACTGCGGTCACGCGGAACT",
        "score": 68.17022311739733,
        "graph_aln":
        "ACAAGGCCATTTGTTCGCATTCTGAAGCCGGCTTCCACCATGGTACTGGGAAACTGTCGGAATATTAAA"
        "TGCCCAATTTTTCCCTCTTTTCTTCTATCCGCAGTATGGACACTGTTTTCCTGAATTTCATTGACAGTT"
        "TAATTTACTGCGGTCACGCGGAACT",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTAT",
        "truncated": True,
        "description": "truncated-alignment-bc-missing-kmers"
    },
    {   # Testing for min distance between correctable SNPs
        # 1st SNP is at position 2+K from beginning, 2nd SNP at position 2+K+K
        "seq":
        "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATACGTTTGATTATCAATTTTGCCGCTTTAACTGG"
        "ATCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATCTTAACAAC"
        "CTCTTTAC",
        "score": 265.608525171,
        "graph_aln":
        "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAACTGG"
        "GTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATCTTAACAAC"
        "CTCTTTAC",
        "read_aln":
        "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATACGTTTGATTATCAATTTTGCCGCTTTAACTGG"
        "ATCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATCTTAACAAC"
        "CTCTTTAC",
        "truncated": False,
        "description": "2 SNPs, one K apart",
    },
    {   # Testing for min distance between correctable SNPs
        # 1st SNP is at position 2+K from beginning, 2nd SNP at position
        # 2+K+K-1
        "seq":
        "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATACCTTTGATTATCAATTTTGCCGCTTTAACTGG"
        "GTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATCTTAACAAC"
        "CTCTTTAC",
        "score": 265.608525171,
        "graph_aln":
        "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAACTGG"
        "GTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATCTTAACAAC"
        "CTCTTTAC",
        "read_aln":
        "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATACGTTTGATTATCAATTTTGCCGCTTTAACTAG"
        "GTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATCTTAACAAC"
        "CTCTTTAC",
        "truncated": False,
        "description": "2 SNPs, K-2 apart",
    }

]


def check_query(aligner, query):
    score, graphAlign, readAlign, trunc = aligner.align(query["seq"])
    print(query["seq"])
    print(graphAlign, query["graph_aln"])
    print(readAlign, query["read_aln"])
    print(trunc, query["truncated"])
    print(score, query["score"])
    assert graphAlign == query["graph_aln"], "\n%r != \n%r" % \
        (graphAlign, query["graph_aln"])
    assert readAlign == query["read_aln"], "\n%r != \n%r" % \
        (readAlign, query["read_aln"])
    eq_(trunc, query["truncated"])
    if query["score"] > 0:
        assert round(score - query["score"], 7) == 0


@pytest.mark.known_failing
@pytest.mark.parametrize('query', queries)
def test_readalign_new(query):
    ch = khmer.Countgraph(32, 1048576, 1)
    aligner = khmer.ReadAligner(ch, 1, 0)
    for seq in ht_seqs:
        ch.consume(seq)

    check_query(aligner, query)


def test_readaligner_load():
    ct = khmer.Countgraph(32, 1048576, 1)
    parameters_json = utils.get_test_data('readaligner-default.json')
    a_aligner = khmer.ReadAligner(ct, 0, 0, filename=parameters_json)
    a_scoring_matrix = a_aligner.scoring_matrix
    a_transition_probabilities = a_aligner.transition_probabilities
    assert a_scoring_matrix[0] == -0.06642736173897607, a_scoring_matrix[0]
    assert a_transition_probabilities[0][0] == -0.021973842014145723, (
        a_transition_probabilities[0][0])

    for seq in ht_seqs:
        ct.consume(seq)

    for query in queries:
        a_aligner.align(query['seq'])

    b_aligner = khmer.ReadAligner(
        ct, 0, 0, transition_probabilities=a_transition_probabilities,
        scoring_matrix=a_scoring_matrix)
    b_scoring_matrix = b_aligner.scoring_matrix
    b_transition_probabilities = b_aligner.transition_probabilities
    assert b_scoring_matrix == a_scoring_matrix, (
        a_scoring_matrix, b_scoring_matrix)
    assert b_transition_probabilities == a_transition_probabilities, (
        a_transition_probabilities, b_transition_probabilities)
