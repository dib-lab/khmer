#
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) Michigan State University, 2009-2015. It is licensed under
# the three-clause BSD license; see LICENSE. Contact: ctb@msu.edu
#
import khmer
import khmer_tst_utils as utils
# from nose.tools import assert_almost_equals


def pretty_compare(a, b):
    print len(a), len(b)

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
        print "".join(line1[i:i+60])
        print "".join(line3[i:i+60])
        print "".join(line2[i:i+60])


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


def test_align_nothing():
    ch = khmer.new_counting_hash(10, 1048576, 1)
    read = "ACCAAGGCTCGAGATTTACC"

    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    score, graphAlign, readAlign, trunc = aligner.align(read)

    print score, graphAlign, readAlign

    assert trunc
    assert len(graphAlign) == 0
    assert len(readAlign) == 0


def test_alignnocov():
    ch = khmer.new_counting_hash(10, 1048576, 1)
    read = "ACCTAGGTTCGACATGTACC"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    ch.consume("ACCTAGGTTCGACATGTACC")
    score, graphAlign, readAlign, trunc = aligner.align(read)

    # should be the same
    eq_(readAlign, 'ACCTAGGTTCGACATGTACC')
    eq_(graphAlign, 'ACCTAGGTTCGACATGTACC')
    assert not trunc


def test_align_middle():
    ch = khmer.new_counting_hash(10, 1048576, 1)
    read = "TCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    ch.consume(read)
    score, graphAlign, readAlign, trunc = aligner.align(read)

    # should be the same
    eq_(readAlign, read)
    eq_(graphAlign, read)
    assert not trunc


def test_align_middle_trunc():
    ch = khmer.new_counting_hash(10, 1048576, 1)
    read = "TCGACAAGTCCTTGACAGATGGGGGG"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")

    # omit suffix from graph
    ch.consume(read[:-5])
    score, graphAlign, readAlign, trunc = aligner.align(read)

    # should not be the same...
    neq_(readAlign, read)
    neq_(graphAlign, read)

    eq_(readAlign, read[:-5])
    eq_(graphAlign, read[:-5])

    # ...but truncated
    assert trunc


def test_align_middle_trunc_2():
    ch = khmer.new_counting_hash(10, 1048576, 1)
    read = "GGGGGGGGGGGGTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AAAAAAAAAAAATCGACAAGTCCTTGACAGAT")

    # omit prefix from graph
    ch.consume(read[12:])
    score, graphAlign, readAlign, trunc = aligner.align(read)

    # here, the alignment must start not at the beginning
    print readAlign
    print graphAlign

    eq_(readAlign, read[12:])
    eq_(graphAlign, read[12:])

    # ...but truncated
    assert trunc


def test_align_fwd_nothing():
    ch = khmer.new_counting_hash(10, 1048576, 1)
    read = "ACCAAGGCTCGAGATTTACC"

    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    score, graphAlign, readAlign, trunc, _ = aligner.align_forward(read)

    print score, graphAlign, readAlign

    assert trunc
    assert len(graphAlign) == 0
    assert len(readAlign) == 0


def test_align_fwd_nocov():
    ch = khmer.new_counting_hash(10, 1048576, 1)
    read = "ACCTAGGTTCGACATGTACC"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    ch.consume("ACCTAGGTTCGACATGTACC")
    score, graphAlign, readAlign, trunc, _ = aligner.align_forward(read)

    # should be the same
    eq_(readAlign, 'ACCTAGGTTCGACATGTACC')
    eq_(graphAlign, 'ACCTAGGTTCGACATGTACC')
    assert not trunc


def test_align_fwd_middle():
    ch = khmer.new_counting_hash(10, 1048576, 1)
    read = "TCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    ch.consume(read)
    score, graphAlign, readAlign, trunc, _ = aligner.align_forward(read)

    # should be the same
    eq_(readAlign, read)
    eq_(graphAlign, read)
    assert not trunc


def test_align_fwd_middle_trunc():
    ch = khmer.new_counting_hash(10, 1048576, 1)
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
    ch = khmer.new_counting_hash(10, 1048576, 1)
    read = "GGGGGGGGGGGGTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(20):
        ch.consume("AAAAAAAAAAAATCGACAAGTCCTTGACAGAT")

    # omit prefix from graph
    ch.consume(read[12:])
    score, graphAlign, readAlign, trunc, _ = aligner.align_forward(read)

    # this will fail, because align_forward chooses the first kmer as the
    # seed.
    assert not readAlign
    assert not graphAlign
    assert trunc


def test_align_fwd_covs_1():
    K = 10

    ch = khmer.new_counting_hash(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("CTCGACAAGTCCTTGACAGAT")
    #           ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    for start in range(0, len(read) - K + 1):
        print ch.get(read[start:start+K]),
    print ''

    assert len(covs) == len(read)
    assert covs[0] == 19
    assert min(covs[1:-K]) == 20, covs
    assert max(covs) == 20, covs


def test_align_fwd_covs_2():
    K = 10

    ch = khmer.new_counting_hash(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("GACGACAAGTCCTTGACAGAT")
    #            ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    print covs, g
    for start in range(0, len(read) - K + 1):
        print ch.get(read[start:start+K]),
    print ''

    assert len(covs) == len(read)
    assert covs[0] == 19
    assert covs[1] == 19
    assert min(covs[2:-K]) == 20, covs
    assert max(covs) == 20, covs


def test_align_fwd_covs_3():
    K = 10

    ch = khmer.new_counting_hash(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("GTAGACAAGTCCTTGACAGAT")
    #             ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    print covs, g
    for start in range(0, len(read) - K + 1):
        print ch.get(read[start:start+K]),
    print ''

    assert len(covs) == len(read)
    assert covs[0] == 19
    assert covs[1] == 19
    assert covs[2] == 19
    assert min(covs[3:-K]) == 20, covs
    assert max(covs) == 20, covs


def test_align_fwd_covs_4():
    K = 10

    ch = khmer.new_counting_hash(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("GTCGACAAGTCCTTGACAGAG")
    #                               ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    print covs, g
    for start in range(0, len(read) - K + 1):
        print ch.get(read[start:start+K]),
    print ''

    assert len(covs) == len(read)
    assert covs[-K] == 19
    assert min(covs[:-K]) == 20, covs
    assert max(covs) == 20, covs


def test_align_fwd_covs_5():
    K = 10

    ch = khmer.new_counting_hash(K, 1048576, 1)
    read = "GTCGACAAGTCCTTGACAGAT"
    aligner = khmer.ReadAligner(ch, 0, 0)
    for i in range(19):
        ch.consume(read)

    ch.consume("GTCGACAAGTCCTTGACAGCT")
    #                              ^
    score, g, r, is_t, covs = aligner.align_forward(read)

    print covs, g
    for start in range(0, len(read) - K + 1):
        print ch.get(read[start:start+K]),
    print ''

    assert len(covs) == len(read)
    assert covs[-K] == 19
    assert covs[-K - 1] == 19
    assert min(covs[:-K - 1]) == 20, covs
    assert max(covs) == 20, covs


def test_simple_readalign():
    return  # DISABLED @CTB
    ch = khmer.new_counting_hash(10, 1048576, 1)
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
    eq_(graphAlign, 'AGCTAGGTTCGACATGTCC-')
    eq_(readAlign, 'ACCTAGGTTCGACAAGTACc')


def test_readalign():
    return  # DISABLED!
    ch = khmer.new_counting_hash(10, 1048576, 1)
    aligner = khmer.ReadAligner(ch, 1, 0)
    for i in range(20):
        ch.consume("AGAGGGAAAGCTAGGTTCGACAAGTCCTTGACAGAT")
    read = "ACCTAGGTTCGACATGTACC"
    #                      ^^            ^  ^

    ch.consume("GCTTTTAAAAAGGTTCGACAAAGGCCCGGG")

    score, graphAlign, readAlign, trunc = aligner.align(read)

    eq_(readAlign, 'ACCTAGGTTCGACATGTACc')
    eq_(graphAlign, 'AGCTAGGTTCGACAAGTCC-')


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
        "CTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTTAACAA"
        "CCTCTTTAC",
        "score": 278.376028204,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCG"
        "CTTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCT"
        "TAACAACCTCTTTAC",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
        "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTT"
        "AACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAA"
        "CTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATCTTAACAA"
        "CCTCTTTAC",
        "score": 271.753976385,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCG"
        "CTTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCT"
        "TAACAACCTCTTTAC",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
        "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTATTGCAATCTT"
        "AACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAAC"
        "TGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTTAACAAC"
        "CTCTTTAC",
        "score": 276.416710585,
        "graph_aln": "TAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
        "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTT"
        "AACAACCTCTTTAC",
        "read_aln": "TAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCT"
        "TTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTTA"
        "ACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TAAATGCGCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAAC"
        "TGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTTAACAAC"
        "CTCTTTAC",
        "score": 269.794658765,
        "graph_aln": "TAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
        "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTT"
        "AACAACCTCTTTAC",
        "read_aln": "TAAATGCGCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCT"
        "TTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAATCTTA"
        "ACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAA",
        "score": 97.5386525659,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTAGATGTTTGATTATCAA",
        "score": 90.9166007464,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTAGATGTTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTATTGATTATCAA",
        "score": 92.9385894977,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGT-TTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTaTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATTGTTTGATTATCAA",
        "score": 84.3383420486,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATAtGTTTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATT-GTTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTATTGATTATCAA",
        "score": 92.9385894977,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGT-TTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTaTTGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTATAGATTATCAA",
        "score": 86.3165376783,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGT-TTGATTATCAA",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTaTAGATTATCAA",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATAATTTTGCCGCTTTAAC"
        "TGGGTCTAGTTTCTACTGCAAACTTTCCACCAACTAGTTTTTCTGCATCCTTTGTTGCAATCTTAACAA"
        "CCTCTTTAC",
        "score": 236.115256507,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAaTT-TtGCC"
        "GCTTTAACTGGGTCT-GTTTCTACTGCAAACTTTCCACCAACAAGTTTTTCTGCATCCTGTGTTGCAAT"
        "CTTAACAACCTCTTTAC",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATAA-TTtT-GCCG"
        "CTTTAACTGGGTCTaGTTTCTACTGCAAACTTTCCACCAACTAGTTTTTCTGCATCCTTTGTTGCAATC"
        "TTAACAACCTCTTTAC",
        "truncated": False
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGAAAATAATTAAAAAAAAAAAAA"
        "AAAAAAAAAAAAAAAAAAAAAAAAAA",
        "score": 44.7543247314,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATatgtt",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTAT-----",
        "truncated": True
    },
    {
        "seq": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGCTTTAA"
        "CTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGAAAAATGTCATCCTGTATTGCAATCTTAACAA"
        "CCTCTTTAC",
        "score": 227.446444943,
        "graph_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCG"
        "CTTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGTtTTTCTG-CATCCTGTGTTGCAATC"
        "TTAACAACCTCTTTAC",
        "read_aln": "TTAAATGCCCAATTTTTCCCTCTTTTCTTCTATATGTTTGATTATCAATTTTGCCGC"
        "TTTAACTGGGTCTGTTTCTACTGCAAACTTTCCACCAACAAGA-AAAATGtCATCCTGTATTGCAATCT"
        "TAACAACCTCTTTAC",
        "truncated": False
    }
]


def test_readalign_new():
    return  # DISABLED
    ch = khmer.new_counting_hash(32, 1048576, 1)
    aligner = khmer.ReadAligner(ch, 1, 0)
    for seq in ht_seqs:
        ch.consume(seq)

    for query in queries:
        score, graphAlign, readAlign, trunc = aligner.align(query["seq"])
        print graphAlign
        print readAlign
        eq_(graphAlign, query["graph_aln"])
        eq_(readAlign, query["read_aln"])
        assert trunc == query["truncated"]
        # assert_almost_equals(score, query["score"])


def test_readaligner_load():
    ct = khmer.new_counting_hash(32, 1048576, 1)
    parameters_json = utils.get_test_data('readaligner-default.json')
    aligner = khmer.ReadAligner(ct, 0, 0, filename=parameters_json)
    assert aligner.get_scoring_matrix()[0] == -0.06642736173897607, (
        aligner.get_scoring_matrix()[0])
    assert aligner.get_transition_probabilities()[0][0] == \
        -0.021973842014145723, aligner.get_transition_probabilities()[0][0]
    for seq in ht_seqs:
        ct.consume(seq)

    for query in queries:
        aligner.align(query['seq'])
