# -*- coding: UTF-8 -*-
# cython: c_string_type=unicode, c_string_encoding=utf8

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from graphs cimport CpHashgraph
from aligner cimport CpReadAligner
from oxli_types cimport *
from math import log2
import json

cdef class ReadAligner:

    default_trans_probs = (
    #    _M,              _Ir,             _Ig,             _Mu,             _Iru,            _Igu
        (log2(0.9848843), log2(0.0000735), log2(0.0000334), log2(0.0150068), log2(0.0000017), log2(0.0000003)),  # M_
        (log2(0.5196194), log2(0.4647955), log2(0.0059060), log2(0.0096792)),  # Ir_
        (log2(0.7611255), log2(0.2294619), log2(0.0072673), log2(0.0021453)),  # Ig_
        (log2(0.0799009), log2(0.0000262), log2(0.0001836), log2(0.9161349), log2(0.0033370), log2(0.0004173)),  # Mu_
        (log2(0.1434529), log2(0.0036995), log2(0.2642928), log2(0.5885548)),  # Iru_
        (log2(0.1384551), log2(0.0431328), log2(0.6362921), log2(0.1821200))  # Igu_
    )

    default_scoring_matrix = [log2(0.955), log2(0.04), log2(0.004), log2(0.001)]

    def __cinit__(self, CpCountgraph count_graph, int trusted_cov_cutoff, double bits_theta, str filename=None):

        if filename is None:
            self._trans_probs = default_trans_probs
            self._scoring_matrix = default_scoring_matrix
        else:
            with open(filename, 'r') as paramfile:
                params = json.load(paramfile)
            self._trans_probs = params['transition_probabilities']
            self._scoring_matrix = params['scoring_matrix']

        cdef CpReadAligner aligner(count_graph, trusted_cov_cutoff, bits_theta, self._scoring_matrix, self._trans_probs)
        self._aligner = aligner

    def align(self, str sequence):
        cdef CpAlignment alignment = self._aligner.Align(sequence)
        return (alignment.graph_alignment, alignment.read_alignment,
                alignment.trusted, alignment.covs, alignment.store,
                alignment.truncated)

    def align_forward(self, str sequence):
        cdef CpAlignment alignment = self._aligner.AlignForward(sequence)
        return (alignment.graph_alignment, alignment.read_alignment,
                alignment.trusted, alignment.covs, alignment.store,
                alignment.truncated)

    @property
    def scoring_matrix(self):
        return self._scoring_matrix

    @property
    def transition_probabilities(self):
        return self._trans_probs
