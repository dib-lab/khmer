# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref

from libcpp.memory cimport make_shared
from libcpp.string cimport string

import json
from math import log

from .graphs cimport Countgraph

cdef class ReadAligner:
    '''
    Sequence to graph aligner.

    ReadAligner uses a Countgraph (the counts of k-mers in the target DNA
    sequences) as an implicit De Bruijn graph. Input DNA sequences are aligned
    to this graph via a paired Hidden Markov Model.

    The HMM is configured upon class instantiation; default paramaters for the
    HMM are provided in 'defaultTransitionProbablitites' and
    'defaultScoringMatrix'.

    The main method is 'align'.
    '''

    defaultTransitionProbabilities = (  # _M, _Ir, _Ig, _Mu, _Iru, _Igu
        (log(0.9848843, 2), log(0.0000735, 2), log(0.0000334, 2),
         log(0.0150068, 2), log(0.0000017, 2), log(0.0000003, 2)),  # M_
        (log(0.5196194, 2), log(0.4647955, 2), log(0.0059060, 2),
         log(0.0096792, 2)),  # Ir_
        (log(0.7611255, 2), log(0.2294619, 2), log(0.0072673, 2),
         log(0.0021453, 2)),  # Ig_
        (log(0.0799009, 2), log(0.0000262, 2), log(0.0001836, 2),
         log(0.9161349, 2), log(0.0033370, 2), log(0.0004173, 2)),  # Mu_
        (log(0.1434529, 2), log(0.0036995, 2), log(0.2642928, 2),
         log(0.5885548, 2)),  # Iru_
        (log(0.1384551, 2), log(0.0431328, 2), log(0.6362921, 2),
         log(0.1821200, 2))  # Igu_
    )

    defaultScoringMatrix = [
        log(0.955, 2), log(0.04, 2), log(0.004, 2), log(0.001, 2)]


    def __cinit__(self, Countgraph count_graph, double trusted_cov_cutoff=2, 
                       double bits_theta=1, *args, **kwargs):
        '''
        Initialize ReadAligner.

        HMM state notation abbreviations:
        M_t - trusted match; M_u - untrusted match
        Ir_t - trusted read insert; Ir_u - untrusted read insert
        Ig_t - trusted graph insert; Ig_u - untrusted graph insert

        Keyword arguments:
        filename - a path to a JSON encoded file providing the scoring matrix
            for the HMM in an entry named 'scoring_matrix' and the transition
            probabilities for the HMM in an entry named
            'transition_probabilities'. If provided the remaining keyword
            arguments are ignored. (default: None)
        scoring_matrix - a list of floats: trusted match, trusted mismatch,
            unstrusted match, untrusted mismatch. (default:
                ReadAligner.defaultScoringMatrix)
        transition_probabilities - A sparse matrix as a tuple of six tuples.
            The inner tuples contain 6, 4, 4, 6, 4, and 4 floats respectively.
            Transition are notated as 'StartState-NextState':
            (
              ( M_t-M_t,  M_t-Ir_t,  M_t-Ig_t,  M_t-M_u,  M_t-Ir_u,  M_t-Ig_u),
              (Ir_t-M_t, Ir_t-Ir_t,            Ir_t-M_u, Ir_t-Ir_u           ),
              (Ig_t-M_t,          , Ig_t-Ig_t, Ig_t-M_u,            Ig_t-Ig_u),
              ( M_u-M_t,  M_u-Ir_t,  M_u-Ig_t,  M_u-M_u,  M_u-Ir_u,  M_u-Ig_u),
              (Ir_u-M_t, Ir_u-Ir_t,            Ir_u-M_u, Ir_u-Ir_u           ),
              (Ig_u-M_t,          , Ig_u-Ig_t, Ig_u-M_u,            Ig_u-Ig_u)
            )
            (default: ReadAligner.defaultTransitionProbabilities)

        '''
        if 'filename' in kwargs:
            with open(kwargs.pop('filename')) as paramfile:
                params = json.load(paramfile)
            scoring_matrix = params['scoring_matrix']
            transition_probabilities = params['transition_probabilities']
        else:
            if 'scoring_matrix' in kwargs:
                scoring_matrix = kwargs.pop('scoring_matrix')
            else:
                scoring_matrix = ReadAligner.defaultScoringMatrix
            if 'transition_probabilities' in kwargs:
                transition_probabilities = kwargs.pop('transition_probabilities')
            else:
                transition_probabilities = \
                    ReadAligner.defaultTransitionProbabilities

        cdef vector[double] _transitions = [x for sublist in transition_probabilities for x in sublist]
        cdef vector[double] _scoring_matrix = scoring_matrix

        self._aln_this = make_shared[CpReadAligner](count_graph._cg_this.get(),
                                                    trusted_cov_cutoff,
                                                    bits_theta,
                                                    &(_scoring_matrix[0]),
                                                    &(_transitions[0]))
        self.graph = count_graph

    def align(self, str sequence, as_dict=False):
        cdef string _sequence = self.graph._valid_sequence(sequence)
        cdef Alignment * aln = deref(self._aln_this).Align(_sequence)
        cdef object score = aln[0].score
        cdef object alignment = aln[0].graph_alignment
        cdef object read_alignment = aln[0].read_alignment
        cdef object truncated = aln[0].truncated

        return score, alignment.upper(), read_alignment.upper(), truncated


    def align_forward(self, str sequence):
        cdef string _sequence = self.graph._valid_sequence(sequence)
        cdef Alignment * aln = deref(self._aln_this).AlignForward(_sequence)
        cdef object score = aln[0].score
        cdef object alignment = aln[0].graph_alignment
        cdef object read_alignment = aln[0].read_alignment
        cdef object truncated = aln[0].truncated
        cdef list covs = aln[0].covs

        return (score, alignment.upper(), read_alignment.upper(),
                truncated, covs)

    @property
    def scoring_matrix(self):
        '''Get the scoring matrix in use.

        Returns a tuple of floats: 
        (trusted_match, trusted_mismatch, untrusted_match, untrusted_mismatch)
        '''
        cdef ScoringMatrix matrix = deref(self._aln_this).getScoringMatrix()
        return [matrix.trusted_match, 
                matrix.trusted_mismatch,
                matrix.untrusted_match, 
                matrix.untrusted_mismatch]

    @property
    def transition_probabilities(self):
        '''Get the transition probabilties in use.

        HMM state notation abbreviations:
            M_t - trusted match; M_u - untrusted match
            Ir_t - trusted read insert; Ir_u - untrusted read insert
            Ig_t - trusted graph insert; Ig_u - untrusted graph insert

        Returns a sparse matrix as a tuple of six tuples.
        The inner tuples contain 6, 4, 4, 6, 4, and 4 floats respectively.
        Transition are notated as 'StartState-NextState':
        (
          ( M_t-M_t,  M_t-Ir_t,  M_t-Ig_t,  M_t-M_u,  M_t-Ir_u,  M_t-Ig_u),
          (Ir_t-M_t, Ir_t-Ir_t,            Ir_t-M_u, Ir_t-Ir_u           ),
          (Ig_t-M_t,          , Ig_t-Ig_t, Ig_t-M_u,            Ig_t-Ig_u),
          ( M_u-M_t,  M_u-Ir_t,  M_u-Ig_t,  M_u-M_u,  M_u-Ir_u,  M_u-Ig_u),
          (Ir_u-M_t, Ir_u-Ir_t,            Ir_u-M_u, Ir_u-Ir_u           ),
          (Ig_u-M_t,          , Ig_u-Ig_t, Ig_u-M_u,            Ig_u-Ig_u)
        )'''
        
        cdef ScoringMatrix matrix = deref(self._aln_this).getScoringMatrix()
        return ((matrix.tsc[0], matrix.tsc[1], matrix.tsc[2],
                 matrix.tsc[3], matrix.tsc[4], matrix.tsc[5]),

                (matrix.tsc[6], matrix.tsc[7], matrix.tsc[8],
                 matrix.tsc[9]), 

                (matrix.tsc[10], matrix.tsc[11], matrix.tsc[12], 
                 matrix.tsc[13]), 

                (matrix.tsc[14], matrix.tsc[15], matrix.tsc[16], 
                 matrix.tsc[17], matrix.tsc[18], matrix.tsc[19]), 
                
                (matrix.tsc[20], matrix.tsc[21], matrix.tsc[22], 
                 matrix.tsc[23]),
                
                (matrix.tsc[24], matrix.tsc[25], matrix.tsc[26],
                 matrix.tsc[27]))

