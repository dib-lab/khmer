from cython.operator cimport dereference as deref

from libc.math cimport log2
from libcpp.memory cimport make_shared
from libcpp.string cimport string

import json
from math import log

from khmer._oxli.graphs cimport Countgraph
from khmer._oxli.utils cimport _flatten_fill, _fill

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

    def __cinit__(self, Countgraph count_graph, BoundedCounterType trusted_cov_cutoff=2, 
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
        cdef double * _transitions = ReadAligner._default_transition_probabilities()
        cdef double * _scoring_matrix = ReadAligner._default_scoring_matrix()

        if 'filename' in kwargs:
            with open(kwargs.pop('filename')) as paramfile:
                params = json.load(paramfile)
            _fill(_scoring_matrix, params['scoring_matrix'])
            _flatten_fill(_transitions, params['transition_probabilities'])
        else:
            if 'scoring_matrix' in kwargs:
                _fill(_scoring_matrix,
                      kwargs.pop('scoring_matrix'))
            if 'transition_probabilities' in kwargs:
                _flatten_fill(_transitions,
                              kwargs.pop('transition_probabilities'))

        self._aln_this = make_shared[CpReadAligner](count_graph._cg_this.get(),
                                                    trusted_cov_cutoff,
                                                    bits_theta,
                                                    &(_scoring_matrix[0]),
                                                    &(_transitions[0]))
        self.graph = count_graph

    def align(self, str sequence):
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

    @staticmethod
    cdef double * _default_transition_probabilities():
        return trans_default

    @staticmethod
    cdef double * _default_scoring_matrix():
        return freq_default

    @staticmethod
    cdef list _format_scoring_matrix(double * sm):
        return [sm[0], sm[1], sm[2], sm[3]]

    @staticmethod
    cdef tuple _format_transition_probabilities(double * tp):
         return ((tp[0], tp[1], tp[2],
                  tp[3], tp[4], tp[5]),
                 (tp[6], tp[7], tp[8],
                  tp[9]), 
                 (tp[10], tp[11], tp[12], 
                  tp[13]), 
                 (tp[14], tp[15], tp[16], 
                  tp[17], tp[18], tp[19]), 
                 (tp[20], tp[21], tp[22], 
                  tp[23]),
                 (tp[24], tp[25], tp[26],
                  tp[27]))

    @property
    def defaultScoringMatrix(self):
        return ReadAligner._format_scoring_matrix(ReadAligner._default_scoring_matrix())

    @property
    def defaultTransitionProbabilities(self):
        return ReadAligner._format_transition_probabilities(ReadAligner._default_transition_probabilities())

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
        return ReadAligner._format_transition_probabilities(matrix.tsc)

