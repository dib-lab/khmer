from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.memory cimport unique_ptr, shared_ptr, weak_ptr
from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libc.stdint cimport uintptr_t as size_t

from khmer._oxli.oxli_types cimport *
from khmer._oxli.graphs cimport CpCountgraph, Countgraph

cdef extern from "oxli/read_aligner.hh" namespace "oxli" nogil:

    cdef enum State:
        MATCH
        INSERT_READ
        INSERT_GRAPH
        MATCH_UNTRUSTED
        INSERT_READ_UNTRUSTED
        INSERT_GRAPH_UNTRUSTED

    cdef enum Transition:
        MM, MIr, MIg, MMu, MIru, MIgu,
        IrM, IrIr, IrMu, IrIru,
        IgM, IgIg, IgMu, IgIgu,
        MuM, MuIr, MuIg, MuMu, MuIru, MuIgu,
        IruM, IruIr, IruMu, IruIru,
        IguM, IguIg, IguMu, IguIgu,
        disallowed

    cdef enum Nucl:
        A, C, G, T

    cdef struct AlignmentNode:
        AlignmentNode(AlignmentNode* _prev, Nucl _emission, size_t _seq_idx,
                      State _state, Transition _trans, HashIntoType _fwd_hash,
                      HashIntoType _rc_hash, size_t _length)
        AlignmentNode * prev
        Nucl base
        size_t seq_idx
        State state
        Transition trans
        HashIntoType fwd_hash
        HashIntoType rc_hash

        double score
        double f_score
        double h_score
        bool trusted
        BoundedCounterType cov

        size_t num_indels
        size_t length

        #bool operator==(const AlignmentNode&)
        #bool operator<(const AlignmentNode&)

    ctypedef struct ScoringMatrix:
        const double trusted_match
        const double trusted_mismatch
        const double untrusted_match
        const double untrusted_mismatch

        const double* tsc

        ScoringMatrix(double trusted_match, double trusted_mismatch,
                      double untrusted_match, double untrusted_mismatch,
                      double* trans)

    cdef struct Alignment:
        string graph_alignment
        string read_alignment
        string trusted
        vector[BoundedCounterType] covs
        double score
        bool truncated

    cdef double trans_default []
    cdef double freq_default []

    cdef cppclass CpReadAligner "oxli::ReadAligner":
        
        Alignment* Align(const string&)
        Alignment* AlignForward(const string&)

        CpReadAligner(CpCountgraph *, BoundedCounterType trusted_cutoff,
                      double bits_theta)

        CpReadAligner(CpCountgraph *, BoundedCounterType trusted_cutoff,
                      double* scoring_matrix, double* transitions)
        ScoringMatrix getScoringMatrix()


cdef class ReadAligner:
    cdef shared_ptr[CpReadAligner] _aln_this
    cdef public Countgraph graph
    @staticmethod
    cdef double * _default_transition_probabilities()
    @staticmethod
    cdef double * _default_scoring_matrix()
    @staticmethod
    cdef list _format_scoring_matrix(double *)
    @staticmethod
    cdef tuple _format_transition_probabilities(double *)
