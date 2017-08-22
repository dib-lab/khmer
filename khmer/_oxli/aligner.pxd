from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from graphs cimport CpHashgraph
from oxli_types cimport *

cdef extern from "oxli/read_aligner.hh" namespace "oxli":

    ctypedef struct CpAlignment "oxli::Alignment":
        string graph_alignment
        string read_alignment
        string trusted
        vector[BoundedCounterType] covs
        double score
        bool truncated

    ctypedef struct CpScoringMatrix "oxli::ScoringMatrix":
        double trusted_match
        double trusted_mismatch
        double untrusted_match
        double untrusted_mismatch
        double* tsc

    cdef cppclass CpReadAligner "oxli::ReadAligner":
        # CpReadAligner(CpHashgraph, BoundedCounterType, double)  # FIXME I don't think this is used
        CpReadAligner(CpHashgraph, BoundedCounterType, double, double*, double*)

        CpAlignment Align(const string)
        CpAlignment AlignForward(const string)
        CpScoringMatrix getScoringMatrix()
