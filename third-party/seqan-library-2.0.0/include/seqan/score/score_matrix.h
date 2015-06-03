// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Code for score matrices with data from files or built-in data.
// ==========================================================================

#ifndef SEQAN_SSCORE_MATRIX_H_
#define SEQAN_SSCORE_MATRIX_H_

// TODO(holtgrew): If the complex type conversions are necessary, a static_cast<> is more C++ and explicit.

namespace seqan {

template <typename TValue, typename TSequenceValue, typename TSpec>
struct ScoringMatrixData_;


template <typename TSequenceValue = AminoAcid, typename TSpec = Default>
struct ScoreMatrix;

/*!
 * @class MatrixScore
 * @headerfile <seqan/score.h>
 * @extends Score
 * @brief A general scoring matrix.
 *
 * @signature template <typename TValue, typename TSeqValue, typename TSpec>
 *            class Score<TValue, ScoreMatrix<[TSeqValue[, TSpec]]> >;
 *
 * @tparam TValue    The score value.
 * @tparam TSeqValue The alphabet type, defaults to AminoAcid.
 * @tparam TSpec     Further specialization, defaults to Default.
 *
 * The TSpec argument can be used to obtain a predefined matrix.
 * Specify one of the following tags:
 *
 * ScoreSpecBlosum30, ScoreSpecBlosum45, ScoreSpecBlosum62, ScoreSpecBlosum80,
 * ScoreSpecPam40, ScoreSpecPam120, ScoreSpecPam200, ScoreSpecPam250, ScoreSpecVtml200.
 *
 * This will internally call @link MatrixScore#setDefaultScoreMatrix setDefaultScoreMatrix@endlink.
 *
 * In order to provide a more user-friendly access to the predefined scoring matrixes, typedefs exist:
 * @link Blosum30 @endlink, @link Blosum45 @endlink,  @link Blosum62 @endlink,
 * @link Blosum80 @endlink, @link Pam40 @endlink,     @link Pam120 @endlink,
 * @link Pam200 @endlink,   @link Pam250 @endlink and @link Vtml200 @endlink.
 *
 * @fn MatrixScore::Score
 * @brief Constructor
 *
 * @signature MatrixScore::Score(gapExtend[, gapOpen]);
 * @signature MatrixScore::Score(fileName, gapExtend[, gapOpen]);
 *
 * @param[in] fileName  Path to load the matrix from, type is <tt>char const *</tt>.
 * @param[in] gapExtend Gap extension score, type is TValue.
 * @param[in] gapOpen   Gap open score, defaults to gapExtend, type is TValue.
 */

template <typename TValue, typename TSequenceValue, typename TSpec>
class Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > {
public:
    // Static computation of the required array size.
    enum {
        VALUE_SIZE = ValueSize<TSequenceValue>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    // The data table.
    TValue data_tab[TAB_SIZE];

    // The gap extension score.
    TValue data_gap_extend;

    // The gap open score.
    TValue data_gap_open;

    explicit Score(TValue _gap_extend = -1)
        : data_gap_extend(_gap_extend),
          data_gap_open(_gap_extend) {
        SEQAN_CHECKPOINT;
        setDefaultScoreMatrix(*this, TSpec());
    }

    Score(TValue _gap_extend, TValue _gap_open)
        : data_gap_extend(_gap_extend), data_gap_open(_gap_open) {
        SEQAN_CHECKPOINT;
        setDefaultScoreMatrix(*this, TSpec());
    }

    explicit Score(char const * filename, TValue _gap_extend = -1)
        : data_gap_extend(_gap_extend), data_gap_open(_gap_extend) {
        SEQAN_CHECKPOINT;
        loadScoreMatrix(*this, filename);
    }

    Score(char const * filename, TValue _gap_extend, TValue _gap_open)
        : data_gap_extend(_gap_extend), data_gap_open(_gap_open) {
        SEQAN_CHECKPOINT;
        loadScoreMatrix(*this, filename);
    }
};


// TODO(holtgrew): Does it make sense to document each Score specialization?  Should dddoc show a list of all specializations of a class?
template <typename TValue, typename TSequenceValue, typename TSpec, typename TVal1, typename TVal2>
inline TValue
score(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc, TVal1 val1, TVal2 val2) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    // TODO(holtgrew): Why not implicit cast?
    unsigned int i = (TSequenceValue) val1;  // conversion TVal1 => TSequenceValue => integral
    unsigned int j = (TSequenceValue) val2;  // conversion TVal2 => TSequenceValue => integral
    return sc.data_tab[i * TScore::VALUE_SIZE + j];
}


/*!
 * @fn MatrixScore#setScore
 * @brief Set the substitution score between to values.
 *
 * @signature void setScore(score, x, y, v);
 *
 * @param[in,out] score The MatrixScore to set the value for.
 * @param[in]     x     The substituted alphabet value.
 * @param[in]     y     The alphabet value to substitute x for.
 * @param[in]     v     The score value to set.
 */

template <typename TValue, typename TSequenceValue, typename TSpec, typename TVal1, typename TVal2, typename T>
inline void
setScore(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, TVal1 val1, TVal2 val2, T score) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    // TODO(holtgrew): Why not implicit cast?
    unsigned int i = (TSequenceValue) val1;  // conversion TVal1 => TSequenceValue => integral
    unsigned int j = (TSequenceValue) val2;  // conversion TVal2 => TSequenceValue => integral
    sc.data_tab[i * TScore::VALUE_SIZE + j] = score;
}

/*!
 * @fn MatrixScore#setDefaultScoreMatrix
 * @brief Set the score matrix of a Score to one of the default matrixes.
 *
 * @signature void setDefaultScoreMatrix(score, tag);
 *
 * @param[in,out] score The MatrixScore to update.
 * @param[in]     tag   The tag to select the default matrix, see description below.
 *
 * @section Remarks
 *
 * The tag must be one of the following:
 * Default, ScoreSpecBlosum30, ScoreSpecBlosum45, ScoreSpecBlosum62, ScoreSpecBlosum80,
 * ScoreSpecPam40, ScoreSpecPam120, ScoreSpecPam200, ScoreSpecPam250, ScoreSpecVtml200.
 *
 * If Default is used for tag then the matrix will be filled with default-constructed TValue values.
 */

template <typename TValue, typename TSequenceValue, typename TSpec, typename TTag>
inline void
setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, TTag) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    TValue const * tab = ScoringMatrixData_<TValue, TSequenceValue, TTag>::getData();
    arrayCopy(tab, tab + TScore::TAB_SIZE, sc.data_tab);
}


template <typename TValue, typename TSequenceValue, typename TSpec>
inline void
setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, Default) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    arrayFill(sc.data_tab, sc.data_tab + TScore::TAB_SIZE, TValue());
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_SSCORE_MATRIX_H_
