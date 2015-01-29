// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SEQAN_SCORE_SCORE_MATRIX_H_
#define SEQAN_SCORE_SCORE_MATRIX_H_

// TODO(holtgrew): If the complex type conversions are necessary, a static_cast<> is more C++ and explicit.

namespace SEQAN_NAMESPACE_MAIN {

template <typename TValue, typename TSequenceValue, typename TSpec>
struct ScoringMatrixData_;


template <typename TSequenceValue = AminoAcid, typename TSpec = Default>
struct ScoreMatrix;


/**
.Spec.Score Matrix:
..cat:Scoring
..summary:A general scoring matrix.
..general:Class.Score
..signature:Score<TValue, ScoreMatrix<TSequenceValue, TSpec> >
..param.TValue:Type of the score values.
...default:$int$
..param.TSequenceValue:Type of alphabet underlying the matrix.
...default:$AminoAcid$
..include:seqan/score.h
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

    /**
.Memfunc.Score Matrix#Score
..cat:Scoring
..summary:Constructor.
..class:Spec.Score Matrix
..signature:Score(gapExtend)
..param.gapExtend:The gap extension penalty.
...remark:TValue
     */
    explicit Score(TValue _gap_extend = -1)
        : data_gap_extend(_gap_extend),
          data_gap_open(_gap_extend) {
        SEQAN_CHECKPOINT;
        setDefaultScoreMatrix(*this, TSpec());
    }

    /**
.Memfunc.Score Matrix#Score
..signature:Score(gapExtend, gapOpen)
..param.gapOpen:The gap open penalty.
...remark:TValue
     */
    Score(TValue _gap_extend, TValue _gap_open)
        : data_gap_extend(_gap_extend), data_gap_open(_gap_open) {
        SEQAN_CHECKPOINT;
        setDefaultScoreMatrix(*this, TSpec());
    }

    /**
.Memfunc.Score Matrix#Score
..signature:Score(filename, gapExtend)
..param.filename:The path to the file to load.
...type:Class.String
..see:Function.loadScoreMatrix
     */
    template <typename TString>
    Score(TString const & filename, TValue _gap_extend = -1)
        : data_gap_extend(_gap_extend), data_gap_open(_gap_extend) {
        SEQAN_CHECKPOINT;
        loadScoreMatrix(*this, filename);
    }

    /**
.Memfunc.Score Matrix#Score
..signature:Score(filename, gapExtend, gapOpen)
     */
    template <typename TString>
    Score(TString const & filename, TValue _gap_extend, TValue _gap_open)
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


/**
.Function.setScore:
..class:Spec.Score Matrix
..cat:Scoring
..summary:Set the substitution score between two values.
..signature:setScore(scoreMatrix, val1, val2, score)
..param.scoreMatrix:
...type:Spec.Score Matrix
..param.val1:First value.
..param.val2:Second value.
..param.score:The value to set the score to.
..include:seqan/score.h
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


/**
.Function.setDefaultScoreMatrix:
..cat:Scoring
..summary:Set the value of the given matrix to the default value.
..signature:setDefaultScoreMatrix(scoreMatrix, tag)
..param.scoreMatrix:The @Spec.Score Matrix@ to set.
...type:Spec.Score Matrix
..param.tag:The tag to specify the matrix.
...type:Shortcut.Blosum30
...type:Shortcut.Blosum62
...type:Shortcut.Blosum80
..include:seqan/score.h
 */
template <typename TValue, typename TSequenceValue, typename TSpec, typename TTag>
inline void
setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, TTag) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    TValue const * tab = ScoringMatrixData_<TValue, TSequenceValue, TTag>::getData();
    arrayCopy(tab, tab + TScore::TAB_SIZE, sc.data_tab);
}


/**
.Function.setDefaultScoreMatrix
..param.tag:
...type:Tag.Default
...remark:If @Tag.Default@, then the matrix will be filled with default constructed $TValue$ values.
..include:seqan/score.h
 */
template <typename TValue, typename TSequenceValue, typename TSpec>
inline void
setDefaultScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, Default) {
    SEQAN_CHECKPOINT;
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;
    arrayFill(sc.data_tab, sc.data_tab + TScore::TAB_SIZE, TValue());
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_SCORE_SCORE_MATRIX_H_
