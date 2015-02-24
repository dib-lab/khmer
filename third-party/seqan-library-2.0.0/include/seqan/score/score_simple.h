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
// Code for the Simple Scoring Schema.
// ==========================================================================

#ifndef SEQAN_SSCORE_SIMPLE_H_
#define SEQAN_SSCORE_SIMPLE_H_

namespace seqan {

/*!
 * @class SimpleScore
 * @extends Score
 * @headerfile <seqan/score.h>
 * @brief Simple scoring scheme that has scores for matches, mismatches, opening gaps and extending gaps.
 *
 * @signature template <typename TValue>
 *            class Score<TValue, Simple>;
 *
 * @tparam TValue The score value to use.
 */

template <typename TValue>
class Score<TValue, Simple> {
public:
    // The score for a match.
    TValue data_match;

    // The score for a mismatch.
    TValue data_mismatch;

    // The gap extension score.
    TValue data_gap_extend;

    // The gap open score.
    TValue data_gap_open;

    Score()
        : data_match(0), data_mismatch(-1), data_gap_extend(-1),
          data_gap_open(-1) {
        SEQAN_CHECKPOINT;
    }

    Score(TValue _match, TValue _mismatch, TValue _gap)
        : data_match(_match), data_mismatch(_mismatch),
          data_gap_extend(_gap), data_gap_open(_gap) {
        SEQAN_CHECKPOINT;
    }

/*!
 * @fn SimpleScore::Score
 * @brief Constructor
 *
 * @signature Score::Score();
 * @signature Score::Score(score);
 * @signature Score::Score(match, mismatch, gap[, gapOpen]);
 *
 * @param[in] score   Other SimpleScore object to copy from.
 * @param[in] match   Match score value, type TValue, default 0.
 * @param[in] msmatch Mismatch score value, type TValue, default -1.
 * @param[in] gap     Gap extension value, type TValue, default -1.
 * @param[in] gapOpen Gap open value (defaults to gap), type TValue.
 */

    Score(TValue _match, TValue _mismatch, TValue _gap_extend, TValue _gap_open)
        : data_match(_match), data_mismatch(_mismatch),
          data_gap_extend(_gap_extend), data_gap_open(_gap_open) {
        SEQAN_CHECKPOINT;
    }
};

/*!
 * @typedef SimpleScoreTypedef SimpleScore
 * @headerfile <seqan/score.h>
 * @brief Simple scoring scheme.
 *
 * @signature typedef Score<int, Simple> SimpleScore;
 */

typedef Score<int, Simple> SimpleScore;

/*!
 * @fn SimpleScore#scoreMatch
 * @brief Match score
 *
 * @signature TValue scoreMatch(score);
 *
 * @param[in] score The SimpleScore scoring scheme.
 *
 * @return TValue The match score.
 */

template <typename TValue, typename TSpec>
inline TValue
scoreMatch(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return me.data_match;
}

/*!
 * @fn SimpleScore#setScoreMatch
 * @brief Set match score.
 *
 * @signature void setScoreMatch(score, value);
 *
 * @param[in,out] score The SimpleScore scoring scheme to set the value for.
 * @param[in]     value The value to set the match score to.
 */

template <typename TValue, typename TSpec>
inline void
setScoreMatch(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_match = value;
}

/*!
 * @fn SimpleScore#scoreMismatch
 * @brief Set mismatch score.
 *
 * @signature TValue scoreMismatch(score);
 *
 * @param[in] score The SimpleScore to query for its mismatch score.
 *
 * @return TValue The mismatch score.
 */

template <typename TValue, typename TSpec>
inline TValue
scoreMismatch(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return me.data_mismatch;
}

/*!
 * @fn SimpleScore#setScoreMismatch
 * @brief Set mismatch score.
 *
 * @signature void setScoreMismatch(score, value);
 *
 * @param[in,out] score The SimpleScore scoring scheme to set the mismatch value for.
 * @param[in]     value The value to set the mismatch score to.
 */

template <typename TValue, typename TSpec>
inline void
setScoreMismatch(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_mismatch = value;
}

/*!
 * @fn SimpleScore#scoreGapExtend
 * @brief Set gap extension score.
 *
 * @signature TValue scoreGapExtend(score);
 *
 * @param[in] score The SimpleScore to query for its gap extension score.
 *
 * @return TValue The gap extension score.
 */

template <typename TValue, typename TSpec>
inline TValue
scoreGapExtend(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return me.data_gap_extend;
}


template <typename TValue, typename TSpec>
inline void
setScoreGapExtend(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_gap_extend = value;
}

/*!
 * @fn SimpleScore#scoreGapOpen
 * @brief Set gap open score.
 *
 * @signature TValue scoreGapOpen(score);
 *
 * @param[in] score The SimpleScore to query for its gap open score.
 *
 * @return TValue The gap open score.
 */

template <typename TValue, typename TSpec>
inline TValue
scoreGapOpen(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return me.data_gap_open;
}


/*!
 * @fn SimpleScore#setScoreGapOpen
 * @brief Set gap open score.
 *
 * @signature void setScoreGapOpen(score, value);
 *
 * @param[in,out] score The SimpleScore scoring scheme to set the gap open value for.
 * @param[in]     value The value to set the gap open score to.
 */

template <typename TValue, typename TSpec>
inline void
setScoreGapOpen(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_gap_open = value;
}

/*!
 * @fn SimpleScore#scoreGap
 * @brief Set gap score.
 *
 * @signature TValue scoreGap(score);
 *
 * @param[in] score The SimpleScore to query for its gap score.
 *
 * @return TValue The gap score.
 */


// TODO(holtgrew): This shortcut/forward should live in score_base.h.
template <typename TValue, typename TSpec>
inline TValue
scoreGap(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return scoreGapExtend(me);
}

/*!
 * @fn SimpleScore#setScoreGap
 * @brief Set gap score.
 *
 * @signature void setScoreGap(score, value);
 *
 * @param[in,out] score The SimpleScore scoring scheme to set the gap value for.
 * @param[in]     value The value to set the gap score to.
 */

template <typename TValue, typename TSpec>
inline void
setScoreGap(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_gap_open = value;
    me.data_gap_extend = value;
}


// TODO(rmaerker): Remove this here!
//template <typename TValue, typename TSpec, typename TVal1, typename TVal2>
//inline TValue
//score(Score<TValue, TSpec> const & me, TVal1 left, TVal2 right) {
//    SEQAN_CHECKPOINT;
//    if (left == right)
//        return scoreMatch(me);
//    else
//        return scoreMismatch(me);
//}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_SSCORE_SIMPLE_H_
