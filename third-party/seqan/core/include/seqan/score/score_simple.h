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
// Code for the Simple Scoring Schema.
// ==========================================================================

#ifndef SEQAN_SCORE_SCORE_SIMPLE_H_
#define SEQAN_SCORE_SCORE_SIMPLE_H_

namespace SEQAN_NAMESPACE_MAIN {

/**
.Spec.Simple Score
..cat:Scoring
..signature:Score<TValue, Simple>
..param.TValue:The value type.
...default:int
..general:Class.Score
..summary:Simple scoring scheme that has scores for matches, mismatches, opening gaps and extending gaps.
..include:seqan/score.h
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

/**
.Memfunc.Score#Score:
..class:Class.Score
..summary:Constructor
..signature:Score<TValue, Simple>()
..signature:Score<TValue, Simple>(score)
..signature:Score<TValue, Simple>(match, mismatch, gap [, gap_open])
..param.score:Other Score object. (copy constructor)
..param.match:TValue object.
...default:0
..param.mismatch:TValue object.
...default:-1
..param.gap:TValue object.
...remarks:The score for a single blank in a gap (linear gap costs).
...default:-1
..param.gap_open:TValue object.
...remarks:The score for the first blank in a gap (affine gap costs).
...default:$gap$
..remarks:
...text:If both gap and gap_open are specified, the total score of a length $n$ gap is $gap_open (n-1)*gap$.
...note:Usually $mismatch$, $gap$, and $gap_open$ are negative values.
 */
    Score(TValue _match, TValue _mismatch, TValue _gap_extend, TValue _gap_open)
        : data_match(_match), data_mismatch(_mismatch),
          data_gap_extend(_gap_extend), data_gap_open(_gap_open) {
        SEQAN_CHECKPOINT;
    }
};


/**
.Shortcut.SimpleScore:
..cat:Scoring
..summary:Simple scoring scheme.
..signature:SimpleScore
..shortcutfor:Spec.Simple Score
...signature:Score<int, Simple>
..see:Spec.Simple Score
*/
typedef Score<int, Simple> SimpleScore;


/**
.Function.scoreMatch:
..class:Class.Score
..cat:Alignments
..summary:Match score.
..signature:scoreMatch(object)
..param.object.type:Class.Score
..returns:Match score.
..see:Function.scoreMismatch
..see:Function.scoreGapExtend
..see:Function.scoreGapOpen
..include:seqan/score.h
*/
template <typename TValue, typename TSpec>
inline TValue
scoreMatch(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return me.data_match;
}

/**
.Function.setScoreMatch:
..class:Class.Score
..cat:Alignments
..summary:Set match score.
..signature:setScoreMatch(object, value)
..param.object.type:Class.Score
..param.value.type:nolink:TValue
..see:Function.scoreMismatch
..see:Function.scoreGapExtend
..see:Function.scoreGapOpen
..include:seqan/score.h
*/
template <typename TValue, typename TSpec>
inline void
setScoreMatch(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_match = value;
}

/**
.Function.scoreMismatch:
..class:Class.Score
..cat:Alignments
..summary:Mismatch score.
..signature:scoreMismatch(object)
..param.object.type:Class.Score
..returns:Mismatch score.
...note:Usually, mismatches have negative scores.
..see:Function.scoreMatch
..see:Function.scoreGapExtend
..see:Function.scoreGapOpen
..include:seqan/score.h
*/
template <typename TValue, typename TSpec>
inline TValue
scoreMismatch(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return me.data_mismatch;
}


/**
.Function.setScoreMismatch:
..class:Class.Score
..cat:Alignments
..summary:Set mismatch score.
..signature:setScoreMismatch(object, value)
..param.object.type:Class.Score
..param.value.type:nolink:TValue
..see:Function.scoreMismatch
..see:Function.scoreGapExtend
..see:Function.scoreGapOpen
..include:seqan/score.h
*/
template <typename TValue, typename TSpec>
inline void
setScoreMismatch(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_mismatch = value;
}


/**
.Function.scoreGapExtend:
..class:Class.Score
..cat:Alignments
..summary:Score for extending gaps.
..signature:scoreGapExtend(object)
..param.object.type:Class.Score
..returns:Score for extending gaps.
...note:Usually, gaps have negative scores.
..see:Function.scoreMismatch
..see:Function.scoreMatch
..see:Function.scoreGapOpen
..include:seqan/score.h
*/
template <typename TValue, typename TSpec>
inline TValue
scoreGapExtend(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return me.data_gap_extend;
}


/**
.Function.setScoreGapExtend:
..class:Class.Score
..cat:Alignments
..summary:Set gap extension score.
..signature:setScoreGapExtend(object, value)
..param.object.type:Class.Score
..param.value.type:nolink:TValue
..see:Function.scoreMismatch
..see:Function.scoreGapExtend
..see:Function.scoreGapOpen
..include:seqan/score.h
*/
template <typename TValue, typename TSpec>
inline void
setScoreGapExtend(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_gap_extend = value;
}


/**
.Function.scoreGapOpen:
..class:Class.Score
..cat:Alignments
..summary:Score for opening a gap.
..signature:scoreGapOpen(object)
..param.object.type:Class.Score
..returns:Score for opening a gap.
...note:Usually, gaps have negative scores.
..see:Function.scoreMismatch
..see:Function.scoreGapExtend
..see:Function.scoreMatch
..include:seqan/score.h
*/
template <typename TValue, typename TSpec>
inline TValue
scoreGapOpen(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return me.data_gap_open;
}


/**
.Function.setScoreGapOpen:
..class:Class.Score
..cat:Alignments
..summary:Set gap opening score.
..signature:setScoreGapOpen(object, value)
..param.object.type:Class.Score
..param.value.type:nolink:TValue
..see:Function.scoreMismatch
..see:Function.scoreGapExtend
..see:Function.scoreGapOpen
..include:seqan/score.h
*/
template <typename TValue, typename TSpec>
inline void
setScoreGapOpen(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_gap_open = value;
}


/**
.Function.scoreGap:
..class:Class.Score
..cat:Alignments
..summary:Score for gaps.
..signature:scoreGapExtend(object)
..param.object.type:Class.Score
..returns:Score for extending gaps.
...note:Usually, gaps have negative scores.
..remarks:This score is used for linear gap costs. For affine gap costs use @Function.scoreGapExtend@ and @Function.scoreGapOpen@ instead.
..see:Function.scoreMismatch
..see:Function.scoreMatch
..see:Function.scoreGapOpen
..see:Function.scoreGapExtend
..include:seqan/score.h
*/
// TODO(holtgrew): This shortcut/forward should live in score_base.h.
template <typename TValue, typename TSpec>
inline TValue
scoreGap(Score<TValue, TSpec> const & me) {
    SEQAN_CHECKPOINT;
    return scoreGapExtend(me);
}


/**
.Function.setScoreGap:
..class:Class.Score
..cat:Alignments
..summary:Set gap opening and extension score.
..signature:setScoreGap(object, value)
..param.object.type:Class.Score
..param.value.type:nolink:TValue
..see:Function.scoreMismatch
..see:Function.scoreGapExtend
..see:Function.scoreGapOpen
..include:seqan/score.h
*/
template <typename TValue, typename TSpec>
inline void
setScoreGap(Score<TValue, TSpec> & me, TValue const & value) {
    SEQAN_CHECKPOINT;
    me.data_gap_open = value;
    me.data_gap_extend = value;
}


// TODO(rmaerker): Remove this here!
/**
.Function.score:
..cat:Scoring
..class:Class.Score
..signature:score(score, value1, value2) [deprecated]
..param.score:A scoring scheme.
...type:Class.Score
..param.value1:Aligned character in the first sequence.
..param.value2:Aligned character in the second sequence.
..status:The version that explicitely requires two characters is deprecated. For using (or defining) custom scoring schemes better use the functions
@Function.score@, @Function.scoreGapOpenVertical@, @Function.scoreGapExtendVertical@, @Function.scoreGapOpenHorizontal@ and @Function.scoreGapExtendHorizontal@ 
as they require character positions and sequences and allow for position-dependent scoring schemes.
..include:seqan/score.h
*/
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

#endif  // SEQAN_SCORE_SCORE_SIMPLE_H_
