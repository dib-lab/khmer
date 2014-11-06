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
// Edit distance score class and supporting code.
// ==========================================================================

#ifndef SEQAN_SCORE_SCORE_EDIT_H_
#define SEQAN_SCORE_SCORE_EDIT_H_

namespace SEQAN_NAMESPACE_MAIN {

/**
.Spec.EditDistance
..cat:Scoring
..summary:Edit distance scoring scheme.
..signature:Score<TValue, EditDistance>
..param.TValue:The value type.
...default:int
..general:Class.Score
..remarks:Semantically equivalent to a default contructed @Spec.Simple Score.Score<int, Simple>@.
..remarks:$EditDistance$ is a synonym for @Tag.LevenshteinDistance@.
..include:seqan/score.h
*/

// TODO(holtgrew): Should EditDistance better live here instead of basic_tag.h?
// EditDistance is defined in basic_tag.h
template <typename TValue>
class Score<TValue, EditDistance> {
public:
    Score() {}
};


/**
.Shortcut.EditDistanceScore:
..cat:Scoring
..summary:Edit distance scoring scheme.
..signature:EditDistanceScore
..shortcutfor:Spec.EditDistance
...signature:Score<int, EditDistance>
..see:Spec.EditDistance
..include:seqan/score.h
*/

typedef Score<int, EditDistance> EditDistanceScore;

template <typename TValue>
inline TValue
scoreMatch(Score<TValue, EditDistance> &) {
    SEQAN_CHECKPOINT;
    return 0;
}


template <typename TValue>
inline TValue
scoreMatch(Score<TValue, EditDistance> const &) {
    SEQAN_CHECKPOINT;
    return 0;
}


template <typename TValue>
inline TValue
scoreMismatch(Score<TValue, EditDistance> &) {
    SEQAN_CHECKPOINT;
    return -1;
}


template <typename TValue>
inline TValue
scoreMismatch(Score<TValue, EditDistance> const &) {
    SEQAN_CHECKPOINT;
    return -1;
}


template <typename TValue>
inline TValue
scoreGapExtend(Score<TValue, EditDistance> &) {
    SEQAN_CHECKPOINT;
    return -1;
}


template <typename TValue>
inline TValue
scoreGapExtend(Score<TValue, EditDistance> const &) {
    SEQAN_CHECKPOINT;
    return -1;
}


template <typename TValue>
inline TValue
scoreGapOpen(Score<TValue, EditDistance> &) {
    SEQAN_CHECKPOINT;
    return -1;
}


template <typename TValue>
inline TValue
scoreGapOpen(Score<TValue, EditDistance> const &) {
    SEQAN_CHECKPOINT;
    return -1;
}

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_SCORE_SCORE_EDIT_H_
