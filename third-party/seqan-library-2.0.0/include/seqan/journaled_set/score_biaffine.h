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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_JOURNALED_SET_SBIAFFINE_H_
#define INCLUDE_SEQAN_JOURNALED_SET_SBIAFFINE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct BiAffine_;
typedef Tag<BiAffine_> const BiAffine;

template <typename TScoreValue>
class Score<TScoreValue, BiAffine>
{
    public:
    TScoreValue _match;
    TScoreValue _mismatch;
    TScoreValue _gapExtendHorizontal;
    TScoreValue _gapOpenHorizontal;
    TScoreValue _gapExtendVertical;
    TScoreValue _gapOpenVertical;


    Score() : _match(0),
              _mismatch(0),
              _gapExtendHorizontal(0),
              _gapOpenHorizontal(0),
              _gapExtendVertical(0),
              _gapOpenVertical(0) {}

    Score(TScoreValue match, TScoreValue mismatch, TScoreValue gap) :
                                    _match(match),
                                    _mismatch(mismatch),
                                    _gapExtendHorizontal(gap),
                                    _gapOpenHorizontal(gap),
                                    _gapExtendVertical(gap),
                                    _gapOpenVertical(gap) {}

    Score(TScoreValue match, TScoreValue mismatch, TScoreValue gapExtend, TScoreValue gapOpen) :
                                    _match(match),
                                    _mismatch(mismatch),
                                    _gapExtendHorizontal(gapExtend),
                                    _gapOpenHorizontal(gapOpen),
                                    _gapExtendVertical(gapExtend),
                                    _gapOpenVertical(gapOpen) {}

    Score(TScoreValue match, TScoreValue mismatch, TScoreValue gapExtendHorizontal, TScoreValue gapOpenHorizontal,
          TScoreValue gapExtendVertical, TScoreValue gapOpenVertical) : _match(match),
                                                                        _mismatch(mismatch),
                                                                        _gapExtendHorizontal(gapExtendHorizontal),
                                                                        _gapOpenHorizontal(gapOpenHorizontal),
                                                                        _gapExtendVertical(gapExtendVertical),
                                                                        _gapOpenVertical(gapOpenVertical) {}

    Score(Score const & other) : _match(other._match),
                                 _mismatch(other._mismatch),
                                 _gapExtendHorizontal(other._gapExtendHorizontal),
                                 _gapOpenHorizontal(other._gapOpenHorizontal),
                                 _gapExtendVertical(other._gapExtendVertical),
                                 _gapOpenVertical(other._gapOpenVertical) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function score()
// ----------------------------------------------------------------------------


template <typename TScoreValue, typename TSeqEntry1, typename TSeqEntry2>
inline TScoreValue
score(Score<TScoreValue, BiAffine> const & me,
      TSeqEntry1 const & seqEntry1,
      TSeqEntry2 const & seqEntry2)
{
    if (seqEntry1 == static_cast<TSeqEntry1>(seqEntry2))
        return me._match;
    return me._mismatch;
}

// ----------------------------------------------------------------------------
// Function scoreGapExtendHorizontal()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSeqEntry1, typename TSeqEntry2>
inline TScoreValue
scoreGapExtendHorizontal(Score<TScoreValue, BiAffine> const & me,
                         TSeqEntry1 const & /*seqEntry1*/,
                         TSeqEntry2 const & /*seqEntry2*/)
{
    return me._gapExtendHorizontal;
}

// ----------------------------------------------------------------------------
// Function scoreGapOpenHorizontal()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSeqEntry1, typename TSeqEntry2>
inline TScoreValue
scoreGapOpenHorizontal(Score<TScoreValue, BiAffine> const & me,
                       TSeqEntry1 const & /*seqEntry1*/,
                       TSeqEntry2 const & /*seqEntry2*/)
{
    return me._gapOpenHorizontal;
}

// ----------------------------------------------------------------------------
// Function scoreGapExtendVertical()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSeqEntry1, typename TSeqEntry2>
inline TScoreValue
scoreGapExtendVertical(Score<TScoreValue, BiAffine> const & me,
                       TSeqEntry1 const & /*seqEntry1*/,
                       TSeqEntry2 const & /*seqEntry2*/)
{
    return me._gapExtendVertical;
}

// ----------------------------------------------------------------------------
// Function scoreGapOpenVertical()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TSeqEntry1, typename TSeqEntry2>
inline TScoreValue
scoreGapOpenVertical(Score<TScoreValue, BiAffine> const & me,
                     TSeqEntry1 const & /*seqEntry1*/,
                     TSeqEntry2 const & /*seqEntry2*/)
{
    return me._gapOpenVertical;
}

// ----------------------------------------------------------------------------
// Function setScoreMatch()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline void
setScoreMatch(Score<TScoreValue, BiAffine> & scoringScheme, TScoreValue const & score)
{
    scoringScheme._match = score;
}

// ----------------------------------------------------------------------------
// Function setScoreMismatch()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline void
setScoreMismatch(Score<TScoreValue, BiAffine>  & scoringScheme, TScoreValue const & score)
{
    scoringScheme._mismatch = score;
}

// ----------------------------------------------------------------------------
// Function setScoreGap()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline void
setScoreGap(Score<TScoreValue, BiAffine> & scoringScheme, TScoreValue const & score)
{
    scoringScheme._gapOpenHorizontal = score;
    scoringScheme._gapOpenVertical = score;
    scoringScheme._gapExtendHorizontal = score;
    scoringScheme._gapExtendVertical = score;
}

// ----------------------------------------------------------------------------
// Function setScoreGapOpen()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline void
setScoreGapOpen(Score<TScoreValue, BiAffine> & scoringScheme, TScoreValue const & score)
{
    scoringScheme._gapOpenHorizontal = score;
    scoringScheme._gapOpenVertical = score;
}

// ----------------------------------------------------------------------------
// Function setScoreGapExtend()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline void
setScoreGapExtend(Score<TScoreValue, BiAffine> & scoringScheme, TScoreValue const & score)
{
    scoringScheme._gapExtendHorizontal = score;
    scoringScheme._gapExtendVertical = score;
}

// ----------------------------------------------------------------------------
// Function setScoreGapOpenHorizontal()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline void
setScoreGapOpenHorizontal(Score<TScoreValue, BiAffine> & scoringScheme, TScoreValue const & score)
{
    scoringScheme._gapOpenHorizontal = score;
}

// ----------------------------------------------------------------------------
// Function setScoreGapExtendHorizontal()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline void
setScoreGapExtendHorizontal(Score<TScoreValue, BiAffine> & scoringScheme, TScoreValue const & score)
{
    scoringScheme._gapExtendHorizontal = score;
}

// ----------------------------------------------------------------------------
// Function setScoreGapOpenVertical()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline void
setScoreGapOpenVertical(Score<TScoreValue, BiAffine> & scoringScheme, TScoreValue const & score)
{
    scoringScheme._gapOpenVertical = score;
}

// ----------------------------------------------------------------------------
// Function setScoreGapExtendVertical()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
inline void
setScoreGapExtendVertical(Score<TScoreValue, BiAffine> & scoringScheme, TScoreValue const & score)
{
    scoringScheme._gapExtendVertical = score;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_SET_SBIAFFINE_H_
