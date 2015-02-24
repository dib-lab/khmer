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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_UNBANDED_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_UNBANDED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Specialization Unbanded LocalAlignmentEnumerator
// ----------------------------------------------------------------------------

struct Unbanded_;
typedef Tag<Unbanded_> Unbanded;

template <typename TScore>
class LocalAlignmentEnumerator<TScore, Unbanded>
{
public:
    typedef typename Value<TScore>::Type TScoreValue_;

    TScore _scoringScheme;
    TScoreValue_ _cutoff;
    LocalAlignmentFinder<TScoreValue_> _finder;

    LocalAlignmentEnumerator(TScore const & score) :
            _scoringScheme(score), _cutoff(0)
    {}

    LocalAlignmentEnumerator(TScore const & score, int cutoff) :
            _scoringScheme(score), _cutoff(cutoff)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getScore()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TScoreSpec>
inline TScoreValue
getScore(LocalAlignmentEnumerator<Score<TScoreValue, TScoreSpec>, Unbanded> const & enumerator)
{
    // TODO(holtgrew): Const correctness!
    return getScore(const_cast<LocalAlignmentFinder<TScoreValue> &>(enumerator._finder));
}

// ----------------------------------------------------------------------------
// Function nextLocalAlignment()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TScoreSpec, typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV>
inline bool
nextLocalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                   Gaps<TSequenceV, TGapsSpecV> & gapsV,
                   LocalAlignmentEnumerator<Score<TScoreValue, TScoreSpec>, Unbanded> & enumerator)
{
    if (enumerator._finder.needReinit)
        return _smithWaterman(gapsH, gapsV, enumerator._finder, enumerator._scoringScheme, enumerator._cutoff);
    else
        return _smithWatermanGetNext(gapsH, gapsV, enumerator._finder, enumerator._scoringScheme, enumerator._cutoff);
}

template <typename TScoreValue, typename TScoreSpec, typename TSequence, typename TAlignSpec>
inline bool
nextLocalAlignment(Align<TSequence, TAlignSpec> & align,
                   LocalAlignmentEnumerator<Score<TScoreValue, TScoreSpec>, Unbanded> & enumerator)
{
    SEQAN_ASSERT_EQ(length(rows(align)), 2u);

    return nextLocalAlignment(row(align, 0), row(align, 1), enumerator);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_UNBANDED_H_
