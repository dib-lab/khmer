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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Test for globalAlignmentScore() implementations that use Hirschberg and
// MyersBitVector, MyersHirschberg.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_SPECIALIZED_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_SPECIALIZED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function globalAlignment()                                      [Hirschberg]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            Hirschberg const & algorithmTag)
{
    SEQAN_ASSERT_EQ(length(rows(align)), 2u);
    return _globalAlignment(row(align, 0), row(align, 1), scoringScheme, algorithmTag);
}

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            Hirschberg const & algorithmTag)
{
    return _globalAlignment(gapsH, gapsV, scoringScheme, algorithmTag);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                                [Myers-Hirschberg]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec>
int globalAlignment(Align<TSequence, TAlignSpec> & align,
                    MyersHirschberg const & algorithmTag)
{
    SEQAN_ASSERT_EQ(length(rows(align)), 2u);
    return _globalAlignment(row(align, 0), row(align, 1), algorithmTag);
}

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV>
int globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                    Gaps<TSequenceV, TGapsSpecV> & gapsV,
                    MyersHirschberg const & algorithmTag)
{
    return _globalAlignment(gapsH, gapsV, algorithmTag);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                                  [MyersBitVector]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec>
int globalAlignment(Align<TSequence, TAlignSpec> & align,
                    MyersBitVector const & algorithmTag)
{
    SEQAN_ASSERT_EQ(length(rows(align)), 2u);
    return _globalAlignment(row(align, 0), row(align, 1), algorithmTag);
}

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV>
int globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                    Gaps<TSequenceV, TGapsSpecV> & gapsV,
                    MyersBitVector const & algorithmTag)
{
    return _globalAlignment(gapsH, gapsV, algorithmTag);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()                                 [Hirschberg]
// ----------------------------------------------------------------------------

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                                 String<TAlphabetV, TSpecV> const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 Hirschberg const & algorithmTag)
{
    Gaps<String<TAlphabetH, TSpecH> const, ArrayGaps> gapsH(seqH);
    Gaps<String<TAlphabetV, TSpecV> const, ArrayGaps> gapsV(seqV);
    return globalAlignment(gapsH, gapsV, scoringScheme, algorithmTag);
}

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 Hirschberg const & algorithmTag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);
    return globalAlignmentScore(strings[0], strings[1], scoringScheme, algorithmTag);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()                           [Myers-Hirschberg]
// ----------------------------------------------------------------------------

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV>
int globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                         String<TAlphabetV, TSpecV> const & seqV,
                         MyersHirschberg const & algorithmTag)
{
    Gaps<String<TAlphabetH, TSpecH> const, ArrayGaps> gapsH(seqH);
    Gaps<String<TAlphabetV, TSpecV> const, ArrayGaps> gapsV(seqV);
    return _globalAlignment(gapsH, gapsV, algorithmTag);
}

template <typename TString, typename TSpec>
int globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                         MyersHirschberg const & algorithmTag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);
    return globalAlignmentScore(strings[0], strings[1], algorithmTag);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()                             [MyersBitVector]
// ----------------------------------------------------------------------------

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV>
int globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                         String<TAlphabetV, TSpecV> const & seqV,
                         MyersBitVector const & algorithmTag)
{
    return _globalAlignmentScore(seqH, seqV, algorithmTag);
}

template <typename TString, typename TSpec>
int globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                         MyersBitVector const & algorithmTag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);
    return _globalAlignmentScore(strings[0], strings[1], algorithmTag);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_SPECIALIZED_H_
