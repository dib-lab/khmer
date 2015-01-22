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
// Global alignment interface for the banded Needleman-Wunsch and Gotoh
// algorithms.
//
// We define the interface functions pretty explicitely (versus just TAlign,
// TFragments etc.) so the candidates the compiler gives when resolution to
// the globalFunction() fails is actually meaningful.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_BANDED_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_BANDED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TScoreValue, typename TSpec>
class Score;
template <typename TSpec>
class Graph;
template <typename TStringSet, typename TCargo, typename TGraphSpec>
struct Alignment;
template <typename TSize, typename TFragmentSpec>
class Fragment;

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
// Function globalAlignment()                                   [banded, Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & algoTag)
{
    typedef Align<TSequence, TAlignSpec> TAlign;
    typedef typename Size<TAlign>::Type  TSize;
    typedef typename Position<TAlign>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> trace;

    // We do not need string ids for this variant and set them to 0u.  They are
    // only required for the Fragment String and the Alignment Graph variant.
    TScoreValue res = _setUpAndRunAlignment(trace, source(row(align, 0)), source(row(align, 1)), scoringScheme,
                                            alignConfig, lowerDiag, upperDiag, algoTag);
    _adaptTraceSegmentsTo(row(align, 0), row(align, 1), trace);
    return res;
}

// Interface without AlignConfig<>.

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(align, scoringScheme, alignConfig, lowerDiag, upperDiag, algoTag);
}

// Interface without algorithm tag.

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            int lowerDiag,
                            int upperDiag)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(align, scoringScheme, alignConfig, lowerDiag, upperDiag, NeedlemanWunsch());
    else
        return globalAlignment(align, scoringScheme, alignConfig, lowerDiag, upperDiag, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(align, scoringScheme, alignConfig, lowerDiag, upperDiag);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                                    [banded, Gaps]
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & algoTag)
{
    typedef typename Size<TSequenceH>::Type TSize;
    typedef typename Position<TSequenceH>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;

    // We do not need string ids for this variant and set them to 0u.  They are
    // only required for the Fragment String and the Alignment Graph variant.
    TScoreValue res = _setUpAndRunAlignment(traceSegments, source(gapsH), source(gapsV), scoringScheme, alignConfig,
                                            lowerDiag, upperDiag, algoTag);
    _adaptTraceSegmentsTo(gapsH, gapsV, traceSegments);
    return res;
}

// Interface without AlignConfig<>.

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, lowerDiag, upperDiag, algoTag);
}

// Interface without algorithm tag.

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            int lowerDiag,
                            int upperDiag)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, lowerDiag, upperDiag, NeedlemanWunsch());
    else
        return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, lowerDiag, upperDiag, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, lowerDiag, upperDiag);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                     [banded, Graph<Alignment<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & algoTag)
{
    typedef Graph<Alignment<TStringSet, TCargo, TGraphSpec> > TGraph;
    typedef typename Position<TGraph>::Type TPosition;
    typedef typename Size<TGraph>::Type TSize;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;

    TScoreValue res = _setUpAndRunAlignment(traceSegments, value(stringSet(alignmentGraph), 0),
                                            value(stringSet(alignmentGraph), 1), scoringScheme, alignConfig, lowerDiag,
                                            upperDiag, algoTag);
    _adaptTraceSegmentsTo(alignmentGraph, positionToId(stringSet(alignmentGraph), 0),
                          positionToId(stringSet(alignmentGraph), 1), traceSegments);
    return res;
}

// Interface without AlignConfig<>.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(alignmentGraph, scoringScheme, alignConfig, lowerDiag, upperDiag, algoTag);
}

// Interface without algorithm tag.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            int lowerDiag,
                            int upperDiag)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(alignmentGraph, scoringScheme, alignConfig, lowerDiag, upperDiag, NeedlemanWunsch());
    else
        return globalAlignment(alignmentGraph, scoringScheme, alignConfig, lowerDiag, upperDiag, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(alignmentGraph, scoringScheme, alignConfig, lowerDiag, upperDiag);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                     [banded, String<Fragment<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & algoTag)
{
    typedef String<Fragment<TSize, TFragmentSpec>, TStringSpec> TFragments;
    typedef typename Position<TFragments>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;

    TScoreValue res = _setUpAndRunAlignment(traceSegments, value(strings, 0), value(strings, 1), scoringScheme,
                                            alignConfig, lowerDiag, upperDiag, algoTag);
    _adaptTraceSegmentsTo(fragmentString, positionToId(strings, 0), positionToId(strings, 1), traceSegments);
    return res;
}

// Interface without AlignConfig<>.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, lowerDiag, upperDiag, algoTag);
}

// Interface without algorithm tag.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            int lowerDiag,
                            int upperDiag)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, lowerDiag, upperDiag, NeedlemanWunsch());
    else
        return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, lowerDiag, upperDiag, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, lowerDiag, upperDiag);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()                        [banded, 2 Strings]
// ----------------------------------------------------------------------------

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                                 String<TAlphabetV, TSpecV> const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                                 int lowerDiag,
                                 int upperDiag,
                                 TAlgoTag const & algoTag)
{
    return _setUpAndRunAlignment(seqH, seqV, scoringScheme, alignConfig, lowerDiag, upperDiag, algoTag);
}

// Interface without AlignConfig<>.

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                                 String<TAlphabetV, TSpecV> const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 int lowerDiag,
                                 int upperDiag,
                                 TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, lowerDiag, upperDiag, algoTag);
}

// Interface without algorithm tag.

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                                 String<TAlphabetV, TSpecV> const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                                 int lowerDiag,
                                 int upperDiag)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, lowerDiag, upperDiag, NeedlemanWunsch());
    else
        return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, lowerDiag, upperDiag, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                                 String<TAlphabetV, TSpecV> const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 int lowerDiag,
                                 int upperDiag)
{
    AlignConfig<> alignConfig;
    return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, lowerDiag, upperDiag);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()                        [banded, StringSet]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                                 int lowerDiag,
                                 int upperDiag,
                                 TAlgoTag const & algoTag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);
    return _setUpAndRunAlignment(strings[0], strings[1], scoringScheme, alignConfig, lowerDiag, upperDiag, algoTag);
}

// Interface without AlignConfig<>.

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 int lowerDiag,
                                 int upperDiag,
                                 TAlgoTag const & algoTag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);

    AlignConfig<> alignConfig;
    return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, lowerDiag, upperDiag, algoTag);
}

// Interface without algorithm tag.

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                                 int lowerDiag,
                                 int upperDiag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);

    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, lowerDiag, upperDiag, NeedlemanWunsch());
    else
        return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, lowerDiag, upperDiag, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 int lowerDiag,
                                 int upperDiag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);

    AlignConfig<> alignConfig;
    return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, lowerDiag, upperDiag);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_BANDED_H_
