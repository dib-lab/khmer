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
// Global alignment interface for journaled strings.
// ==========================================================================

#ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNAL_ALIGNMENT_INTERFACE_H_
#define INCLUDE_SEQAN_JOURNALED_SET_JOURNAL_ALIGNMENT_INTERFACE_H_

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
// Function globalAlignment()                  [unbanded, String<Journaled<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TValue, typename THostSpec, typename TJournalSpec,
          typename TBuffSpec, typename TReference, typename TSource,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journaledString,
                            TReference const & reference,
                            TSource const & source,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                            TAlgoTag const & /*tag*/)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Position<TJournalString>::Type TPosition;
    typedef typename Size<TJournalString>::Type TSize;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPConfig;

    String<TTraceSegment> traceSegments;
    DPScoutState_<Default> dpScoutState;
    // We need to do that in order to build journal strings from two simple sequences.
    TScoreValue res = _setUpAndRunAlignment(traceSegments, dpScoutState, reference, source, scoringScheme, TDPConfig(),
                                            IsSameType<TAlgoTag, NeedlemanWunsch>());
    _adaptTraceSegmentsTo(journaledString, reference, source, traceSegments);
    return res;
}

// Interface without AlignConfig<>.

template <typename TValue, typename THostSpec, typename TJournalSpec,
          typename TBuffSpec, typename TReference, typename TSource,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journaledString,
                            TReference const & reference,
                            TSource const & source,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(journaledString, reference, source, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec,
          typename TReference, typename TSource,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journaledString,
                            TReference const & reference,
                            TSource const & source,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    typedef typename SequenceEntryForScore<Score<TScoreValue, TScoreSpec>, TReference>::Type TReferenceEntryForScore;
    typedef typename SequenceEntryForScore<Score<TScoreValue, TScoreSpec>, TSource>::Type TSourceEntryForScore;

    TReferenceEntryForScore entryRef = sequenceEntryForScore(scoringScheme, reference, 0);
    TSourceEntryForScore entrySource = sequenceEntryForScore(scoringScheme, source, 0);
    if (scoreGapOpenHorizontal(scoringScheme, entryRef, entrySource) ==
        scoreGapExtendHorizontal(scoringScheme, entryRef, entrySource) ||
        scoreGapOpenVertical(scoringScheme, entryRef, entrySource) ==
        scoreGapExtendVertical(scoringScheme, entryRef, entrySource))
        return globalAlignment(journaledString, reference, source, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(journaledString, reference, source, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TValue, typename THostSpec, typename TJournalSpec,
          typename TBuffSpec, typename TReference, typename TSource,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journaledString,
                            TReference const & reference,
                            TSource const & source,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(journaledString, reference, source, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                    [banded, String<Journaled<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TValue, typename THostSpec, typename TJournalSpec,
          typename TBuffSpec, typename TReference, typename TSource,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journaledString,
                            TReference const & reference,
                            TSource const & source,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & /*algoTag*/)
{
    typedef String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > TJournalString;
    typedef typename Position<TJournalString>::Type TPosition;
    typedef typename Size<TJournalString>::Type TSize;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOn>, TFreeEndGaps, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > TDPConfig;

    String<TTraceSegment> traceSegments;
    DPScoutState_<Default> dpScoutState;
    // We need to do that in order to build journal strings from two simple sequences.
    TScoreValue res = _setUpAndRunAlignment(traceSegments, dpScoutState, reference, source, scoringScheme,
                                            TDPConfig(lowerDiag, upperDiag),
                                            IsSameType<TAlgoTag, NeedlemanWunsch>());
    _adaptTraceSegmentsTo(journaledString, reference, source, traceSegments);
    return res;
}

// Interface without AlignConfig<>.

template <typename TValue, typename THostSpec, typename TJournalSpec,
          typename TBuffSpec, typename TReference, typename TSource,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journaledString,
                            TReference const & reference,
                            TSource const & source,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(journaledString, reference, source, scoringScheme, alignConfig, lowerDiag, upperDiag,
                           algoTag);
}

// Interface without algorithm tag.

template <typename TValue, typename THostSpec, typename TJournalSpec, typename TBuffSpec,
          typename TReference, typename TSource,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journaledString,
                            TReference const & reference,
                            TSource const & source,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            int lowerDiag,
                            int upperDiag)
{
    typedef typename SequenceEntryForScore<Score<TScoreValue, TScoreSpec>, TReference>::Type TReferenceEntryForScore;
    typedef typename SequenceEntryForScore<Score<TScoreValue, TScoreSpec>, TSource>::Type TSourceEntryForScore;

    TReferenceEntryForScore entryRef = sequenceEntryForScore(scoringScheme, reference, 0);
    TSourceEntryForScore entrySource = sequenceEntryForScore(scoringScheme, source, 0);
    if (scoreGapOpenHorizontal(scoringScheme, entryRef, entrySource) ==
        scoreGapExtendHorizontal(scoringScheme, entryRef, entrySource) ||
        scoreGapOpenVertical(scoringScheme, entryRef, entrySource) ==
        scoreGapExtendVertical(scoringScheme, entryRef, entrySource))
        return globalAlignment(journaledString, reference, source, scoringScheme, alignConfig, lowerDiag, upperDiag,
                               NeedlemanWunsch());
    else
        return globalAlignment(journaledString, reference, source, scoringScheme, alignConfig, lowerDiag, upperDiag,
                               Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TValue, typename THostSpec, typename TJournalSpec,
          typename TBuffSpec, typename TReference, typename TSource,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(String<TValue, Journaled<THostSpec, TJournalSpec, TBuffSpec> > & journaledString,
                            TReference const & reference,
                            TSource const & source,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            int lowerDiag,
                            int upperDiag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(journaledString, reference, source, scoringScheme, alignConfig, lowerDiag, upperDiag);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_JOURNALED_SET_JOURNAL_ALIGNMENT_INTERFACE_H_
