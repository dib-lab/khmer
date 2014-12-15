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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_SETUP_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_SETUP_H_

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

// ----------------------------------------------------------------------------
// SubstituteAlignConfig_
// ----------------------------------------------------------------------------

template <typename TAlignConfig>
struct SubstituteAlignConfig_;

// 0000
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, false, false, false, TSpec> >
{
    typedef FreeEndGaps_<False, False, False, False> Type;
};

// 0001
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, false, false, true, TSpec> >
{
    typedef FreeEndGaps_<False, False, True, False> Type;
};

// 0010
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, false, true, false, TSpec> >
{
    typedef FreeEndGaps_<False, False, False, True> Type;
};


// 0011
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, false, true, true, TSpec> >
{
    typedef FreeEndGaps_<False, False, True, True> Type;
};


// 0100
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, true, false, false, TSpec> >
{
    typedef FreeEndGaps_<False, True, False, False> Type;
};


// 0101
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, true, false, true, TSpec> >
{
    typedef FreeEndGaps_<False, True, True, False> Type;
};


// 0110
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, true, true, false, TSpec> >
{
    typedef FreeEndGaps_<False, True, False, True> Type;
};


// 0111
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<false, true, true, true, TSpec> >
{
    typedef FreeEndGaps_<False, True, True, True> Type;
};


// 1000
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, false, false, false, TSpec> >
{
    typedef FreeEndGaps_<True, False, False, False> Type;
};


// 1001
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, false, false, true, TSpec> >
{
    typedef FreeEndGaps_<True, False, True, False> Type;
};


// 1010
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, false, true, false, TSpec> >
{
    typedef FreeEndGaps_<True, False, False, True> Type;
};


// 1011
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, false, true, true, TSpec> >
{
    typedef FreeEndGaps_<True, False, True, True> Type;
};


// 1100
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, true, false, false, TSpec> >
{
    typedef FreeEndGaps_<True, True, False, False> Type;
};


// 1101
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, true, false, true, TSpec> >
{
    typedef FreeEndGaps_<True, True, True, False> Type;
};


// 1110
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, true, true, false, TSpec> >
{
    typedef FreeEndGaps_<True, True, False, True> Type;
};

// 1111
template <typename TSpec>
struct SubstituteAlignConfig_<AlignConfig<true, true, true, true, TSpec> >
{
    typedef FreeEndGaps_<True, True, True, True> Type;
};

// ----------------------------------------------------------------------------
// SetUpAlignmentProfile
// ----------------------------------------------------------------------------

template <typename TAlgoTag, typename TAlignConfig, typename TTraceSwitch, typename TGapCosts>
struct SetupAlignmentProfile_;

// Profile for Needleman-Wunsch algorithm.
template <typename TAlignConfig, typename TGapCosts, typename TTraceSwitch>
struct SetupAlignmentProfile_<NeedlemanWunsch, TAlignConfig, TGapCosts, TTraceSwitch>
{
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps_;
    typedef DPProfile_<GlobalAlignment_<TFreeEndGaps_>, LinearGaps, TTraceSwitch> Type;
};

// Profile for Gotoh algorithm.
template <typename TAlignConfig, typename TGapCosts, typename TTraceSwitch>
struct SetupAlignmentProfile_<Gotoh, TAlignConfig, TGapCosts, TTraceSwitch>
{
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps_;
    typedef DPProfile_<GlobalAlignment_<TFreeEndGaps_>, AffineGaps, TTraceSwitch> Type;
};

// Profile for Smith-Waterman algorithm.
template <typename TAlignConfig, typename TGapCosts, typename TTraceSwitch>
struct SetupAlignmentProfile_<SmithWaterman, TAlignConfig, TGapCosts, TTraceSwitch>
{
    typedef DPProfile_<LocalAlignment_<>, TGapCosts, TTraceSwitch> Type;
};

// Profile for Waterman-Eggert algorithm
template <typename TAlignConfig, typename TGapCosts, typename TTraceSwitch>
struct SetupAlignmentProfile_<WatermanEggert, TAlignConfig, TGapCosts, TTraceSwitch>
{
    typedef DPProfile_<LocalAlignment_<SuboptimalAlignment>, TGapCosts, TracebackOn<TracebackConfig_<SingleTrace, GapsLeft> > > Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _setUpAndRunAlignment()                                  [Unbanded]
// ----------------------------------------------------------------------------

// Interface without AlignConfig.
template <typename TDPScoutStateSpec, typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV,
          typename TScoreValue, typename TScoreSpec, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlgoTag const &,
                      TGapsTag const &)
{
    typedef Score<TScoreValue, TScoreSpec> TScoringScheme;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceH>::Type TSequenceHEntry;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceV>::Type TSequenceVEntry;

    SEQAN_ASSERT_GEQ(length(seqH), 1u);
    SEQAN_ASSERT_GEQ(length(seqV), 1u);

    TSequenceHEntry seqHEntry = sequenceEntryForScore(scoringScheme, seqH, 0);
    TSequenceVEntry seqVEntry = sequenceEntryForScore(scoringScheme, seqV, 0);

    if (scoreGapExtendHorizontal(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenHorizontal(scoringScheme, seqHEntry, seqVEntry) ||
        scoreGapExtendVertical(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenVertical(scoringScheme, seqHEntry, seqVEntry))
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, AlignConfig<>, AffineGaps, TracebackOn<TGapsTag> >::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOff>(), TDPProfile());
    }
    else
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, AlignConfig<>, LinearGaps, TracebackOn<TGapsTag> >::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOff>(), TDPProfile());
    }
}

template <typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV, typename TScoreValue,
          typename TScoreSpec, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlgoTag const & algoTag,
                      TGapsTag const & gapsTag)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(traceSegments, noState, seqH, seqV, scoringScheme, algoTag, gapsTag);
}

template <typename TDPScoutStateSpec, typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV,
          typename TScoreValue, typename TScoreSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlgoTag const & algoTag)
{
    return _setUpAndRunAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

template <typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV, typename TScoreValue,
          typename TScoreSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlgoTag const & algoTag)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(traceSegments, noState, seqH, seqV, scoringScheme, algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

// Interface with AlignConfig.
template <typename TDPScoutStateSpec, typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV,
          typename TScoreValue, typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom, typename TACSpec,
          typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const &,
                      TAlgoTag const &,
                      TGapsTag const &)
{
    typedef Score<TScoreValue, TScoreSpec> TScoringScheme;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceH>::Type TSequenceHEntry;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceV>::Type TSequenceVEntry;
    typedef AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> TAlignConfig;

    SEQAN_ASSERT_GEQ(length(seqH), 1u);
    SEQAN_ASSERT_GEQ(length(seqV), 1u);

    TSequenceHEntry seqHEntry = sequenceEntryForScore(scoringScheme, seqH, 0);
    TSequenceVEntry seqVEntry = sequenceEntryForScore(scoringScheme, seqV, 0);

    if (scoreGapExtendHorizontal(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenHorizontal(scoringScheme, seqHEntry, seqVEntry) ||
        scoreGapExtendVertical(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenVertical(scoringScheme, seqHEntry, seqVEntry))
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, TAlignConfig, AffineGaps, TracebackOn<TGapsTag> >::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOff>(), TDPProfile());
    }
    else
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, TAlignConfig, LinearGaps, TracebackOn<TGapsTag> >::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOff>(), TDPProfile());
    }
}

template <typename TDPScoutStateSpec, typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV,
          typename TScoreValue, typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom, typename TACSpec,
          typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const & alignConfig,
                      TAlgoTag const & algoTag)
{
    return _setUpAndRunAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, alignConfig, algoTag,
                                 TracebackConfig_<SingleTrace, GapsLeft>());
}

template <typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV, typename TScoreValue,
          typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom, typename TACSpec,
          typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const & alignConfig,
                      TAlgoTag const & algoTag,
                      TGapsTag const & gapsTag)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(traceSegments, noState, seqH, seqV, scoringScheme, alignConfig, algoTag, gapsTag);
}

template <typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV, typename TScoreValue,
          typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom, typename TACSpec,
          typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const & alignConfig,
                      TAlgoTag const & algoTag)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(traceSegments, noState, seqH, seqV, scoringScheme, alignConfig, algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

// Interface without AlignConfig and with traceback disabled.
template <typename TDPScoutStateSpec, typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlgoTag const &,
                      TGapsTag const & /*unused*/)
{
    typedef String<TAlphabetH, TSpecH> const TSequenceH;
    typedef String<TAlphabetV, TSpecV> const TSequenceV;
    typedef Score<TScoreValue, TScoreSpec> TScoringScheme;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceH>::Type TSequenceHEntry;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceV>::Type TSequenceVEntry;

    SEQAN_ASSERT_GEQ(length(seqH), 1u);
    SEQAN_ASSERT_GEQ(length(seqV), 1u);

    String<TraceSegment_<unsigned, unsigned> > traceSegments;

    TSequenceHEntry seqHEntry = sequenceEntryForScore(scoringScheme, seqH, 0);
    TSequenceVEntry seqVEntry = sequenceEntryForScore(scoringScheme, seqV, 0);

    if (scoreGapExtendHorizontal(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenHorizontal(scoringScheme, seqHEntry, seqVEntry) ||
        scoreGapExtendVertical(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenVertical(scoringScheme, seqHEntry, seqVEntry))
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, AlignConfig<>, AffineGaps, TracebackOff>::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOff>(), TDPProfile());
    }
    else
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, AlignConfig<>, LinearGaps, TracebackOff>::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOff>(), TDPProfile());
    }
}

template <typename TDPScoutStateSpec, typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlgoTag const & algoTag)
{
    return _setUpAndRunAlignment(dpScoutState, seqH, seqV, scoringScheme, algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

template <typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV, typename TScoreValue,
          typename TScoreSpec, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlgoTag const &,
                      TGapsTag const & gapsTag)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(noState, seqH, seqV, scoringScheme, gapsTag);
}

template <typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV, typename TScoreValue,
          typename TScoreSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlgoTag const &)
{
    // Note that GapsLeft could be nothing, is unused in callee without traceback.
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(noState, seqH, seqV, scoringScheme, TracebackConfig_<SingleTrace, GapsLeft>());
}

// Interface with AlignConfig and with traceback disabled.
template <typename TDPScoutStateSpec, typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom,
          typename TACSpec, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const &,
                      TAlgoTag const &,
                      TGapsTag const & /*unused*/)
{
    typedef AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> TAlignConfig;
    typedef String<TAlphabetH, TSpecH> const TSequenceH;
    typedef String<TAlphabetV, TSpecV> const TSequenceV;
    typedef Score<TScoreValue, TScoreSpec> TScoringScheme;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceH>::Type TSequenceHEntry;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceV>::Type TSequenceVEntry;

    SEQAN_ASSERT_GEQ(length(seqH), 1u);
    SEQAN_ASSERT_GEQ(length(seqV), 1u);

    String<TraceSegment_<unsigned, unsigned> > traceSegments;

    TSequenceHEntry seqHEntry = sequenceEntryForScore(scoringScheme, seqH, 0);
    TSequenceVEntry seqVEntry = sequenceEntryForScore(scoringScheme, seqV, 0);

    if (scoreGapExtendHorizontal(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenHorizontal(scoringScheme, seqHEntry, seqVEntry) ||
        scoreGapExtendVertical(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenVertical(scoringScheme, seqHEntry, seqVEntry))
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, TAlignConfig, AffineGaps, TracebackOff>::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOff>(), TDPProfile());
    }
    else
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, TAlignConfig, LinearGaps, TracebackOff>::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOff>(), TDPProfile());
    }
}

template <typename TDPScoutStateSpec, typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom,
          typename TACSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const & alignConfig,
                      TAlgoTag const &)
{
    // Note that GapsLeft could be nothing, is unused in callee without traceback.
    return _setUpAndRunAlignment(dpScoutState, seqH, seqV, scoringScheme, alignConfig, TracebackConfig_<SingleTrace, GapsLeft>());
}

template <typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV, typename TScoreValue,
          typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom, typename TACSpec, typename TAlgoTag,
          typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const & alignConfig,
                      TAlgoTag const & algoTag,
                      TGapsTag const & gapsTag)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(noState, seqH, seqV, scoringScheme, alignConfig, algoTag, gapsTag);
}

template <typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV, typename TScoreValue,
          typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom, typename TACSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const & alignConfig,
                      TAlgoTag const & algoTag)
{
    DPScoutState_<Default> noState;
    // Note that GapsLeft could be nothing, is unused in callee without traceback.
    return _setUpAndRunAlignment(noState, seqH, seqV, scoringScheme, alignConfig, algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

// ----------------------------------------------------------------------------
// Function _setUpAndRunAlignment()                                    [Banded]
// ----------------------------------------------------------------------------

// Interface without AlignConfig.
template <typename TDPScoutStateSpec, typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV,
          typename TScoreValue, typename TScoreSpec, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const &,
                      TGapsTag const &)
{
    typedef Score<TScoreValue, TScoreSpec> TScoringScheme;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceH>::Type TSequenceHEntry;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceV>::Type TSequenceVEntry;
    
    SEQAN_ASSERT_GEQ(length(seqH), 1u);
    SEQAN_ASSERT_GEQ(length(seqV), 1u);

    TSequenceHEntry seqHEntry = sequenceEntryForScore(scoringScheme, seqH, 0);
    TSequenceVEntry seqVEntry = sequenceEntryForScore(scoringScheme, seqV, 0);

    if (scoreGapExtendHorizontal(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenHorizontal(scoringScheme, seqHEntry, seqVEntry) ||
        scoreGapExtendVertical(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenVertical(scoringScheme, seqHEntry, seqVEntry))
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, AlignConfig<>, AffineGaps, TracebackOn<TGapsTag> >::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOn>(lowerDiagonal, upperDiagonal), TDPProfile());
    }
    else
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, AlignConfig<>, LinearGaps, TracebackOn<TGapsTag> >::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOn>(lowerDiagonal, upperDiagonal), TDPProfile());
    }
}

template <typename TDPScoutStateSpec, typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV,
          typename TScoreValue, typename TScoreSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag)
{
    // Note that GapsLeft could be nothing, is unused in callee without traceback.
    return _setUpAndRunAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, lowerDiagonal, upperDiagonal,
                                 algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

template <typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV, typename TScoreValue,
          typename TScoreSpec, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag,
                      TGapsTag const & /*ignored*/)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(traceSegments, noState, seqH, seqV, scoringScheme, lowerDiagonal, upperDiagonal, algoTag);
}

template <typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV, typename TScoreValue,
          typename TScoreSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag)
{
    // Note that GapsLeft could be nothing, is unused in callee without traceback.
    return _setUpAndRunAlignment(traceSegments, seqH, seqV, scoringScheme, lowerDiagonal, upperDiagonal, algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

// Interface with AlignConfig.
template <typename TDPScoutStateSpec, typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV,
          typename TScoreValue, typename TScoreSpec, typename TAlignConfig, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlignConfig const &,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const &,
                      TGapsTag const &)
{
    typedef Score<TScoreValue, TScoreSpec> TScoringScheme;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceH>::Type TSequenceHEntry;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceV>::Type TSequenceVEntry;

    SEQAN_ASSERT_GEQ(length(seqH), 1u);
    SEQAN_ASSERT_GEQ(length(seqV), 1u);

    TSequenceHEntry seqHEntry = sequenceEntryForScore(scoringScheme, seqH, 0);
    TSequenceVEntry seqVEntry = sequenceEntryForScore(scoringScheme, seqV, 0);

    if (scoreGapExtendHorizontal(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenHorizontal(scoringScheme, seqHEntry, seqVEntry) ||
        scoreGapExtendVertical(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenVertical(scoringScheme, seqHEntry, seqVEntry))
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, TAlignConfig, AffineGaps, TracebackOn<TGapsTag> >::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOn>(lowerDiagonal, upperDiagonal), TDPProfile());
    }
    else
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, TAlignConfig, LinearGaps, TracebackOn<TGapsTag> >::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOn>(lowerDiagonal, upperDiagonal), TDPProfile());
    }
}

template <typename TDPScoutStateSpec, typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV,
          typename TScoreValue, typename TScoreSpec, typename TAlignConfig, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlignConfig const & alignConfig,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag)
{
    return _setUpAndRunAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, alignConfig,
                                 lowerDiagonal, upperDiagonal, algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}


template <typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV, typename TScoreValue,
          typename TScoreSpec, typename TAlignConfig, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlignConfig const & alignConfig,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag,
                      TGapsTag const & gapsTag)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(traceSegments, noState, seqH, seqV, scoringScheme, alignConfig, lowerDiagonal,
                                 upperDiagonal, algoTag, gapsTag);
}

template <typename TTraceSegment, typename TSpec, typename TSequenceH, typename TSequenceV, typename TScoreValue,
          typename TScoreSpec, typename TAlignConfig, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TTraceSegment, TSpec> & traceSegments,
                      TSequenceH const & seqH,
                      TSequenceV const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      TAlignConfig const & alignConfig,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag)
{
    return _setUpAndRunAlignment(traceSegments, seqH, seqV, scoringScheme, alignConfig, lowerDiagonal, upperDiagonal,
                                 algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

// Interface without AlignConfig and with traceback disabled.
template <typename TDPScoutStateSpec, typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const &,
                      TGapsTag const &)
{
    typedef String<TAlphabetH, TSpecH> const TSequenceH;
    typedef String<TAlphabetV, TSpecV> const TSequenceV;
    typedef Score<TScoreValue, TScoreSpec> TScoringScheme;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceH>::Type TSequenceHEntry;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceV>::Type TSequenceVEntry;

    SEQAN_ASSERT_GEQ(length(seqH), 1u);
    SEQAN_ASSERT_GEQ(length(seqV), 1u);

    String<TraceSegment_<unsigned, unsigned> > traceSegments;

    TSequenceHEntry seqHEntry = sequenceEntryForScore(scoringScheme, seqH, 0);
    TSequenceVEntry seqVEntry = sequenceEntryForScore(scoringScheme, seqV, 0);

    if (scoreGapExtendHorizontal(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenHorizontal(scoringScheme, seqHEntry, seqVEntry) ||
        scoreGapExtendVertical(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenVertical(scoringScheme, seqHEntry, seqVEntry))
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, AlignConfig<>, AffineGaps, TracebackOff>::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOn>(lowerDiagonal, upperDiagonal), TDPProfile());
    }
    else
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, AlignConfig<>, LinearGaps, TracebackOff>::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOn>(lowerDiagonal, upperDiagonal), TDPProfile());
    }
}

template <typename TDPScoutStateSpec, typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag)
{
    // Note that GapsLeft could be nothing, is unused in callee without traceback.
    return _setUpAndRunAlignment(dpScoutState, seqH, seqV, scoringScheme, lowerDiagonal, upperDiagonal, algoTag,
                                 TracebackConfig_<SingleTrace, GapsLeft>());
}

template <typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV, typename TScoreValue,
          typename TScoreSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(noState, seqH, seqV, scoringScheme, lowerDiagonal, upperDiagonal, algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

// Interface with AlignConfig and with traceback disabled.
template <typename TDPScoutStateSpec, typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom,
          typename TACSpec, typename TAlgoTag, typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const &,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const &,
                      TGapsTag const & /*ignored*/)
{
    typedef AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> TAlignConfig;
    typedef String<TAlphabetH, TSpecH> const TSequenceH;
    typedef String<TAlphabetV, TSpecV> const TSequenceV;
    typedef Score<TScoreValue, TScoreSpec> TScoringScheme;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceH>::Type TSequenceHEntry;
    typedef typename SequenceEntryForScore<TScoringScheme, TSequenceV>::Type TSequenceVEntry;

    SEQAN_ASSERT_GEQ(length(seqH), 1u);
    SEQAN_ASSERT_GEQ(length(seqV), 1u);

    String<TraceSegment_<unsigned, unsigned> > traceSegments;

    TSequenceHEntry seqHEntry = sequenceEntryForScore(scoringScheme, seqH, 0);
    TSequenceVEntry seqVEntry = sequenceEntryForScore(scoringScheme, seqV, 0);

    if (scoreGapExtendHorizontal(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenHorizontal(scoringScheme, seqHEntry, seqVEntry) ||
        scoreGapExtendVertical(scoringScheme, seqHEntry, seqVEntry) !=
        scoreGapOpenVertical(scoringScheme, seqHEntry, seqVEntry))
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, TAlignConfig, AffineGaps, TracebackOff>::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOn>(lowerDiagonal, upperDiagonal), TDPProfile());
    }
    else
    {
        typedef typename SetupAlignmentProfile_<TAlgoTag, TAlignConfig, LinearGaps, TracebackOff>::Type TDPProfile;
        return _computeAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, DPBand_<BandOn>(lowerDiagonal, upperDiagonal), TDPProfile());
    }
}

template <typename TDPScoutStateSpec, typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom,
          typename TACSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(DPScoutState_<TDPScoutStateSpec> & dpScoutState,
                      String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const & alignConfig,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag)
{
    // Note that GapsLeft could be nothing, is unused in callee without traceback.
    return _setUpAndRunAlignment(dpScoutState, seqH, seqV, scoringScheme, alignConfig, lowerDiagonal, upperDiagonal,
                                 algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

template <typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV, typename TScoreValue,
          typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom, typename TACSpec, typename TAlgoTag,
          typename TGapsTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const & alignConfig,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag,
                      TGapsTag const & /*ignored*/)
{
    DPScoutState_<Default> noState;
    return _setUpAndRunAlignment(noState, seqH, seqV, scoringScheme, alignConfig, lowerDiagonal, upperDiagonal,
                                 algoTag);
}

template <typename TAlphabetH, typename TSpecH, typename TAlphabetV, typename TSpecV, typename TScoreValue,
          typename TScoreSpec, bool TTop, bool TRight, bool TLeft, bool TBottom, typename TACSpec, typename TAlgoTag>
typename Value<Score<TScoreValue, TScoreSpec> >::Type
_setUpAndRunAlignment(String<TAlphabetH, TSpecH> const & seqH,
                      String<TAlphabetV, TSpecV> const & seqV,
                      Score<TScoreValue, TScoreSpec> const & scoringScheme,
                      AlignConfig<TTop, TRight, TLeft, TBottom, TACSpec> const & alignConfig,
                      int lowerDiagonal,
                      int upperDiagonal,
                      TAlgoTag const & algoTag)
{
    // Note that GapsLeft could be nothing, is unused in callee without traceback.
    return _setUpAndRunAlignment(seqH, seqV, scoringScheme, alignConfig, lowerDiagonal, upperDiagonal, algoTag, TracebackConfig_<SingleTrace, GapsLeft>());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_SETUP_H_
