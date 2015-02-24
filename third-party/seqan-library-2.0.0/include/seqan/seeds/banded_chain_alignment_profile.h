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
// The adaption of the meta profile and the setup for the banded chain
// alignment.
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_META_INFO_H_
#define INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_META_INFO_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag BandedChainInitialDPMatrix
// ----------------------------------------------------------------------------

// Specifies the first dp matrix that is computed.
struct BandedChainInitialDPMatrix_;
typedef Tag<BandedChainInitialDPMatrix_> BandedChainInitialDPMatrix;

// ----------------------------------------------------------------------------
// Tag BandedChainInnerDPMatrix
// ----------------------------------------------------------------------------

// Specifies the inner dp matrix that is computed.
struct BandedChainInnerDPMatrix_;
typedef Tag<BandedChainInnerDPMatrix_> BandedChainInnerDPMatrix;

// ----------------------------------------------------------------------------
// Tag BandedChainLastDPMatrix
// ----------------------------------------------------------------------------

// Specifies the lasr dp matrix that is computed.
struct BandedChainFinalDPMatrix_;
typedef Tag<BandedChainFinalDPMatrix_> BandedChainFinalDPMatrix;

// ----------------------------------------------------------------------------
// Struct BandedChainAlignment
// ----------------------------------------------------------------------------

// Algorithm tag for the banded chain alignments. The first type determines the
// free end-gaps and the second one which dp matrix is currently considered (initial, inner, final)
template <typename TSpec = FreeEndGaps_<>, typename TDPMatrixLocation = BandedChainInnerDPMatrix>
struct BandedChainAlignment_{};

// ----------------------------------------------------------------------------
// Class DPMetaColumn                                              [FullColumn]
// ----------------------------------------------------------------------------

template <typename TFreeEndGaps, typename TDPMatrixLocation, typename TGapCosts, typename TTraceback, typename TColumnType>
struct DPMetaColumn_<DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TTraceback>,
                    MetaColumnDescriptor<TColumnType, FullColumn> >
{
    // If InitialColumn -> replaced, replaced, replaced
    // If InnerColumn -> replaced, All, All
    // If FinalColumn -> replaced, All, All

    typedef RecursionDirectionZero TRecursionTypeFirstCell_;
    typedef RecursionDirectionAll TRecursionTypeInnerCell_;
    typedef RecursionDirectionAll TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

// ----------------------------------------------------------------------------
// Class DPMetaColumn                                        [PartialColumnTop]
// ----------------------------------------------------------------------------

template <typename TFreeEndGaps, typename TDPMatrixLocation, typename TGapCosts, typename TTraceback, typename TColumnType>
struct DPMetaColumn_<DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TTraceback>,
                    MetaColumnDescriptor<TColumnType, PartialColumnTop> >
{
    // If InitialColumn -> replaced, replaced, replaced
    // If InnerColumn -> replaced, All, LowerBand
    // If FinalColumn -> replaced, All, LowerBand

    typedef RecursionDirectionZero  TRecursionTypeFirstCell_;
    typedef RecursionDirectionAll TRecursionTypeInnerCell_;
    typedef RecursionDirectionLowerDiagonal TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

// ----------------------------------------------------------------------------
// Class DPMetaColumn                                      [PartialColumnMiddle]
// ----------------------------------------------------------------------------

template <typename TFreeEndGaps, typename TDPMatrixLocation, typename TGapCosts, typename TTraceback, typename TColumnType>
struct DPMetaColumn_<DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TTraceback>,
                    MetaColumnDescriptor<TColumnType, PartialColumnMiddle> >
{
    // If InitialColumn -> replaced, replaced, replaced
    // If InnerColumn -> UpperDiagonal, All, LowerDiagonal
    // If FinalColumn -> UpperDiagonal, All, LowerDiagonal

    typedef RecursionDirectionUpperDiagonal TRecursionTypeFirstCell_;
    typedef RecursionDirectionAll TRecursionTypeInnerCell_;
    typedef RecursionDirectionLowerDiagonal TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

// ----------------------------------------------------------------------------
// Class DPMetaColumn                                      [PartialColumnBottom]
// ----------------------------------------------------------------------------

template <typename TFreeEndGaps, typename TDPMatrixLocation, typename TGapCosts, typename TTraceback, typename TColumnType>
struct DPMetaColumn_<DPProfile_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TGapCosts, TTraceback>,
                    MetaColumnDescriptor<TColumnType, PartialColumnBottom> >
{
    // If InitialColumn -> replaced, replaced, replaced
    // If InnerColumn -> UpperDiagonal, All, All
    // If FinalColumn -> UpperDiagonal, All, All

    typedef RecursionDirectionUpperDiagonal TRecursionTypeFirstCell_;
    typedef RecursionDirectionAll TRecursionTypeInnerCell_;
    typedef RecursionDirectionAll TRecursionTypeLastCell_;

    typedef DPMetaCell_<TRecursionTypeFirstCell_, True> TFirstCell_;
    typedef DPMetaCell_<TRecursionTypeInnerCell_, True> TInnerCell_;
    typedef DPMetaCell_<TRecursionTypeLastCell_, True> TLastCell_;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction IsGlobalAlignment_
// ----------------------------------------------------------------------------

template <typename TSpec, typename TDPMatrixLocation>
struct IsGlobalAlignment_<BandedChainAlignment_<TSpec, TDPMatrixLocation> > : False{};

// ----------------------------------------------------------------------------
// Metafunction IsFreeEndGap_
// ----------------------------------------------------------------------------

template <typename TFreeEndGaps, typename TDPMatrixLocation, typename TDPSide>
struct IsFreeEndGap_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation>, TDPSide> :
    IsFreeEndGap_<TFreeEndGaps, TDPSide>
{};

template <typename TFreeEndGaps, typename TDPMatrixLocation, typename TDPSide>
struct IsFreeEndGap_<BandedChainAlignment_<TFreeEndGaps, TDPMatrixLocation> const, TDPSide> :
    IsFreeEndGap_<TFreeEndGaps const, TDPSide>
{};


// ----------------------------------------------------------------------------
// Metafunction SetupBandedChainAlignmentProfile_
// ----------------------------------------------------------------------------

// Profile for BandedChainAlignment algorithm
template <typename TAlignConfig, typename TGapCosts, typename TGapsPlacement>
struct SetupBandedChainAlignmentProfile_
{
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps_;
    typedef DPProfile_<BandedChainAlignment_<TFreeEndGaps_, BandedChainInnerDPMatrix>, TGapCosts, TracebackOn<TracebackConfig_<CompleteTrace,TGapsPlacement> > > Type;
};


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _setupAndRunBandedChainAlignment()
// ----------------------------------------------------------------------------

template <typename TTraceSegment, typename TStringSetSpec, typename TSeeds, typename TSequenceH, typename TSequenceV,
          typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap, bool TFirstRow, bool TFirstColumn,
          bool TLastColumn, bool TLastRow, typename TACSpec, typename TGapsPlacement>
inline TScoreValue
_setupAndRunBandedChainAlignment(StringSet<String<TTraceSegment>, TStringSetSpec> & globalTraceSet,
                                 TSeeds const & seedSet,
                                 TSequenceH const & seqH,
                                 TSequenceV const & seqV,
                                 Score<TScoreValue, TScoreSpecAnchor> const & scoringSchemeAnchor,
                                 Score<TScoreValue, TScoreSpecGap> const & scoringSchemeGap,
                                 AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const &,
                                 unsigned bandExtension,
                                 TGapsPlacement const &)
{
    //typedef typename Position<TSequenceH const>::Type TPosH;
    //typedef typename Position<TSequenceV const>::Type TPosV;
    typedef AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> TAlignConfig;
    typedef Score<TScoreValue, TScoreSpecAnchor> TScoringSchemeAnchor;
    typedef Score<TScoreValue, TScoreSpecGap> TScoringSchemeGap;
    typedef typename SequenceEntryForScore<TScoringSchemeAnchor, TSequenceH>::Type TSequenceHEntryAnchor;
    typedef typename SequenceEntryForScore<TScoringSchemeAnchor, TSequenceV>::Type TSequenceVEntryAnchor;
    typedef typename SequenceEntryForScore<TScoringSchemeGap, TSequenceH>::Type TSequenceHEntryGap;
    typedef typename SequenceEntryForScore<TScoringSchemeGap, TSequenceV>::Type TSequenceVEntryGap;

    TSequenceHEntryAnchor seqHEntryAnchor = sequenceEntryForScore(scoringSchemeAnchor, seqH, 0);
    TSequenceVEntryAnchor seqVEntryAnchor = sequenceEntryForScore(scoringSchemeAnchor, seqV, 0);
    TSequenceHEntryGap seqHEntryGap = sequenceEntryForScore(scoringSchemeGap, seqH, 0);
    TSequenceVEntryGap seqVEntryGap = sequenceEntryForScore(scoringSchemeGap, seqV, 0);

    if (scoreGapExtendHorizontal(scoringSchemeAnchor, seqHEntryAnchor, seqVEntryAnchor) !=
        scoreGapOpenHorizontal(scoringSchemeAnchor, seqHEntryAnchor, seqVEntryAnchor) ||
        scoreGapExtendVertical(scoringSchemeAnchor, seqHEntryAnchor, seqVEntryAnchor) !=
        scoreGapOpenVertical(scoringSchemeAnchor, seqHEntryAnchor, seqVEntryAnchor) ||
        scoreGapExtendHorizontal(scoringSchemeGap, seqHEntryGap, seqVEntryGap) !=
        scoreGapOpenHorizontal(scoringSchemeGap, seqHEntryGap, seqVEntryGap) ||
        scoreGapExtendVertical(scoringSchemeGap, seqHEntryGap, seqVEntryGap) !=
        scoreGapOpenVertical(scoringSchemeGap, seqHEntryGap, seqVEntryGap))
    {
        typedef typename SetupBandedChainAlignmentProfile_<TAlignConfig, AffineGaps, TGapsPlacement>::Type TDPProfile;
        return _computeAlignment(globalTraceSet, seedSet, seqH, seqV, scoringSchemeAnchor, scoringSchemeGap, bandExtension,
                                  TDPProfile());
    }
    else
    {
        typedef typename SetupBandedChainAlignmentProfile_<TAlignConfig, LinearGaps, TGapsPlacement>::Type TDPProfile;
        return _computeAlignment(globalTraceSet, seedSet, seqH, seqV, scoringSchemeAnchor, scoringSchemeGap, bandExtension,
                                  TDPProfile());
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_META_INFO_H_
