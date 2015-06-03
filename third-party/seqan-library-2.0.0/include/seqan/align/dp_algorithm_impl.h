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
// Implements the core of the dp algorithms.
// This is the crucial part of the refactoring of the alignment algorithms.
// It implements - at the moment only a column wise approach - the core
// loop structure for all alignment profiles. We generally differ between an
// unbanded alignment which is very easy, a banded alignment and a special
// case of the banded alignment the Hamming distance, where upper diagonal
// equals lower diagonal.
//
// The unbanded alignment:
// The computation of the unbanded alignment is divided into three parts.
// In the following we refer to a track as the part where the inner loop
// is iterating (in case of column wise navigation a track is equivalent
// to a column).
// First we compute the initial track. Afterwards we continue with all
// inner tracks of the dp matrix and in the end we compute the last track
// separately. This is because all three types have a special property that
// is different from the other track types.
// Each track itself is further divided into three parts, namely the first cell
// the inner cell and the last cell, corresponding to the initial row,
// all inner rows and the last row of a typical dp matrix. This partition of
// the dp matrix allows us to easily change the behavior of different cells
// according to the chosen dp profile at compile time.
// See alignment_dp_meta_info.h to learn about the different meta objects
// that manage the characteristics of each cell of a particular track type.
//
// The banded alignment:
// In the banded alignment we generally divide the dp matrix into the same
// partition as for the unbanded alignment. The only difference is that we,
// additionally add a specific information of how the current track is
// located within the dp matrix. Since we only consider a band we do not
// necessarily span over the full matrix size for a particular column.
// We distinguish between the locations: PartialColumnTop,
// PartialColumnMiddle, PartialColumnBottom and FullColumn (which is the
// default location for unbanded alignments). Each location of the column
// implies a different composition of the cells contained within a
// particular track. Thus, we are able to set different recursion
// directions and tracking informations for each cell independent from the
// runtime. The only difference is that the outer for-loop (iterating over
// the tracks) is split up into three loops. The first loop then only
// iterates over these tracks that are located at the top of the matrix.
// The second for-loop iterates over the tracks that either are of type
// PartialColumnMiddle or FullColumn (wide bands, where the upper diagonal
// begins behind the track where the lower diagonal crosses the last row of
// the dp matrix). And the last for-loop iterates over the tail of the band
// which is located at the PartialColumnBottom.
//
// The Hamming distance:
// In the special case where upper diagonal equals lower diagonal we only
// have to parse one diagonal of the matrix so we have a special
// implementation for that, though it works for all dp profiles.
//
// Restricitons:
// At the moment we have implemented a restriction such that not all bands
// are accepted. If the dp profile consists of the standard global alignment
// algorithm (NeedlemanWunsch or Gotoh), the band is required to go through
// the sink and the source of the dp matrix. If this is not given the
// alignment algorithm is aborted and the score MinValue<TScoreValue>::VALUE
// is returned.
// There are no further restrictions.
//
// GapCosts:
// Another detail of the new module is the selection of the gap functions,
// which is also now part of the compile time configuration. Whenever an
// algorithm is implemented it would automatically work for both gap
// functions (linear gaps and affine gaps).
//
// Navigation:
// It is possible to a certain degree to change the behavior of how to parse
// through the dp matrix. Using the new navigators one can implement
// different behaviors for different matrices. At the moment we only support
// column wise navigation for full and sparse score matrices and for full
// traceback matrices. Another detail of this navigators comes into account,
// when we want to compute only the score. We actually create a navigator
// for the dp matrix but implemented it this way that it gets never actually
// called when the traceback is disabled. Thus we do not store the traceback
// matrix if it is not necessary.
//
// Traceback:
// The traceback is now implemented as a single function that is used by all
// alignment profiles. Here we prefer the diagonal direction before the
// vertical before the horizontal direction.
// All tracebacks are first stored within the String<TraceSegment> object
// and afterwards, when the traceback is finished adapted to its given
// target object such as Align, Graph, Fragments, etc.
//
// Tracing:
// We use now an object called DPScout to keep track of the maximal score.
// This object scouts for the best value and can be overloaded to implement
// different strategies of how the score should be traced. Togehter with
// the meta_info file it only traces thus cells that are allowed to be
// traced given the current dp profile. Since this is also a compile time
// property we do not need to track every cell for the global alignment,
// while we do in the local alignment.
//
// Structure:
// The sequences within the matrix are marked as horizontal and vertical
// sequence to determine there orientation within the matrix.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_ALGORITHM_IMPL_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_DP_ALGORITHM_IMPL_H_

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
// Function _checkBandProperties()
// ----------------------------------------------------------------------------

// Checks whether the chosen band fits the dp profile.
template <typename TSequenceH, typename TSequenceV, typename TAlignmentProfile>
inline bool _checkBandProperties(TSequenceH const & /*seqH*/,
                                 TSequenceV const & /*seqV*/,
                                 DPBandConfig<BandOff> const & /*band*/,
                                 TAlignmentProfile const & /*alignProfile*/)
{
    return true;
}

template <typename TSequenceH, typename TSequenceV, typename TAlignmentProfile>
inline bool _checkBandProperties(TSequenceH const & seqH,
                                 TSequenceV const & seqV,
                                 DPBandConfig<BandOn> const & band,
                                 TAlignmentProfile const & /*alignProfile*/)
{
    typedef typename MakeSigned<typename Size<TSequenceH>::Type>::Type TSignedSize;

    // Check if the intersection between band and DP matrix is empty.
    if (upperDiagonal(band) < (0 - static_cast<TSignedSize>(length(seqV))) ||
        lowerDiagonal(band) > static_cast<TSignedSize>(length(seqH)))
    {
        return false;
    }

    // If the band begins before the beginning of the horizontal sequence
    // then check if free end-gaps are enabled at the beginning of the vertical sequence.
    if (upperDiagonal(band) < 0 && !IsFreeEndGap_<TAlignmentProfile, DPFirstColumn>::VALUE)
        return false;

    // If the band begins before the beginning of the vertical sequence
    // then check if free end-gaps are enabled at the beginning of the horizontal sequence.
    if (lowerDiagonal(band) > 0 && !IsFreeEndGap_<TAlignmentProfile, DPFirstRow>::VALUE)
        return false;

    // If the band ends behind the end of the vertical sequence
    // then check if free end-gaps are enabled at the end of the horizontal sequence.
    if (upperDiagonal(band) + static_cast<TSignedSize>(length(seqV)) < static_cast<TSignedSize>(length(seqH)) &&
        !IsFreeEndGap_<TAlignmentProfile, DPLastRow>::VALUE)
    {
        return false;
    }

    // If the band ends behind the end of the horizontal sequence
    // then check if free end-gaps are enabled at the end of the vertical sequence.
    if (lowerDiagonal(band) + static_cast<TSignedSize>(length(seqV)) > static_cast<TSignedSize>(length(seqH)) &&
        !IsFreeEndGap_<TAlignmentProfile, DPLastColumn>::VALUE)
    {
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------
// Function _invalidDPSettings()
// ----------------------------------------------------------------------------


// Checks if the settings for the dp algorithm are valid.
// Returns true if they are valid, false otherwise.
template <typename TSequenceH, typename TSequenceV, typename TBand, typename TAlignmentProfile>
inline bool _isValidDPSettings(TSequenceH const & seqH,
                               TSequenceV const & seqV,
                               TBand const & band,
                               TAlignmentProfile const & alignProfile)
{
    // Check if the sequences are empty.
    if (empty(seqH) || empty(seqV))
    {
        return false;
    }

    return _checkBandProperties(seqH, seqV, band, alignProfile);
}

// ----------------------------------------------------------------------------
// Function _isBandEnabled()
// ----------------------------------------------------------------------------

// Returns true if a band is selected, otherwise false.
template <typename TBandSpec>
inline bool
_isBandEnabled(DPBandConfig<TBandSpec> const & /*band*/)
{
    return IsSameType<TBandSpec, BandOn>::VALUE;
}

// ----------------------------------------------------------------------------
// Function _computeCell()
// ----------------------------------------------------------------------------

// Computes the score and tracks it if enabled.
template <typename TDPScout, typename TTraceMatrixNavigator, typename TScoreValue, typename TGapCosts,
          typename TSequenceHValue, typename TSequenceVValue, typename TScoringScheme, typename TColumnDescriptor,
          typename TCellDescriptor, typename TDPProfile>
inline void
_computeCell(TDPScout & scout,
             TTraceMatrixNavigator & traceMatrixNavigator,
             DPCell_<TScoreValue, TGapCosts> & activeCell,
             DPCell_<TScoreValue, TGapCosts> const & previousDiagonal,
             DPCell_<TScoreValue, TGapCosts> const & previousHorizontal,
             DPCell_<TScoreValue, TGapCosts> const & previousVertical,
             TSequenceHValue const & seqHVal,
             TSequenceVValue const & seqVVal,
             TScoringScheme const & scoringScheme,
             TColumnDescriptor const &,
             TCellDescriptor const &,   // One of FirstCell, InnerCell or LastCell.
             TDPProfile const &)
{
    typedef DPMetaColumn_<TDPProfile, TColumnDescriptor> TMetaColumn;
    assignValue(traceMatrixNavigator,
                _computeScore(activeCell, previousDiagonal, previousHorizontal, previousVertical, seqHVal, seqVVal,
                              scoringScheme, typename RecursionDirection_<TMetaColumn, TCellDescriptor>::Type(),
                              TDPProfile()));

    if (TrackingEnabled_<TMetaColumn, TCellDescriptor>::VALUE)
    {
        typedef typename IsSameType<typename TColumnDescriptor::TColumnProperty,
            DPFinalColumn>::Type TIsLastColumn;
        typedef typename And<IsSameType<TCellDescriptor, LastCell>, Or<
            IsSameType<typename TColumnDescriptor::TLocation,
            PartialColumnBottom>,IsSameType<
            typename TColumnDescriptor::TLocation, FullColumn> > >::Type
            TIsLastRow;
        _scoutBestScore(scout, activeCell, traceMatrixNavigator,
                        TIsLastColumn(), TIsLastRow());
    }
}

// ----------------------------------------------------------------------------
// Function _computeTrack()
// ----------------------------------------------------------------------------

// Computes one track of the dp algorithm. A track is defined as the area that is filled by the inner loop and
// iterated by the outer loop. For the column-wise navigation the track is equivalent with the column.
template <typename TDPScout, typename TDPScoreMatrixNavigator, typename TDPTraceMatrixNavigator, typename TSeqHValue,
          typename TSeqVValue, typename TSeqVIterator, typename TScoringScheme, typename TColumnDescriptor, typename TDPProfile>
inline void
_computeTrack(TDPScout & scout,
              TDPScoreMatrixNavigator & dpScoreMatrixNavigator,
              TDPTraceMatrixNavigator & dpTraceMatrixNavigator,
              TSeqHValue const & seqHValue,
              TSeqVValue const & seqVValue,
              TSeqVIterator const & seqBegin,
              TSeqVIterator const & seqEnd,
              TScoringScheme const & scoringScheme,
              TColumnDescriptor const &,
              TDPProfile const &)
{
    // Set the iterator to the begin of the track.
    _goNextCell(dpScoreMatrixNavigator, TColumnDescriptor(), FirstCell());
    _goNextCell(dpTraceMatrixNavigator, TColumnDescriptor(), FirstCell());

    // Compute the first cell.
    _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
                 previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
                 previousCellVertical(dpScoreMatrixNavigator), seqHValue, seqVValue, scoringScheme,
                 TColumnDescriptor(), FirstCell(), TDPProfile());
//    std::cerr << _scoreOfCell(value(dpScoreMatrixNavigator)) << " \t";
//    std::cerr << _scoreOfCell(value(dpScoreMatrixNavigator)) << "(" << _horizontalScoreOfCell(value(dpScoreMatrixNavigator)) << "," << _verticalScoreOfCell(value(dpScoreMatrixNavigator)) << ")" << " \t";

    TSeqVIterator iter = seqBegin;
    TSeqVIterator itEnd = (seqEnd - 1);
    // Compute the inner cells of the current track.
    for (; iter != itEnd; ++iter)
    {
        // Set the iterator to the next cell within the track.
        _goNextCell(dpScoreMatrixNavigator, TColumnDescriptor(), InnerCell());
        _goNextCell(dpTraceMatrixNavigator, TColumnDescriptor(), InnerCell());
        // Compute the inner cell.
        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
                     previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
                     previousCellVertical(dpScoreMatrixNavigator), seqHValue,
                     sequenceEntryForScore(scoringScheme, container(iter), position(iter)), scoringScheme,
                     TColumnDescriptor(), InnerCell(), TDPProfile());
//        std::cerr << _scoreOfCell(value(dpScoreMatrixNavigator)) << " \t";
//        std::cerr << _scoreOfCell(value(dpScoreMatrixNavigator)) << "(" << _horizontalScoreOfCell(value(dpScoreMatrixNavigator)) << "," << _verticalScoreOfCell(value(dpScoreMatrixNavigator)) << ")" << " \t";
    }
    // Set the iterator to the last cell of the track.
    _goNextCell(dpScoreMatrixNavigator, TColumnDescriptor(), LastCell());
    _goNextCell(dpTraceMatrixNavigator, TColumnDescriptor(), LastCell());
    // Compute the last cell.
    _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
                 previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
                 previousCellVertical(dpScoreMatrixNavigator), seqHValue,
                 sequenceEntryForScore(scoringScheme, container(iter), position(iter)), scoringScheme,
                 TColumnDescriptor(), LastCell(), TDPProfile());

//    std::cerr << _scoreOfCell(value(dpScoreMatrixNavigator)) << std::endl;
//    std::cerr << _scoreOfCell(value(dpScoreMatrixNavigator)) << "(" << _horizontalScoreOfCell(value(dpScoreMatrixNavigator)) << "," << _verticalScoreOfCell(value(dpScoreMatrixNavigator)) << ")" << std::endl;
}

// TODO(rmaerker): Debug code!
//template <typename TDPScout, typename TDPScoreMatrixNavigator, typename TDPTraceMatrixNavigator, typename TSeqHValue,
//            typename TSeqVValue, typename TSeqVIterator, typename TScoringScheme, typename TDPProfile,
//            typename TColumnDescriptor>
//inline void
//_computeTrack(TDPScout & scout,
//              TDPScoreMatrixNavigator & dpScoreMatrixNavigator,
//              TDPTraceMatrixNavigator & dpTraceMatrixNavigator,
//              TSeqHValue const & seqHValue,
//              TSeqVValue const & seqVValue,
//              TSeqVIterator const & seqBegin,
//              TSeqVIterator const & seqEnd,
//              TScoringScheme const & scoringScheme,
//              TColumnDescriptor const &,
//              TDPProfile const &,
//              unsigned col,
//              unsigned row,
//              String<std::string> & testMatrix)
//{
//    // Set the iterator to the begin of the track.
//    _goNextCell(dpScoreMatrixNavigator, TColumnDescriptor(), FirstCell());
//    _goNextCell(dpTraceMatrixNavigator, TColumnDescriptor(), FirstCell());
//
//    // Compute the first cell.
//    _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
//            previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
//            previousCellVertical(dpScoreMatrixNavigator), seqHValue, seqVValue, scoringScheme,
//            TColumnDescriptor(), FirstCell(), TDPProfile());
//    // TODO(rmaerker): remove debug code
////    std::cout << _scoreOfCell(value(dpScoreMatrixNavigator)) << "\t";
//    std::stringstream stream;
//    stream << _scoreOfCell(value(dpScoreMatrixNavigator));
////    stream << col + row;
//    testMatrix[col + row] = stream.str();
//    ++row;
//
//    TSeqVIterator iter = seqBegin;
//    TSeqVIterator itEnd = (seqEnd - 1);
//    // Compute the inner cells of the current track.
//    for (; iter != itEnd; ++iter, ++row)  // He will out of scope....
//    {
//      // Set the iterator to the next cell within the track.
//      _goNextCell(dpScoreMatrixNavigator, TColumnDescriptor(), InnerCell());
//      _goNextCell(dpTraceMatrixNavigator, TColumnDescriptor(), InnerCell());
//      // Compute the inner cell.
//        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
//                previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
//                previousCellVertical(dpScoreMatrixNavigator), seqHValue, value(iter), scoringScheme,
//                TColumnDescriptor(), InnerCell(), TDPProfile());
//        // TODO(rmaerker): remove debug code
////        std::cout << _scoreOfCell(value(dpScoreMatrixNavigator)) << "\t";
//        stream.str("");
//        stream << _scoreOfCell(value(dpScoreMatrixNavigator));
////        stream << col + row;
//        testMatrix[col + row] = stream.str();
//    }
//    // Set the iterator to the last cell of the track.
//    _goNextCell(dpScoreMatrixNavigator, TColumnDescriptor(), LastCell());
//    _goNextCell(dpTraceMatrixNavigator, TColumnDescriptor(), LastCell());
//    // Compute the last cell.
//    _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
//            previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
//            previousCellVertical(dpScoreMatrixNavigator), seqHValue, value(iter), scoringScheme,
//            TColumnDescriptor(), LastCell(), TDPProfile());
//    // TODO(rmaerker): remove debug code
////    std::cout << _scoreOfCell(value(dpScoreMatrixNavigator)) << "\n";
//    stream.str("");
//        stream << _scoreOfCell(value(dpScoreMatrixNavigator));
////        stream << col + row;
//    testMatrix[col + row] = stream.str();
//}

// ----------------------------------------------------------------------------
// Function _computeUnbandedAlignmentHelperTerminate()
// ----------------------------------------------------------------------------

template <typename TDPCell, typename TSpec>
inline bool //TODO(C++11) constexpr
_computeAlignmentHelperCheckTerminate(DPScout_<TDPCell, TSpec > const & /**/)
{
    return false;
}

template <typename TDPCell, typename TSpec>
inline bool
_computeAlignmentHelperCheckTerminate(DPScout_<TDPCell,Terminator_<TSpec> > const & s)
{
    return _terminationCriteriumIsMet(s);
}


// ----------------------------------------------------------------------------
// Function _computeUnbandedAlignment()
// ----------------------------------------------------------------------------

// Computes the standard DP-algorithm.
template <typename TDPScout, typename TDPScoreMatrixNavigator, typename TDPTraceMatrixNavigator, typename TSequenceH,
          typename TSequenceV, typename TScoringScheme, typename TAlignmentAlgo, typename TGapCosts, typename TTraceFlag>
inline void
_computeUnbandedAlignment(TDPScout & scout,
                          TDPScoreMatrixNavigator & dpScoreMatrixNavigator,
                          TDPTraceMatrixNavigator & dpTraceMatrixNavigator,
                          TSequenceH const & seqH,
                          TSequenceV const & seqV,
                          TScoringScheme const & scoringScheme,
                          DPProfile_<TAlignmentAlgo, TGapCosts, TTraceFlag> const & dpProfile)
{
    typedef typename Iterator<TSequenceH const, Rooted>::Type TConstSeqHIterator;
    typedef typename Iterator<TSequenceV const, Rooted>::Type TConstSeqVIterator;

    // ============================================================================
    // PREPROCESSING
    // ============================================================================

    TConstSeqVIterator seqVBegin = begin(seqV, Rooted());
    TConstSeqVIterator seqVEnd = end(seqV, Rooted());

    SEQAN_ASSERT_GT(length(seqH), 0u);
    SEQAN_ASSERT_GT(length(seqV), 0u);
    _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                  sequenceEntryForScore(scoringScheme, seqH, 0),
                  sequenceEntryForScore(scoringScheme, seqV, 0),
                  seqVBegin, seqVEnd, scoringScheme,
                  MetaColumnDescriptor<DPInitialColumn, FullColumn>(), dpProfile);

    // ============================================================================
    // MAIN DP
    // ============================================================================

    TConstSeqHIterator seqHIter = begin(seqH, Rooted());
    TConstSeqHIterator seqHIterEnd = end(seqH, Rooted()) - 1;

    for (; seqHIter != seqHIterEnd; ++seqHIter)
    {
        _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                      sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                      sequenceEntryForScore(scoringScheme, seqV, 0),
                      seqVBegin, seqVEnd, scoringScheme,
                      MetaColumnDescriptor<DPInnerColumn, FullColumn>(), dpProfile);

        if (_computeAlignmentHelperCheckTerminate(scout))
        {
                return;
        }
    }

    // ============================================================================
    // POSTPROCESSING
    // ============================================================================

    _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                  sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                  sequenceEntryForScore(scoringScheme, seqV, 0),
                  seqVBegin, seqVEnd, scoringScheme,
                  MetaColumnDescriptor<DPFinalColumn, FullColumn>(), dpProfile);

    // If we compute only the single option. we need to check if there are other possibilities at the end.
    // Traceback only from Diagonal, but could also come from vertical or horizontal.

//    for (unsigned i = 0; i < length(recMatrix); ++i)
//    {
//        std::cout << recMatrix[i]._score << "\t";
//    }
//    std::cout << std::endl;
}

// ----------------------------------------------------------------------------
// Function _computeBandedAlignment()
// ----------------------------------------------------------------------------

// Computes the banded DP-algorithm.
template <typename TDPScout, typename TDPScoreMatrixNavigator, typename TDPTraceMatrixNavigator, typename TSequenceH,
          typename TSequenceV, typename TScoringScheme, typename TBand, typename TAlignmentAlgo, typename TGapCosts,
          typename TTraceFlag>
inline void
_computeBandedAlignment(TDPScout & scout,
                        TDPScoreMatrixNavigator & dpScoreMatrixNavigator,
                        TDPTraceMatrixNavigator & dpTraceMatrixNavigator,
                        TSequenceH const & seqH,
                        TSequenceV const & seqV,
                        TScoringScheme const & scoringScheme,
                        TBand const & band,
                        DPProfile_<TAlignmentAlgo, TGapCosts, TTraceFlag> const & dpProfile)
{
    typedef DPProfile_<TAlignmentAlgo, TGapCosts, TTraceFlag> TDPProfile;
    typedef typename MakeSigned<typename Size<TSequenceH>::Type>::Type TSignedSizeSeqH;
    typedef typename MakeSigned<typename Size<TSequenceV>::Type>::Type TSignedSizeSeqV;
    typedef typename Iterator<TSequenceH const, Rooted>::Type TConstSeqHIterator;
    typedef typename Iterator<TSequenceV const, Rooted>::Type TConstSeqVIterator;


    // Now we have the problem of not knowing when we are in the last cell.

    // ============================================================================
    // PREPROCESSING
    // ============================================================================
    TSignedSizeSeqH seqHlength = static_cast<TSignedSizeSeqH>(length(seqH));
    TSignedSizeSeqH seqVlength = static_cast<TSignedSizeSeqV>(length(seqV));

    TConstSeqVIterator seqVBegin = begin(seqV, Rooted()) - _min(0, 1 + upperDiagonal(band));
    TConstSeqVIterator seqVEnd = begin(seqV, Rooted()) - _min(0, _max(-seqVlength, lowerDiagonal(band)));

    // We have to distinguish two band sizes. Some which spans the whole matrix in between and thus who not.
    // This can be distinguished, if UpperDiagonal > length(seqV) + LowerDiagonal

    // We start at least at the first position of the horizontal sequence or wherever the lower diagonal begins first.
    TConstSeqHIterator seqHIterBegin = begin(seqH, Rooted()) + _max(0, _min(seqHlength - 1, lowerDiagonal(band)));

    // The horizontal initial phase ends after the upper diagonal but at most after the horizontal sequence, or there is no horizontal initialization phase.
    TConstSeqHIterator seqHIterEndColumnTop = begin(seqH, Rooted()) + _min(seqHlength - 1, _max(0, upperDiagonal(band)));

    // The middle band phase ends after the lower diagonal crosses the bottom of the alignment matrix or after the horizontal sequence if it is smaller.
    TConstSeqHIterator seqHIterEndColumnMiddle = begin(seqH, Rooted()) + _min(seqHlength - 1, _max(0, seqVlength + lowerDiagonal(band)));
    // Swap the two iterators if we are in a band that spans over the full column.
    if (upperDiagonal(band) > seqVlength + lowerDiagonal(band))
        std::swap(seqHIterEndColumnTop, seqHIterEndColumnMiddle);

    // The bottom band phase ends after the upper diagonal of the band crosses the bottom of the matrix or after the horizontal sequence if it is smaller.
    TConstSeqHIterator seqHIterEndColumnBottom = begin(seqH, Rooted()) + _max(0, _min(seqHlength,
                                                                                      upperDiagonal(band) + seqVlength) - 1);

    // The Initial column can be PartialColumnTop which is given if the upper diagonal is >= 0,
    // otherwise it only can be PartialColumnMiddle or PartialColumnBottom depending where the lower diagonal is.

    // Check for single initialization cells in InitialColumn and FinalColumn.
    if (seqHIterBegin == end(seqH, Rooted()) - 1)
    {
        // Set the iterator to the begin of the track.
        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
        // Only one cell
        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
                     previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
                     previousCellVertical(dpScoreMatrixNavigator),
                     sequenceEntryForScore(scoringScheme, seqH, position(seqHIterBegin)),
                     sequenceEntryForScore(scoringScheme, seqV, 0), scoringScheme,
                     MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell(), TDPProfile());
        // we might need to additionally track this point.
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop> >, FirstCell>::VALUE)
            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator, True(), False());
        return;
    }
    if (seqHIterEndColumnBottom == begin(seqH, Rooted()))
    {
        // Set the iterator to the begin of the track.
        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell());
        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell());
        // Only one cell
        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
                     previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
                     previousCellVertical(dpScoreMatrixNavigator),
                     sequenceEntryForScore(scoringScheme, seqH, 0),
                     sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin)), scoringScheme,
                     MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell(), TDPProfile());
        // We might need to additionally track this point.
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom> >, LastCell>::VALUE)
            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator, False(), True());
        return;
    }

    if (upperDiagonal(band) < 0)
    {
        ++seqVBegin;
        if (lowerDiagonal(band) > -seqVlength)
            _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                          sequenceEntryForScore(scoringScheme, seqH, 0),
                          sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin) - 1),
                          seqVBegin, seqVEnd, scoringScheme,
                          MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), dpProfile);
        else
            _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                          sequenceEntryForScore(scoringScheme, seqH, 0),
                          sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin) - 1),
                          seqVBegin, seqVEnd, scoringScheme,
                          MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), dpProfile);
    }
    else if (lowerDiagonal(band) >= 0)
    {
        // Set the iterator to the begin of the track.
        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
        // Should we not just compute the cell?
        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
                     previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
                     previousCellVertical(dpScoreMatrixNavigator),
                     sequenceEntryForScore(scoringScheme, seqH, position(seqHIterBegin)),
                     sequenceEntryForScore(scoringScheme, seqV, 0),
                     scoringScheme,
                     MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell(), TDPProfile());
        // we might need to additionally track this point.
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop> >, FirstCell>::VALUE)
            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator, False(), False());
    }
    else  // Upper diagonal >= 0 and lower Diagonal < 0
    if (lowerDiagonal(band) <= -seqVlength)      // The band is bounded by the top and bottom of the matrix.
        _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                      sequenceEntryForScore(scoringScheme, seqH, 0),
                      sequenceEntryForScore(scoringScheme, seqV, 0),
                      seqVBegin, seqVEnd, scoringScheme,
                      MetaColumnDescriptor<DPInitialColumn, FullColumn>(), dpProfile);
    else       // The band is bounded by the top but not the bottom of the matrix.
        _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                      sequenceEntryForScore(scoringScheme, seqH, 0),
                      sequenceEntryForScore(scoringScheme, seqV, 0),
                      seqVBegin, seqVEnd, scoringScheme,
                      MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), dpProfile);

    if (_computeAlignmentHelperCheckTerminate(scout))
    {
            return;
    }

    // ============================================================================
    // MAIN DP
    // ============================================================================

    TConstSeqHIterator seqHIter = seqHIterBegin;
    // Compute the first part of the band, where the band is bounded by the top but not by the bottom of the matrix.
    for (; seqHIter != seqHIterEndColumnTop; ++seqHIter)
    {
        ++seqVEnd;
        _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                      sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                      sequenceEntryForScore(scoringScheme, seqV, 0),
                      seqVBegin, seqVEnd, scoringScheme,
                      MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), dpProfile);
        if (_computeAlignmentHelperCheckTerminate(scout))
        {
                return;
        }
    }

    // TODO(rmaerker): Check if putting the if-statement before the actual algorithm can speedup the code.
    // Check whether the band spans over the full column or not at some point.
    if (upperDiagonal(band) > seqVlength + lowerDiagonal(band))
    {
        // Compute the second part of the band, where the band is bounded by the top and the bottom of the matrix.
        // We might want to track the current cell here, since this is the first cell that crosses the bottom but is
        // not part of the FullColumn tracks.
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, FullColumn> >, LastCell>::VALUE)
            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator, False(), True());
        for (; seqHIter != seqHIterEndColumnMiddle; ++seqHIter)
        {
            _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                          sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                          sequenceEntryForScore(scoringScheme, seqV, 0),
                          seqVBegin, seqVEnd, scoringScheme,
                          MetaColumnDescriptor<DPInnerColumn, FullColumn>(), dpProfile);

            if (_computeAlignmentHelperCheckTerminate(scout))
            {
                    return;
            }
        }
    }
    else  // Compute the second part of the band, where the band is not bounded by the top and bottom of the matrix
    {
        for (; seqHIter != seqHIterEndColumnMiddle; ++seqHIter)
        {
            ++seqVBegin;
            ++seqVEnd;
            _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                          sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                          sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin) - 1),
                          seqVBegin, seqVEnd, scoringScheme,
                          MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), dpProfile);
            if (_computeAlignmentHelperCheckTerminate(scout))
            {
                    return;
            }
        }   // We might want to track the current cell here, since this is the first cell that crosses the bottom.
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom> >, LastCell>::VALUE)
        {
            if (lowerDiagonal(band) + seqVlength < seqHlength)
            {
                _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator, False(), True());
            }
        }

    }
    // Compute the third part of the band, where the band, is bounded by the bottom but not by the top of the matrix.
    for (; seqHIter != seqHIterEndColumnBottom; ++seqHIter)
    {
        ++seqVBegin;
        _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                      sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                      sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin) - 1),
                      seqVBegin, seqVEnd, scoringScheme,
                      MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), dpProfile);
        if (_computeAlignmentHelperCheckTerminate(scout))
        {
                return;
        }
    }

    // ============================================================================
    // POSTPROCESSING
    // ============================================================================

    // Check where the last track of the column is located.
    if (seqHIter - begin(seqH, Rooted()) < seqHlength - 1)  // Case 1: The band ends before the final column is reached.
    {
        // Set the iterator to the begin of the track.
        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell());
        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell());

        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
                     previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
                     previousCellVertical(dpScoreMatrixNavigator),
                     sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                     sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin)),
                     scoringScheme,
                     MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell(), TDPProfile());
        // We might need to additionally track this point.
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom> >, LastCell>::VALUE)
            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator, False(), True());
    }
    else if (seqHIter == end(seqH, Rooted()) - 1) // Case 2: The band ends somewhere in the final column of the matrix.
    {
        // Case2a: The band ends in the last cell of the final column.
        if (upperDiagonal(band) == seqHlength - seqVlength)
        {
            // Set the iterator to the begin of the track.
            _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell());
            _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell());

            _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
                         previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
                         previousCellVertical(dpScoreMatrixNavigator),
                         sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                         sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin)),
                         scoringScheme,
                         MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell(), TDPProfile());
            // we might need to additionally track this point.
            if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom> >, LastCell>::VALUE)
                _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator, True(), True());
        }
        else  // Case2b: At least two cells intersect between the band and the matrix in the final column of the matrix.
        {
            if (upperDiagonal(band) >= seqHlength)  // The band is bounded by the top of the matrix only or by the top and the bottom.
            {
                if (lowerDiagonal(band) + seqVlength > seqHlength) // The band is bounded by the top of the matrix
                {
                    ++seqVEnd;
                    _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                                  sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                                  sequenceEntryForScore(scoringScheme, seqV, 0),
                                  seqVBegin, seqVEnd, scoringScheme,
                                  MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), dpProfile);
                }
                else  // The band is bounded by the top and the bottom of the matrix.
                {
                    if (lowerDiagonal(band) + seqVlength + 1 > seqHlength)  // We have to go into the last cell.
                    {
                        ++seqVEnd;
                        _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                                      sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                                      sequenceEntryForScore(scoringScheme, seqV, 0),
                                      seqVBegin, seqVEnd, scoringScheme,
                                      MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), dpProfile);
                        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, FullColumn> >, LastCell>::VALUE)
                            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator, True(), True());
                    }
                    else
                        _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                                      sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                                      sequenceEntryForScore(scoringScheme, seqV, 0),
                                      seqVBegin, seqVEnd, scoringScheme,
                                      MetaColumnDescriptor<DPFinalColumn, FullColumn>(), dpProfile);
                }

            }
            else  // The band is bounded by bottom of matrix or completely unbounded.
            {
                ++seqVBegin;
                if (lowerDiagonal(band) + seqVlength <= seqHlength)  // The band is bounded by the bottom of the matrix.
                {
                    if (lowerDiagonal(band) + seqVlength == seqHlength)  // We have to go into the last cell.
                    {
                        ++seqVEnd;
                        _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                                      sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                                      sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin) - 1),
                                      seqVBegin, seqVEnd, scoringScheme,
                                      MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), dpProfile);
                        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom> >, LastCell>::VALUE)
                            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator, True(), True());
                    }
                    else
                    {
                        _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                                      sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                                      sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin) - 1),
                                      seqVBegin, seqVEnd, scoringScheme,
                                      MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), dpProfile);
                    }
                }
                else  // The band is unbounded by the matrix.
                {
                    ++seqVEnd;
                    _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
                                  sequenceEntryForScore(scoringScheme, seqH, position(seqHIter)),
                                  sequenceEntryForScore(scoringScheme, seqV, position(seqVBegin) - 1),
                                  seqVBegin, seqVEnd, scoringScheme,
                                  MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), dpProfile);
                }

            }
        }
    }
}

// TODO(rmaerker): This is denbug code only.
//template <typename TDPScout, typename TDPScoreMatrixNavigator, typename TDPTraceMatrixNavigator, typename TSequenceH,
//    typename TSequenceV, typename TScoringScheme, typename TBand, typename TAlignmentAlgo, typename TGapCosts,
//    typename TTraceFlag>
//inline void
//_debugBandedAlignment(TDPScout & scout,
//                        TDPScoreMatrixNavigator & dpScoreMatrixNavigator,
//                        TDPTraceMatrixNavigator & dpTraceMatrixNavigator,
//                        TSequenceH const & seqH,
//                        TSequenceV const & seqV,
//                        TScoringScheme const & scoringScheme,
//                        TBand const & band,
//                        DPProfile_<TAlignmentAlgo, TGapCosts, TTraceFlag> const & dpProfile)
//{
//    typedef DPProfile_<TAlignmentAlgo, TGapCosts, TTraceFlag> TDPProfile;
//    typedef typename Value<TSequenceH>::Type TSeqHValue;
//    typedef typename Value<TSequenceV>::Type TSeqVValue;
//    typedef typename Iterator<TSequenceH const, Standard>::Type TConstSeqHIterator;
//    typedef typename Iterator<TSequenceV const, Standard>::Type TConstSeqVIterator;
//
//
//    String<std::string> testMatrix;
//    resize(testMatrix, (length(seqH) + 1) * (length(seqV)+1));
//    // Now we have the problem of not knowing when we are in the last cell.
//
//    // INITIALIZATION
//    TConstSeqVIterator seqVBegin = begin(seqV, Standard()) - _min(0, 1+upperDiagonal(band));
//    TConstSeqVIterator seqVEnd = begin(seqV, Standard()) - _min(0, _max(-static_cast<int>(length(seqV)), lowerDiagonal(band)));
//
////    std::cout << "Begin Pos: " << seqVBegin - begin(seqV) << "\n";
////    std::cout << "End Pos: " << seqVEnd - begin(seqV) << "\n";
//
//    // We have to distinguish two band sizes. Some which spans the whole matrix in between and thus who not.
//    // This can be distinguished, if UpperDiagonal > length(seqV) + LowerDiagonal
//
//    // We start at the least at the first sequence or wherever the lower diagonal begins first.
//    TConstSeqHIterator seqHIterBegin = begin(seqH, Standard()) + _max(0, _min(static_cast<int>(length(seqH) - 1), lowerDiagonal(band)));
//    // TODO(rmaerker): Cehck if this assertion is correct.
////    SEQAN_ASSERT_NEQ(seqHIterBegin, end(seqH, Standard()));  // The iterator never points to the end of the horizontal sequence.
//
//    // The horizontal initial phase ends after the upper diagonal but at most after the horizontal sequence, or there is no horizontal initialization phase.
//    TConstSeqHIterator seqHIterEndColumnTop = begin(seqH, Standard()) + _min(static_cast<int>(length(seqH))-1, _max(0, upperDiagonal(band)));
//
//    // The middle band phase ends after the lower diagonal crosses the bottom of the alignment matrix or after the horizontal sequence if it is smaller.
//    TConstSeqHIterator seqHIterEndColumnMiddle = begin(seqH, Standard()) + _min(static_cast<int>(length(seqH))-1, _max(0,static_cast<int>(length(seqV)) + lowerDiagonal(band)));
//    // Swap the two iterators if we are in a band that spans over the full column.
//    if (upperDiagonal(band) > static_cast<int>(length(seqV)) + lowerDiagonal(band))
//        std::swap(seqHIterEndColumnTop, seqHIterEndColumnMiddle);
//
//    // The bottom band phase ends after the upper diagonal of the band crosses the bottom of the matrix or after the horizontal sequence if it is smaller.
//    TConstSeqHIterator seqHIterEndColumnBottom = begin(seqH, Standard()) + _max(0, _min(static_cast<int>(length(seqH)),
//                                                 upperDiagonal(band) + static_cast<int>(length(seqV))) -1);
//
////    std::cout << "seqHIterBegin Pos H: " << seqHIterBegin - begin(seqH) << "\n";
////    std::cout << "seqHIterEndColumnTop Pos H: " << seqHIterEndColumnTop - begin(seqH) << "\n";
////    std::cout << "seqHIterEndColumnMiddle Pos H: " << seqHIterEndColumnMiddle - begin(seqH) << "\n";
////    std::cout << "seqHIterEndColumnBottom Pos H: " << seqHIterEndColumnBottom - begin(seqH) << "\n";
////
////    std::cout << "_activeColIterator Pos: " << dpScoreMatrixNavigator._activeColIterator - begin(*dpScoreMatrixNavigator._ptrDataContainer) << "\n";
////    std::cout << "_prevColIterator Pos: " << dpScoreMatrixNavigator._prevColIterator - begin(*dpScoreMatrixNavigator._ptrDataContainer) << "\n";
//
//    // The Initial column can be PartialColumnTop which is given if the upper diagonal is >= 0,
//    // otherwise it only can be PartialColumnMiddle or PartialColumnBottom depending where the lower diagonal is.
//
//    // Have to check for single initialization cells in InitialColumn and FinalColumn.
//    if (seqHIterBegin == end(seqH)-1)
//    {
//        // Set the iterator to the begin of the track.
//        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
//        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
//        // Only one cell
//        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
//                previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
//                previousCellVertical(dpScoreMatrixNavigator), value(seqHIterBegin), TSeqVValue(), scoringScheme,
//                MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell(), TDPProfile());
////      std::cout << _scoreOfCell(value(dpScoreMatrixNavigator)) << "\n";
//      std::stringstream stream;
//      stream << _scoreOfCell(value(dpScoreMatrixNavigator));
//      testMatrix[length(seqH) * length(seqV)] = stream.str();
//        // we might need to additionally track this point.
//        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, PartialColumnTop> >, FirstCell>::VALUE)
//            _scoutBestScore(scout, _scoreOfCell(value(dpScoreMatrixNavigator)), dpTraceMatrixNavigator);
//        return;
//    }
//    if (seqHIterEndColumnBottom == begin(seqH))
//    {
//        // Set the iterator to the begin of the track.
//        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell());
//        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell());
//        // Only one cell
//        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
//                previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
//                previousCellVertical(dpScoreMatrixNavigator), TSeqHValue(), value(seqVBegin), scoringScheme,
//                MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), FirstCell(), TDPProfile());
////      std::cout << _scoreOfCell(value(dpScoreMatrixNavigator)) << "\n";
//      std::stringstream stream;
//      stream << _scoreOfCell(value(dpScoreMatrixNavigator));
//      testMatrix[length(seqH) * length(seqV)] = stream.str();
//        // We might need to additionally track this point.
//        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom> >, LastCell>::VALUE)
//            _scoutBestScore(scout, _scoreOfCell(value(dpScoreMatrixNavigator)), dpTraceMatrixNavigator);
//        return;
//    }
//
//    if (upperDiagonal(band) < 0)
//    {
//        ++seqVBegin;
//        if (lowerDiagonal(band) > -static_cast<int>(length(seqV)))
//          _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                  TSeqHValue(), value(seqVBegin-1), seqVBegin, seqVEnd, scoringScheme,
//                  MetaColumnDescriptor<DPInitialColumn, PartialColumnMiddle>(), dpProfile, 0, seqVBegin - begin(seqV) + 1, testMatrix);
//        else
//          _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                  TSeqHValue(), value(seqVBegin-1), seqVBegin, seqVEnd, scoringScheme,
//                  MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom>(), dpProfile, 0, seqVBegin - begin(seqV) + 1, testMatrix);
//    }
//    else if (lowerDiagonal(band) >= 0)
//    {
//        // Set the iterator to the begin of the track.
//        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
//        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell());
//        // Should we not just compute the cell?
//        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
//                previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
//                previousCellVertical(dpScoreMatrixNavigator), value(seqHIterBegin), TSeqVValue(), scoringScheme,
//                MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), FirstCell(), TDPProfile());
//        // we might need to additionally track this point.
////      std::cout << _scoreOfCell(value(dpScoreMatrixNavigator)) << "\n";
//      int col = lowerDiagonal(band);
//      std::stringstream stream;
//      stream << _scoreOfCell(value(dpScoreMatrixNavigator));
//      testMatrix[col * (length(seqV) + 1)] = stream.str();
//        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, PartialColumnTop> >, FirstCell>::VALUE)
//            _scoutBestScore(scout, _scoreOfCell(value(dpScoreMatrixNavigator)), dpTraceMatrixNavigator);
//    }
//    else  // Upper diagonal >= 0 and lower Diagonal < 0
//        if (lowerDiagonal(band) <= -static_cast<int>(length(seqV)))  // The band is bounded by the top and bottom of the matrix.
//          _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                  TSeqHValue(), TSeqVValue(), seqVBegin, seqVEnd, scoringScheme,
//                  MetaColumnDescriptor<DPInitialColumn, FullColumn>(), dpProfile, 0, 0, testMatrix);
//        else   // The band is bounded by the top but not the bottom of the matrix.
//          _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                  TSeqHValue(), TSeqVValue(), seqVBegin, seqVEnd, scoringScheme,
//                  MetaColumnDescriptor<DPInitialColumn, PartialColumnTop>(), dpProfile, 0,0, testMatrix);
//
//
//    // RECURSION
//    TConstSeqHIterator seqHIter = seqHIterBegin;
//    // Compute the first part of the band, where the band is bounded by the top but not by the bottom of the matrix.
//    for (;seqHIter != seqHIterEndColumnTop; ++seqHIter)
//    {
////      std::cout << value(seqHIter) << ":\t";
//        ++seqVEnd;
//      _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//              value(seqHIter), TSeqVValue(), seqVBegin, seqVEnd, scoringScheme,
//              MetaColumnDescriptor<DPInnerColumn, PartialColumnTop>(), dpProfile, ((seqHIter - begin(seqH)) + 1)*(length(seqV) +1), 0, testMatrix);
//    }
//
//    // Check whether the band spans over the full column or not at some point.
//    if (upperDiagonal(band) > static_cast<int>(length(seqV)) + lowerDiagonal(band))
//    {
//        // Compute the second part of the band, where the band is bounded by the top and the bottom of the matrix.
//        // We might want to track the current cell here, since this is the first cell that crosses the bottom.
//        // TODO(rmaerker): We have to check if the initial column/ row can be tracked!!!
//        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, FullColumn> >, LastCell>::VALUE)
//                _scoutBestScore(scout, _scoreOfCell(value(dpScoreMatrixNavigator)), dpTraceMatrixNavigator);
//        for (;seqHIter != seqHIterEndColumnMiddle; ++seqHIter)
//        {
//          _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                  value(seqHIter), TSeqVValue(), seqVBegin, seqVEnd, scoringScheme,
//                  MetaColumnDescriptor<DPInnerColumn, FullColumn>(), dpProfile, (seqHIter - begin(seqH) + 1) * (length(seqV) + 1), 0, testMatrix);
//        }
//    }
//    else // Compute the second part of the band, where the band is not bounded by the top and bottom of the matrix
//    {
//        for (;seqHIter != seqHIterEndColumnMiddle; ++seqHIter)
//        {
//            ++seqVBegin;
//            ++seqVEnd;
//          _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                  value(seqHIter), value(seqVBegin-1), seqVBegin, seqVEnd, scoringScheme,
//                  MetaColumnDescriptor<DPInnerColumn, PartialColumnMiddle>(), dpProfile, (seqHIter - begin(seqH) +1) * (length(seqV)+1), seqVBegin - begin(seqV), testMatrix);
//        }   // We might want to track the current cell here, since this is the first cell that crosses the bottom.
//        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom> >, LastCell>::VALUE)
//        {
//            // TODO(rmaerker): This is only a hot fix.
//            if (lowerDiagonal(band) + length(seqV) < length(seqH))
//                _scoutBestScore(scout, _scoreOfCell(value(dpScoreMatrixNavigator)), dpTraceMatrixNavigator);
//        }
//    }
//    // Compute the third part of the band, where the band, is bounded by the bottom but not by the top of the matrix.
//    for (;seqHIter != seqHIterEndColumnBottom; ++seqHIter)
//    {
//        ++seqVBegin;
//      _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//              value(seqHIter), value(seqVBegin-1), seqVBegin, seqVEnd, scoringScheme,
//              MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), dpProfile, (seqHIter - begin(seqH) + 1) * (length(seqV) + 1), seqVBegin - begin(seqV), testMatrix);
//    }
//
//    // Where ends the last cell?
//    if(seqHIter - begin(seqH) < static_cast<int>(length(seqH))-1)  // Case 1: The band ends before the final column is reached.
//    {
//        // Set the iterator to the begin of the track.
//        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell());
//        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell());
//
//        _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
//                previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
//                previousCellVertical(dpScoreMatrixNavigator), value(seqHIter), value(seqVBegin), scoringScheme,
//                MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom>(), FirstCell(), TDPProfile());
////      std::cout << _scoreOfCell(value(dpScoreMatrixNavigator)) << "\n";
//      std::stringstream stream;
//      stream << _scoreOfCell(value(dpScoreMatrixNavigator));
//      testMatrix[(seqHIter - begin(seqH) + 1) * (length(seqV) + 1) + length(seqV)] = stream.str();
//        // We might need to additionally track this point.
//        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, PartialColumnBottom> >, LastCell>::VALUE)
//            _scoutBestScore(scout, _scoreOfCell(value(dpScoreMatrixNavigator)), dpTraceMatrixNavigator);
//    }
//    else if(seqHIter == end(seqH)-1)  // Case 2: The band ends somewhere in the final column of the matrix.
//    {
//        if (upperDiagonal(band) == static_cast<int>(length(seqH))-static_cast<int>(length(seqV)))  // Case2a: The band ends in the last cell of the final column.  // He should be here....
//        {
//            // Set the iterator to the begin of the track.
//            _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell());
//            _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell());
//
//            _computeCell(scout, dpTraceMatrixNavigator, value(dpScoreMatrixNavigator),
//                    previousCellDiagonal(dpScoreMatrixNavigator), previousCellHorizontal(dpScoreMatrixNavigator),
//                    previousCellVertical(dpScoreMatrixNavigator), value(seqHIter), value(seqVBegin), scoringScheme,
//                    MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), FirstCell(), TDPProfile());
////          std::cout << _scoreOfCell(value(dpScoreMatrixNavigator)) << "\n";
//          std::stringstream stream;
//          stream << _scoreOfCell(value(dpScoreMatrixNavigator));
//          testMatrix[(length(seqH) * (length(seqV) + 1)) + length(seqV)] = stream.str();
//            // we might need to additionally track this point.
//            if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom> >, LastCell>::VALUE)
//                _scoutBestScore(scout, _scoreOfCell(value(dpScoreMatrixNavigator)), dpTraceMatrixNavigator);
//        }
//        else  // Case2b: At least two cells intersect between the band and the matrix in the final column of the matrix.
//        {
//            if (upperDiagonal(band) >= static_cast<int>(length(seqH)))  // The band is bounded by the top of the matrix only or by the top and the bottom.
//            {
//                if (lowerDiagonal(band) + static_cast<int>(length(seqV)) > static_cast<int>(length(seqH))) // The band is bounded by the top of the matrix
//                {
//                    ++seqVEnd;
//                  _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                          value(seqHIter), TSeqVValue(), seqVBegin, seqVEnd, scoringScheme,
//                          MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), dpProfile, length(seqH) * (length(seqV) + 1), 0, testMatrix);
//                }
//                else  // The band is bounded by the top and the bottom of the matrix.
//                {
//                    if (lowerDiagonal(band) + static_cast<int>(length(seqV)) + 1 > static_cast<int>(length(seqH)) )
//                    {
//                        ++seqVEnd;
//                      _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                              value(seqHIter), TSeqVValue(), seqVBegin, seqVEnd, scoringScheme,
//                              MetaColumnDescriptor<DPFinalColumn, PartialColumnTop>(), dpProfile, length(seqH) * (length(seqV) + 1), 0, testMatrix);
//                        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, FullColumn> >, LastCell>::VALUE)
//                            _scoutBestScore(scout, _scoreOfCell(value(dpScoreMatrixNavigator)), dpTraceMatrixNavigator);
//                    }
//                    else
//                      _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                              value(seqHIter), TSeqVValue(), seqVBegin, seqVEnd, scoringScheme,
//                              MetaColumnDescriptor<DPFinalColumn, FullColumn>(), dpProfile, length(seqH) * (length(seqV) +1), 0, testMatrix);
//
//                }
//
//            }
//            else  // The band is bounded by bottom of matrix or completely unbounded.
//            {
//                ++seqVBegin;
//                if (lowerDiagonal(band) + length(seqV) <= length(seqH)) // The band is bounded by the bottom of the matrix.
//                {
//                    if (lowerDiagonal(band) + length(seqV) == length(seqH))
//                    {
//                        ++seqVEnd;
//                      _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                          value(seqHIter), value(seqVBegin-1), seqVBegin, seqVEnd, scoringScheme,
//                          MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), dpProfile, length(seqH) * (length(seqV)+1), seqVBegin - begin(seqV), testMatrix);
//                      if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom> >, LastCell>::VALUE)
//                          _scoutBestScore(scout, _scoreOfCell(value(dpScoreMatrixNavigator)), dpTraceMatrixNavigator);
//                    }
//                    else
//                    {
//                      _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                          value(seqHIter), value(seqVBegin-1), seqVBegin, seqVEnd, scoringScheme,
//                          MetaColumnDescriptor<DPFinalColumn, PartialColumnBottom>(), dpProfile, length(seqH) * (length(seqV) + 1), seqVBegin - begin(seqV), testMatrix);
//                    }
//                }
//                else  // The band is unbounded by the matrix.
//                {
//                    ++seqVEnd;
//                  _computeTrack(scout, dpScoreMatrixNavigator, dpTraceMatrixNavigator,
//                      value(seqHIter), value(seqVBegin-1), seqVBegin, seqVEnd, scoringScheme,
//                      MetaColumnDescriptor<DPFinalColumn, PartialColumnMiddle>(), dpProfile, length(seqH) * (length(seqV)+1), seqVBegin - begin(seqV), testMatrix);
//                }
//
//            }
//        }
//    }
//
//    for (unsigned i = 0; i <= length(seqV); ++i)
//    {
//      for (unsigned j = 0; j <= length(seqH); ++j)
//      {
//          unsigned pos = j * (length(seqV) + 1) + i;
//          //if (!testMatrix[pos].empty())
//              std::cout << testMatrix[pos] << "\t" << std::flush;
//      }
//      std::cout << std::endl;
//    }
//
//}

// ----------------------------------------------------------------------------
// Function _computeHammingDistance()
// ----------------------------------------------------------------------------

// Computes the Hamming-Distance if the band-size is 1.
template <typename TDPScout, typename TDPScoreMatrixNavigator, typename TDPTraceMatrixNavigator, typename TSequenceH,
          typename TSequenceV, typename TScoringScheme, typename TBand, typename TAlignmentAlgo, typename TGapCosts,
          typename TTraceFlag>
inline void
_computeHammingDistance(TDPScout & scout,
                        TDPScoreMatrixNavigator & dpScoreMatrixNavigator,
                        TDPTraceMatrixNavigator & dpTraceMatrixNavigator,
                        TSequenceH const & seqH,
                        TSequenceV const & seqV,
                        TScoringScheme const & scoringScheme,
                        TBand const & band,
                        DPProfile_<TAlignmentAlgo, TGapCosts, TTraceFlag> const &)
{
    typedef typename MakeSigned<typename Size<TSequenceH const>::Type>::Type TSignedSizeSeqH;
    typedef typename MakeSigned<typename Size<TSequenceV const>::Type>::Type TSignedSizeSeqV;
    typedef typename Iterator<TSequenceH const, Rooted>::Type TConstSeqHIterator;
    typedef typename Iterator<TSequenceV const, Rooted>::Type TConstSeqVIterator;
    typedef typename Value<TDPScoreMatrixNavigator>::Type TDPCell;
    typedef DPProfile_<TAlignmentAlgo, TGapCosts, TTraceFlag> TDPProfile;

    // ============================================================================
    // PREPROCESSING
    // ============================================================================

    TSignedSizeSeqH seqHlength = static_cast<TSignedSizeSeqH>(length(seqH));
    TSignedSizeSeqH seqVlength = static_cast<TSignedSizeSeqV>(length(seqV));

    TConstSeqHIterator itH = begin(seqH, Rooted()) + _max(0, _min(seqHlength - 1, upperDiagonal(band)));
    TConstSeqHIterator itHEnd = begin(seqH, Rooted()) + _min(seqHlength - 1, upperDiagonal(band) + seqVlength);

    TConstSeqVIterator itV = begin(seqV, Rooted()) + _max(0, _min(seqVlength - 1, -lowerDiagonal(band)));
    TConstSeqVIterator itVEnd = begin(seqV, Rooted()) + _min(seqVlength - 1, lowerDiagonal(band) + seqHlength);

    assignValue(dpTraceMatrixNavigator,
                _computeScore(value(dpScoreMatrixNavigator), TDPCell(), TDPCell(), TDPCell(),
                              sequenceEntryForScore(scoringScheme, seqH, position(itH)),
                              sequenceEntryForScore(scoringScheme, seqV, position(itV)),
                              scoringScheme, RecursionDirectionZero(), TDPProfile()));

    if (upperDiagonal(band) < 0)
    {
        if (upperDiagonal(band) == -seqVlength)
        {
            if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInitialColumn, FullColumn> >, LastCell>::VALUE)
                _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator);
            return;
        }
        else
        {
            if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInitialColumn, FullColumn> >, InnerCell>::VALUE)
                _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator);
        }
    }
    else if (lowerDiagonal(band) > 0)
    {
        if (lowerDiagonal(band) == seqHlength)
        {
            if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, FullColumn> >, FirstCell>::VALUE)
                _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator);
            return;
        }
        else
        {
            if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, FullColumn> >, FirstCell>::VALUE)
                _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator);
        }
    }
    else
    {
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInitialColumn, FullColumn> >, FirstCell>::VALUE)
            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator);
    }

    TDPCell prevDiagonal = value(dpScoreMatrixNavigator);

    // ============================================================================
    // MAIN DP
    // ============================================================================

    while (itH != itHEnd && itV != itVEnd)
    {
        _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), FirstCell());
        _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), FirstCell());
        assignValue(dpTraceMatrixNavigator,
                    _computeScore(value(dpScoreMatrixNavigator), prevDiagonal, TDPCell(), TDPCell(),
                                  sequenceEntryForScore(scoringScheme, seqH, position(itH)),
                                  sequenceEntryForScore(scoringScheme, seqV, position(itV)),
                                  scoringScheme, RecursionDirectionDiagonal(), TDPProfile()));
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, FullColumn> >, InnerCell>::VALUE)
            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator);
        prevDiagonal = value(dpScoreMatrixNavigator);
        ++itH;
        ++itV;
    }

    // ============================================================================
    // POSTPROCESSING
    // ============================================================================

    _goNextCell(dpScoreMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), FirstCell());
    _goNextCell(dpTraceMatrixNavigator, MetaColumnDescriptor<DPInnerColumn, FullColumn>(), FirstCell());

    assignValue(dpTraceMatrixNavigator,
                _computeScore(value(dpScoreMatrixNavigator), prevDiagonal, TDPCell(), TDPCell(),
                              sequenceEntryForScore(scoringScheme, seqH, position(itH)),
                              sequenceEntryForScore(scoringScheme, seqV, position(itV)),
                              scoringScheme, RecursionDirectionDiagonal(), TDPProfile()));

    if (itH == itHEnd)
    {
        if (itV == itVEnd)   // Is in the last cell of final column
        {
            if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, FullColumn> >, LastCell>::VALUE)
                _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator);
        }
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPInnerColumn, FullColumn> >, LastCell>::VALUE)
            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator);
    }
    else
    {
        if (TrackingEnabled_<DPMetaColumn_<TDPProfile, MetaColumnDescriptor<DPFinalColumn, FullColumn> >, InnerCell>::VALUE)
            _scoutBestScore(scout, value(dpScoreMatrixNavigator), dpTraceMatrixNavigator);
    }
}

// ----------------------------------------------------------------------------
// Function _printTracebackMatrix()
// ----------------------------------------------------------------------------

template <typename TTraceMatrix>
void _printTracebackMatrix(TTraceMatrix & dpTraceMatrix)
{
    typedef typename Size<TTraceMatrix>::Type TSize;
    TSize dimH = length(dpTraceMatrix, +DPMatrixDimension_::HORIZONTAL);
    TSize dimV = length(dpTraceMatrix, +DPMatrixDimension_::VERTICAL);

    for (unsigned row = 0; row < dimV; ++row)
    {
        for (unsigned column = 0; column < dimH; ++column)
        {
            std::cout << _translateTraceValue(value(dpTraceMatrix, row + column * dimV)) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template <typename TTraceNavigator, typename TScoreValue, typename TDPScoutSpec>
inline void
_correctTraceValue(TTraceNavigator &,
                   DPScout_<DPCell_<TScoreValue, LinearGaps>, TDPScoutSpec> const &)
{
    // Nothing to do.
}

template <typename TTraceNavigator, typename TScoreValue, typename TDPScoutSpec>
inline void
_correctTraceValue(TTraceNavigator & traceNavigator,
                   DPScout_<DPCell_<TScoreValue, AffineGaps>, TDPScoutSpec>  const & dpScout)
{
    _setToPosition(traceNavigator, maxHostPosition(dpScout));
    if (_verticalScoreOfCell(dpScout._maxScore) == _scoreOfCell(dpScout._maxScore))
    {
        value(traceNavigator) &= ~TraceBitMap_::DIAGONAL;
        value(traceNavigator) |= TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
    }
    else if (_horizontalScoreOfCell(dpScout._maxScore) == _scoreOfCell(dpScout._maxScore))
    {
        value(traceNavigator) &= ~TraceBitMap_::DIAGONAL;
        value(traceNavigator) |= TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
    }
}

template <typename TTraceNavigator, typename TScoreValue, typename TDPScoutSpec>
inline void
_correctTraceValue(TTraceNavigator & traceNavigator,
                   DPScout_<DPCell_<TScoreValue, DynamicGaps>, TDPScoutSpec>  const & dpScout)
{
    _setToPosition(traceNavigator, maxHostPosition(dpScout));
    if (isGapExtension(dpScout._maxScore, DynamicGapExtensionVertical()))
    {
        value(traceNavigator) &= ~TraceBitMap_::DIAGONAL;
        value(traceNavigator) |= TraceBitMap_::MAX_FROM_VERTICAL_MATRIX;
    }
    else if (isGapExtension(dpScout._maxScore, DynamicGapExtensionHorizontal()))
    {
        value(traceNavigator) &= ~TraceBitMap_::DIAGONAL;
        value(traceNavigator) |= TraceBitMap_::MAX_FROM_HORIZONTAL_MATRIX;
    }
}

// ----------------------------------------------------------------------------
// Function _computeAligmnment()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TGapScheme, typename TTraceTarget, typename TScoutState, typename TSequenceH, typename TSequenceV,
          typename TScoreScheme, typename TBandSwitch, typename TAlignmentAlgorithm, typename TTraceFlag>
inline typename Value<TScoreScheme>::Type
_computeAlignment(DPContext<TScoreValue, TGapScheme> & dpContext,
                  TTraceTarget & traceSegments,
                  TScoutState & scoutState,
                  TSequenceH const & seqH,
                  TSequenceV const & seqV,
                  TScoreScheme const & scoreScheme,
                  DPBandConfig<TBandSwitch> const & band,
                  DPProfile_<TAlignmentAlgorithm, TGapScheme, TTraceFlag> const & dpProfile)
{
    typedef typename GetDPScoreMatrix<DPContext<TScoreValue, TGapScheme> >::Type TDPScoreMatrixHost;
    typedef typename Value<TDPScoreMatrixHost>::Type TDPScoreValue;

    typedef typename GetDPTraceMatrix<DPContext<TScoreValue, TGapScheme> >::Type TDPTraceMatrixHost;
    typedef typename Value<TDPTraceMatrixHost>::Type TTraceValue;

    typedef typename DefaultScoreMatrixSpec_<TAlignmentAlgorithm>::Type TScoreMatrixSpec;

    typedef DPMatrix_<TDPScoreValue, TScoreMatrixSpec> TDPScoreMatrix;
    typedef DPMatrix_<TTraceValue, FullDPMatrix> TDPTraceMatrix;

    typedef DPMatrixNavigator_<TDPScoreMatrix, DPScoreMatrix, NavigateColumnWise> TDPScoreMatrixNavigator;
    typedef DPMatrixNavigator_<TDPTraceMatrix, DPTraceMatrix<TTraceFlag>, NavigateColumnWise> TDPTraceMatrixNavigator;

    typedef typename ScoutSpecForAlignmentAlgorithm_<TAlignmentAlgorithm>::Type TDPScoutSpec;
    typedef DPScout_<TDPScoreValue, TDPScoutSpec> TDPScout;

    // Check if current dp settings are valid. If not return infinity value for dp score value.
    if (!_isValidDPSettings(seqH, seqV, band, dpProfile))
        return MinValue<TScoreValue>::VALUE;

    TDPScoreMatrix dpScoreMatrix;
    TDPTraceMatrix dpTraceMatrix;

    // TODO(rmaerker): Check whether the matrix allocation can be reduced if upperDiagonal < 0?
    setLength(dpScoreMatrix, +DPMatrixDimension_::HORIZONTAL, length(seqH) + 1 - std::max(0, lowerDiagonal(band)));
    setLength(dpTraceMatrix, +DPMatrixDimension_::HORIZONTAL, length(seqH) + 1 - std::max(0, lowerDiagonal(band)));

    if (IsSameType<TBandSwitch, BandOff>::VALUE)
    {
        setLength(dpScoreMatrix, +DPMatrixDimension_::VERTICAL, length(seqV) + 1);
        setLength(dpTraceMatrix, +DPMatrixDimension_::VERTICAL, length(seqV) + 1);
    }
    else
    {
        int bandSize = _min(static_cast<int>(length(seqH)), upperDiagonal(band)) - _max(lowerDiagonal(band), -static_cast<int>(length(seqV))) + 1;
        setLength(dpScoreMatrix, +DPMatrixDimension_::VERTICAL, _min(static_cast<int>(length(seqV)) + 1, bandSize));
        setLength(dpTraceMatrix, +DPMatrixDimension_::VERTICAL, _min(static_cast<int>(length(seqV)) + 1, bandSize));
    }

    // We set the host to the score matrix and the dp matrix.
    setHost(dpScoreMatrix, getDpScoreMatrix(dpContext));
    setHost(dpTraceMatrix, getDpTraceMatrix(dpContext));

    resize(dpScoreMatrix);
    // We do not need to allocate the memory for the trace matrix if the traceback is disabled.
    if (IsTracebackEnabled_<TTraceFlag>::VALUE)
        resize(dpTraceMatrix);

    TDPScoreMatrixNavigator dpScoreMatrixNavigator;
    TDPTraceMatrixNavigator dpTraceMatrixNavigator;

    _init(dpScoreMatrixNavigator, dpScoreMatrix, band);
    _init(dpTraceMatrixNavigator, dpTraceMatrix, band);

    TDPScout dpScout(scoutState);

    // Execute the alignment.
    if (!_isBandEnabled(band))
        _computeUnbandedAlignment(dpScout, dpScoreMatrixNavigator, dpTraceMatrixNavigator, seqH, seqV, scoreScheme, dpProfile);
    else if (upperDiagonal(band) == lowerDiagonal(band))
        _computeHammingDistance(dpScout, dpScoreMatrixNavigator, dpTraceMatrixNavigator, seqH, seqV, scoreScheme, band, dpProfile);
    else
        _computeBandedAlignment(dpScout, dpScoreMatrixNavigator, dpTraceMatrixNavigator, seqH, seqV, scoreScheme,
                                band, dpProfile);

    if (IsSameType<TTraceFlag, TracebackOff>::VALUE)
        return maxScore(dpScout);

    if (IsSingleTrace_<TTraceFlag>::VALUE)
    {
        // Check if max was found at the bottom right corner of the matrix.
        // This is also true if in last row, and last column
//        if ((maxHostPosition(dpScout) + 1) == (end(dpTraceMatrix) - begin(dpTraceMatrix)))

//            maxHostPosition(dpScout); // We only have the trace value not the score value.
            _correctTraceValue(dpTraceMatrixNavigator, dpScout);
    }

//    _printTracebackMatrix(dpTraceMatrix);
    _computeTraceback(traceSegments, dpTraceMatrixNavigator, dpScout, seqH, seqV, band, dpProfile);

    return maxScore(dpScout);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_DP_ALGORITHM_IMPL_H_
