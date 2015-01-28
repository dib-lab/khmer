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
// The navigator for the full score dp-matrix. We need two iterators over the
// current column and the previous column. We also store the three neighboring
// cells needed for the recursion formula.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPMatrixNavigator                          [FullDPMatrix, ScoreMatrix]
// ----------------------------------------------------------------------------

// The navigator for the score matrix.
//
// This navigator runs on a FullDPMatrix while it navigates column wise.
template <typename TValue>
class DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise>
{
public:
    typedef  DPMatrix_<TValue, FullDPMatrix> TDPMatrix_;
    typedef typename Pointer_<TDPMatrix_>::Type TDPMatrixPointer_;
    typedef typename Iterator<TDPMatrix_, Standard>::Type TDPMatrixIterator;

    TDPMatrixPointer_ _ptrDataContainer;        // Pointer to the matrix this navigator is working on.
    int _laneLeap;                              // Stores the jump to the next column
    TDPMatrixIterator _activeColIterator;       // The active column iterator.
    TDPMatrixIterator _prevColIterator;         // The previous column iterator.
    TValue _prevCellDiagonal;                   // The previous diagonal cell
    TValue _prevCellHorizontal;                 // The previous Horizontal cell
    TValue _prevCellVertical;                   // The previous Vertical cell



    DPMatrixNavigator_() :
        _ptrDataContainer(TDPMatrixPointer_(0)),
        _laneLeap(0),
        _activeColIterator(),
        _prevColIterator(),
        _prevCellDiagonal(),
        _prevCellHorizontal(),
        _prevCellVertical()
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _init()
// ----------------------------------------------------------------------------

// Initializes the navigator for an unbanded alignment.
template <typename TValue>
inline void
_init(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & navigator,
      DPMatrix_<TValue, FullDPMatrix> & dpMatrix,
      DPBand_<BandOff> const &)
{
    navigator._ptrDataContainer = &dpMatrix;
    navigator._activeColIterator = begin(dpMatrix, Standard());
    navigator._prevColIterator = navigator._activeColIterator - _dataFactors(dpMatrix)[DPMatrixDimension_::HORIZONTAL];
    navigator._laneLeap = 1;
}

// Initializes the navigator for a banded alignment.
template <typename TValue>
inline void
_init(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & navigator,
      DPMatrix_<TValue, FullDPMatrix> & dpMatrix,
      DPBand_<BandOn> const & band)
{
    typedef typename Size<DPMatrix_<TValue, FullDPMatrix> >::Type TMatrixSize;
    typedef typename MakeSigned<TMatrixSize>::Type TSignedSize;
    navigator._ptrDataContainer = &dpMatrix;


    // Band begins within the first row.
    if (lowerDiagonal(band) >= 0)
    {
        navigator._laneLeap = _min(length(dpMatrix, DPMatrixDimension_::VERTICAL), bandSize(band));
        navigator._activeColIterator = begin(dpMatrix, Standard()) + _dataLengths(dpMatrix)[DPMatrixDimension_::VERTICAL] - 1;
    }
    else if (upperDiagonal(band) <= 0)  // Band begins within the first column.
    {
        navigator._laneLeap = 1;
        navigator._activeColIterator = begin(dpMatrix, Standard());
    }
    else  // Band intersects with the point of origin.
    {
        TMatrixSize lengthVertical = length(dpMatrix, DPMatrixDimension_::VERTICAL);
        int lastPos = _max(-static_cast<TSignedSize>(lengthVertical - 1), lowerDiagonal(band));
        navigator._laneLeap = lengthVertical + lastPos;
        navigator._activeColIterator = begin(dpMatrix, Standard()) + navigator._laneLeap - 1;
    }
    // Set previous iterator to same position, one column left.
    navigator._prevColIterator = navigator._activeColIterator - _dataFactors(dpMatrix)[DPMatrixDimension_::HORIZONTAL];
}

// ----------------------------------------------------------------------------
// Function _goNextCell()                          [DPInitialColumn, FirstCell]
// ----------------------------------------------------------------------------

// In the initial column we don't need to do anything because, the navigagtor is already initialized.
template <typename TValue>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, PartialColumnTop> const &,
            FirstCell const &)
{
    // no-op
}

template <typename TValue>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            FirstCell const &)
{
    // no-op
}

template <typename TValue, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            FirstCell const &)
{
    // no-op
}

// ----------------------------------------------------------------------------
// Function _goNextCell()                         [PartialColumnTop, FirstCell]
// ----------------------------------------------------------------------------

// We are in the banded case, where the band crosses the first row.
// The left cell of the active cell is not valid, beacause we only can come from horizontal direction.
// The lower left cell of the active cell is the horizontal direction.
template <typename TValue, typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, PartialColumnTop> const &,
            FirstCell const &)
{
    --dpNavigator._laneLeap;
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevColIterator += dpNavigator._laneLeap;
    dpNavigator._prevCellHorizontal = value(++dpNavigator._prevColIterator);
}

// ----------------------------------------------------------------------------
// Function _goNextCell()                               [FullColumn, FirstCell]
// ----------------------------------------------------------------------------

// We are in the unbanded case or in the middle phase of the wide band.
// The left cell of the active cell represents horizontal direction.
template <typename TValue, typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            FirstCell const &)
{
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevColIterator += dpNavigator._laneLeap;
    dpNavigator._prevCellHorizontal = value(dpNavigator._prevColIterator);
}

// ----------------------------------------------------------------------------
// Function _goNextCell() [PartialColumnMiddle, PartialColumnBottom, FirstCell]
// ----------------------------------------------------------------------------

// We are in the banded case.
// The left cell of the active cell represents diagonal direction. The lower left diagonal represents the horizontal direction.

template <typename TValue, typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            FirstCell const &)
{
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevColIterator += dpNavigator._laneLeap;
    dpNavigator._prevCellDiagonal = value(dpNavigator._prevColIterator);
    dpNavigator._prevCellHorizontal = value(++dpNavigator._prevColIterator);
}

// ----------------------------------------------------------------------------
// Function _goNextCell                            [DPInitialColumn, InnerCell]
// ----------------------------------------------------------------------------

// If we are in the initial column, we only need to represent the vertical direction.
// But we still have to update the previous column iterator.
template <typename TValue, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            InnerCell const &)
{
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
    ++dpNavigator._prevColIterator;  // Do we have to increase the prevColIterator....
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                  [AnyColumn, InnerCell]
// ----------------------------------------------------------------------------

// For any other column type and location we can use the same navigation procedure.
template <typename TValue, typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            InnerCell const &)
{
    dpNavigator._prevCellDiagonal = dpNavigator._prevCellHorizontal;
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    dpNavigator._prevCellHorizontal = value(++dpNavigator._prevColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                             [DPInitialColumn, LastCell]
// ----------------------------------------------------------------------------

// If we are in the initial column we only need to represent the vertical direction.
// But we still have to update the previous column iterator.
template <typename TValue, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            LastCell const &)
{
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
    ++dpNavigator._prevColIterator;
}

// We need this function to avoid ambiguous function calls.
template <typename TValue>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom> const &,
            LastCell const &)
{
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
    ++dpNavigator._prevColIterator;
}

// We need this function to avoid ambiguous function calls.
template <typename TValue>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            LastCell const &)
{
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
    ++dpNavigator._prevColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                  [FullColumn, LastCell]
// ----------------------------------------------------------------------------

// If we are in a full column the values correspond to standard dp directions.
template <typename TValue, typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            LastCell const &)
{
    dpNavigator._prevCellDiagonal = dpNavigator._prevCellHorizontal;
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    dpNavigator._prevCellHorizontal = value(++dpNavigator._prevColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                         [PartialColumnBottom, LastCell]
// ----------------------------------------------------------------------------

// If we are in banded case and are the band crosses the last row, we have to update
// the additional leap for the current track.
template <typename TValue, typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, PartialColumnBottom> const &,
            LastCell const &)
{
    dpNavigator._prevCellDiagonal = dpNavigator._prevCellHorizontal;
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    dpNavigator._prevCellHorizontal = value(++dpNavigator._prevColIterator);
    ++dpNavigator._activeColIterator;
    ++dpNavigator._laneLeap;
}

// ----------------------------------------------------------------------------
// Function _goNextCell      [PartialColumnTop & PartialColumnBottom, LastCell]
// ----------------------------------------------------------------------------

// If we are in the banded case the left cell of the active represents the diagonal direction.
template <typename TValue, typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, FullDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            LastCell const &)
{
    dpNavigator._prevCellDiagonal = dpNavigator._prevCellHorizontal;
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function previousCellDiagonal()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> >::Type
previousCellDiagonal(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> & dpNavigator)
{
    return dpNavigator._prevCellDiagonal;
}

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const>::Type
previousCellDiagonal(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const & dpNavigator)
{
    return dpNavigator._prevCellDiagonal;
}

// ----------------------------------------------------------------------------
// Function previousCellHorizontal()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> >::Type
previousCellHorizontal(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> & dpNavigator)
{
    return dpNavigator._prevCellHorizontal;
}

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const>::Type
previousCellHorizontal(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const & dpNavigator)
{
    return dpNavigator._prevCellHorizontal;
}

// ----------------------------------------------------------------------------
// Function previousCellVertical()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> >::Type
previousCellVertical(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> & dpNavigator)
{
    return dpNavigator._prevCellVertical;
}

template <typename TDPMatrix, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const>::Type
previousCellVertical(DPMatrixNavigator_<TDPMatrix, DPScoreMatrix, TNavigationSpec> const & dpNavigator)
{
    return dpNavigator._prevCellVertical;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_H_
