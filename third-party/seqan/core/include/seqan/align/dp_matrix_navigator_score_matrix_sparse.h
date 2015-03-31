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
// The navigator for the sparse score dp-matrix. This class also provides an
// iterator for the active and the previous column. It stores the neighbouring
// cells needed for the recursion formula.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_SPARSE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_SPARSE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPMatrixNavigator                        [SparseDPMatrix, ScoreMatrix]
// ----------------------------------------------------------------------------

// Specialization of the score matrix navigator for a sparse dp matrix.
template <typename TValue>
class DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise>
{
public:
    typedef  DPMatrix_<TValue, SparseDPMatrix> TDPMatrix_;
    typedef typename Pointer_<TDPMatrix_>::Type TDPMatrixPointer_;
    typedef typename Iterator<TDPMatrix_, Standard>::Type TDPMatrixIterator;

    TDPMatrixPointer_ _ptrDataContainer;        // Pointer to the underlying matrix to navigate on.
    int _laneLeap;                              // The distance to leap when going to the next column.
    TDPMatrixIterator _activeColIterator;       // The iterator over the active column.
    TDPMatrixIterator _prevColIterator;         // The iterator over the previous column.
    TValue _prevCellDiagonal;                   // The previous value in diagonal direction.
    TValue _prevCellHorizontal;                 // The previous value in horizontal direction.
    TValue _prevCellVertical;                   // The previous value in vertical direction.



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


// Initializes the navigator for unbanded alignments
template <typename TValue>
inline void
_init(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & navigator,
      DPMatrix_<TValue, SparseDPMatrix> & dpMatrix,
      DPBand_<BandOff> const &)
{
    navigator._ptrDataContainer = &dpMatrix;
    navigator._activeColIterator = begin(dpMatrix, Standard());
    navigator._prevColIterator = navigator._activeColIterator;
    navigator._laneLeap = 1 - _dataLengths(dpMatrix)[DPMatrixDimension_::VERTICAL];
}

// Initializes the navigator for banded alignments
template <typename TValue>
inline void
_init(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & navigator,
      DPMatrix_<TValue, SparseDPMatrix> & dpMatrix,
      DPBand_<BandOn> const & band)
{
    typedef DPMatrix_<TValue, SparseDPMatrix> TSparseDPMatrix;
    typedef typename Size<TSparseDPMatrix>::Type TSize;
    typedef typename MakeSigned<TSize>::Type TSignedSize;
    navigator._ptrDataContainer = &dpMatrix;

    // Band begins within the first row.
    if (lowerDiagonal(band) >= 0)
    {
        navigator._laneLeap = 0;
        navigator._activeColIterator = begin(dpMatrix, Standard()) + length(dpMatrix, DPMatrixDimension_::VERTICAL) - 1;
    }
    else if (upperDiagonal(band) <= 0) // Band begins within the first column
    {
        navigator._laneLeap = 1 - _dataLengths(dpMatrix)[DPMatrixDimension_::VERTICAL];
        navigator._activeColIterator = begin(dpMatrix, Standard());
    }
    else  // Band intersects with the point of origin.
    {
        navigator._laneLeap = _max(lowerDiagonal(band), 1 - static_cast<TSignedSize>(length(dpMatrix, DPMatrixDimension_::VERTICAL)));
        navigator._activeColIterator = begin(dpMatrix, Standard()) + length(dpMatrix, DPMatrixDimension_::VERTICAL) + navigator._laneLeap - 1;
    }
    navigator._prevColIterator = navigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell()        [DPInitialColumn, PartialColumnTop, FirstCell]
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, PartialColumnTop> const &,
            FirstCell const &)
{
    // no-op
}

// ----------------------------------------------------------------------------
// Function _goNextCell()              [DPInitialColumn, FullColumn, FirstCell]
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            FirstCell const &)
{
    // no-op
}

// ----------------------------------------------------------------------------
// Function _goNextCell()                          [DPInitialColumn, FirstCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & /*dpNavigator*/,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            FirstCell const &)
{
    // no-op
}

// ----------------------------------------------------------------------------
// Function _goNextCell()                         [PartialColumnTop, FirstCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, PartialColumnTop> const &,
            FirstCell const &)
{
    --dpNavigator._laneLeap;
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevColIterator = dpNavigator._activeColIterator;
    dpNavigator._prevCellHorizontal = value(++dpNavigator._prevColIterator);
}

// ----------------------------------------------------------------------------
// Function _goNextCell()                               [FullColumn, FirstCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            FirstCell const &)
{
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevCellHorizontal = value(dpNavigator._activeColIterator);
}

// ----------------------------------------------------------------------------
// Function _goNextCell()                                           [FirstCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            FirstCell const &)
{
    dpNavigator._activeColIterator += dpNavigator._laneLeap;
    dpNavigator._prevColIterator = dpNavigator._activeColIterator;
    dpNavigator._prevCellDiagonal = value(dpNavigator._prevColIterator);
    dpNavigator._prevCellHorizontal = value(++dpNavigator._prevColIterator);
}

// ----------------------------------------------------------------------------
// Function _goNextCell                            [DPInitialColumn, InnerCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            InnerCell const &)
{
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                [DPInitialColumn, FullColumn, InnerCell]
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            InnerCell const &)
{
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                             [InnerCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            InnerCell const &)
{
    dpNavigator._prevCellDiagonal = dpNavigator._prevCellHorizontal;
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    dpNavigator._prevCellHorizontal = value(++dpNavigator._prevColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                 [FullColumn, InnerCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            InnerCell const &)
{
    dpNavigator._prevCellDiagonal = dpNavigator._prevCellHorizontal;
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    dpNavigator._prevCellHorizontal = value(++dpNavigator._activeColIterator);
}

// ----------------------------------------------------------------------------
// Function _goNextCell                             [DPInitialColumn, LastCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, TColumnLocation> const &,
            LastCell const &)
{
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell        [DPInitialColumn, PartialColumnBottom, LastCell]
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, PartialColumnBottom> const &,
            LastCell const &)
{
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                 [DPInitialColumn, FullColumn, LastCell]
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<DPInitialColumn, FullColumn> const &,
            LastCell const &)
{
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                                              [LastCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnType, typename TColumnLocation>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, TColumnLocation> const &,
            LastCell const &)
{
    dpNavigator._prevCellDiagonal = dpNavigator._prevCellHorizontal;
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    ++dpNavigator._activeColIterator;
}

// ----------------------------------------------------------------------------
// Function _goNextCell                         [PartialColumnBottom, LastCell]
// ----------------------------------------------------------------------------

template <typename TValue, typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
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
// Function _goNextCell                                  [FullColumn, LastCell]
// ----------------------------------------------------------------------------


template <typename TValue, typename TColumnType>
inline void
_goNextCell(DPMatrixNavigator_<DPMatrix_<TValue, SparseDPMatrix>, DPScoreMatrix, NavigateColumnWise> & dpNavigator,
            MetaColumnDescriptor<TColumnType, FullColumn> const &,
            LastCell const &)
{
    dpNavigator._prevCellDiagonal = dpNavigator._prevCellHorizontal;
    dpNavigator._prevCellVertical = value(dpNavigator._activeColIterator);
    dpNavigator._prevCellHorizontal = value(++dpNavigator._activeColIterator);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_SCORE_SPARSE_H_
