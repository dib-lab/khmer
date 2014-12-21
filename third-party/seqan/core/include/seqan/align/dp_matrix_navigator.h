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
// Implements the DPMatrixNavigator.
// This is the parent class for all navigators that parse over a dp-matrix.
// This class facilitates the correct navigation through the dp-matrix for
// all kind of standard dp-algorithms. At the moment there only exists the
// NavigateColumnWise but this can be complemented by other navigation
// structures like anti-diagonals or in tiles.
// The Navigator can be specialized with three parameters. The first one is
// the used specialization of the dp-matrix (FullDPMatrix or SparseDPMatrix).
// the second parameter decides if it is a navigator for a score matrix or
// a trace matrix. And the last parameter determines the sort of navigation.
// ==========================================================================

// TODO(holtgrew): Documentation in this header necessary or internal only?

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag DPScoreMatrix
// ----------------------------------------------------------------------------

// Used to select a navigator for the score matrix.
struct DPScoreMatrix_;
typedef Tag<DPScoreMatrix_> DPScoreMatrix;

// ----------------------------------------------------------------------------
// Tag DPTraceMatrix
// ----------------------------------------------------------------------------

// Used to select a navigator for the trace matrix.
template <typename TTraceFlag>
struct DPTraceMatrix {};

// ----------------------------------------------------------------------------
// Tag NavigateColumnWise
// ----------------------------------------------------------------------------

// Facilitates column wise navigation through the dp-matrix.
struct NavigateColumnWise_;
typedef Tag<NavigateColumnWise_> NavigateColumnWise;

// ----------------------------------------------------------------------------
// Class DPMatrixNavigator_
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec = NavigateColumnWise>
class DPMatrixNavigator_;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
struct Value<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> >
{
    typedef typename Value<TDPMatrix>::Type Type;
};

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
struct Value<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const>
{
    typedef typename Value<TDPMatrix const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
struct Reference<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> >
{
    typedef typename Reference<TDPMatrix>::Type Type;
};

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
struct Reference<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const>
{
    typedef typename Reference<TDPMatrix const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
struct Container<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> >
{
    typedef TDPMatrix Type;
};

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
struct Container<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const>
{
    typedef TDPMatrix const Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec, typename TValue>
inline void
assignValue(DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> & dpNavigator,
            TValue const & element)
{
    assignValue(dpNavigator._activeColIterator, element);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> >::Type
value(DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> & dpNavigator)
{
    return value(dpNavigator._activeColIterator);
}

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const>::Type
value(DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const & dpNavigator)
{
    return value(dpNavigator._activeColIterator);
}

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec, typename TCoordinateV, typename TCoordinateH>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> >::Type
value(DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> & dpNavigator,
      TCoordinateH const & coordinateV,
      TCoordinateV const & coordinateH)
{
    return value(container(dpNavigator), coordinateV, coordinateH);
}

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec, typename TCoordinateV, typename TCoordinateH>
inline typename Reference<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const>::Type
value(DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const & dpNavigator,
      TCoordinateV const & coordinateV,
      TCoordinateH const & coordinateH)
{
    return value(container(dpNavigator), coordinateV, coordinateH);
}

// ----------------------------------------------------------------------------
// Function coordinate()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
inline typename DPMatrixDimension_::TValue
coordinate(DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const & dpNavigator,
           typename DPMatrixDimension_::TValue const & dimension)
{
    // Simply delegate to coordinate of underlying matrix.
    return coordinate(value(dpNavigator._ptrDataContainer),
                      dpNavigator._activeColIterator - begin(*dpNavigator._ptrDataContainer, Standard()), dimension);
}

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
inline typename Container<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> >::Type &
container(DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> & dpNavigator)
{
    return *dpNavigator._ptrDataContainer;
}

template <typename TDPMatrix, typename TDPMatrixType, typename TNavigationSpec>
inline typename Container<DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const>::Type &
container(DPMatrixNavigator_<TDPMatrix, TDPMatrixType, TNavigationSpec> const & dpNavigator)
{
    return *dpNavigator._ptrDataContainer;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_NAVIGATOR_H_
