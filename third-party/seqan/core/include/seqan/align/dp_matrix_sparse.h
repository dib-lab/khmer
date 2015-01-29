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

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_SPARSE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_SPARSE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class DPMatrix                                              [SparseDPMatrix]
// ----------------------------------------------------------------------------

template <typename TValue>
class DPMatrix_<TValue, SparseDPMatrix>
{
public:

    typedef Matrix<TValue, 2> THost;

    Holder<THost>   _dataHost;  // The host containing the actual matrix.

    DPMatrix_() :
        _dataHost()
    {
        create(_dataHost);
    }

    DPMatrix_(DPMatrix_ const & other) :
        _dataHost(other._dataHost) {}

    ~DPMatrix_() {}

    DPMatrix_ & operator=(DPMatrix_ const & other)
    {
        if (this != &other)
        {
            _dataHost = other._dataHost;
        }
        return *this;
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
resize(DPMatrix_<TValue, SparseDPMatrix> & dpMatrix)
{
    typedef DPMatrix_<TValue, SparseDPMatrix> TDPMatrix;
    typedef typename Size<TDPMatrix>::Type TSize;

    TSize _dimVertical = length(dpMatrix, DPMatrixDimension_::VERTICAL);

    if (_dimVertical > 0)
        resize(host(dpMatrix), _dimVertical, Exact());
}

template <typename TValue>
inline void
resize(DPMatrix_<TValue, SparseDPMatrix> & dpMatrix,
       TValue const & fillValue)
{
    typedef DPMatrix_<TValue, SparseDPMatrix> TDPMatrix;
    typedef typename Size<TDPMatrix>::Type TSize;

    TSize _dimVertical = length(dpMatrix, DPMatrixDimension_::VERTICAL);

    if (_dimVertical > 0)
        resize(host(dpMatrix), _dimVertical, fillValue, Exact());
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TValue, typename TPositionV, typename TPositionH>
inline typename Reference<DPMatrix_<TValue, SparseDPMatrix> >::Type
value(DPMatrix_<TValue, SparseDPMatrix> & dpMatrix,
      TPositionV const & posV,
      TPositionH const &)
{
    return value(dpMatrix, posV);
}

template <typename TValue, typename TPositionV, typename TPositionH>
inline typename Reference<DPMatrix_<TValue, SparseDPMatrix> const>::Type
value(DPMatrix_<TValue, SparseDPMatrix> const & dpMatrix,
      TPositionV const & posV,
      TPositionH const &)
{
    return value(dpMatrix, posV);
}

// ----------------------------------------------------------------------------
// Function coordinate()
// ----------------------------------------------------------------------------


template <typename TValue, typename TPosition>
inline typename Position<DPMatrix_<TValue, SparseDPMatrix> >::Type
coordinate(DPMatrix_<TValue, SparseDPMatrix> const & /*dpMatrix*/,
           TPosition hostPos,
           typename DPMatrixDimension_::TValue dimension)
{
    SEQAN_ASSERT(_checkCorrectDimension(dimension));

    if (dimension == DPMatrixDimension_::VERTICAL)
        return hostPos;

    return 0u;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_SPARSE_H_
