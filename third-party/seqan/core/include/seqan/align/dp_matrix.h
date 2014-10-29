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
// Author: Renï¿½ Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// This file implements the class DPMatrix and its specialization
// FullDPMatrix. The DPMatrix is a wrapper class around the Matrix<TValue,2>
// class. Thus we can implement different specializations for the dp-matrix
// that are used through the different dp-algorithms. For example, we need
// a full dp matrix to store the traceback or the score for the Waterman-
// Eggert algorithm, while for the other dp-algorithms we only need one
// column vector to compute the scores. The default dp-matrix specialization
// can be selected using a special meta-function.
// ==========================================================================

// TODO(holtgrew): Documentation in this header necessary or internal only?

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TAlgorithm>
struct DefaultScoreMatrixSpec_;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag SparseDPMatrix
// ----------------------------------------------------------------------------

struct SparseDPMatrix_;
typedef Tag<SparseDPMatrix_> SparseDPMatrix;

// ----------------------------------------------------------------------------
// Tag FullDPMatrix
// ----------------------------------------------------------------------------

struct FullDPMatrix_;
typedef Tag<FullDPMatrix_> FullDPMatrix;


// ----------------------------------------------------------------------------
// Enum DPMatrixDimension
// ----------------------------------------------------------------------------

// Used to globally specify the correct dimension and the correct size of
// dimension for the dp matrix.
struct DPMatrixDimension_
{
    typedef unsigned int TValue;

    static const TValue VERTICAL = 0u;
    static const TValue HORIZONTAL = 1u;
    static const TValue DIMENSION = 2u;
};

// ----------------------------------------------------------------------------
// Class DPMatrix_
// ----------------------------------------------------------------------------

// The dp matrix used as a score matrix and as a trace-back matrix.
template <typename TValue, typename TMatrixSpec>
class DPMatrix_
{};


// Default dp matrix implementation stores all cells of the dp matrix in the
// underlying two-dimensional matrix.
template <typename TValue>
class DPMatrix_<TValue, FullDPMatrix>
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

// ----------------------------------------------------------------------------
// Metafunction DefaultScoreMatrixSpec_
// ----------------------------------------------------------------------------

// This meta-function determines the default specialization of dp matrix
// based on the given algorithm tag.
template <typename TAlgorithm>
struct DefaultScoreMatrixSpec_
{
    typedef SparseDPMatrix Type;
};

// TODO(rmaerker): Move to WatermanEggert implementation?
template <>
struct DefaultScoreMatrixSpec_<LocalAlignment_<WatermanEggert> >
{
    typedef FullDPMatrix Type;
};

// ----------------------------------------------------------------------------
// Metafunction _DataHost
// ----------------------------------------------------------------------------

// Returns the type of the underlying matrix.
template <typename TDPMatrix>
struct _DataHost {};

template <typename TValue, typename TMatrixSpec>
struct _DataHost<DPMatrix_<TValue, TMatrixSpec> >
{
    typedef DPMatrix_<TValue, TMatrixSpec> TDPMatrix_;
    typedef typename TDPMatrix_::THost Type;
};

template <typename TValue, typename TMatrixSpec>
struct _DataHost<DPMatrix_<TValue, TMatrixSpec> const>
{
    typedef DPMatrix_<TValue, TMatrixSpec> TDPMatrix_;
    typedef typename TDPMatrix_::THost const Type;
};

// ----------------------------------------------------------------------------
// Metafunction _SizeArr
// ----------------------------------------------------------------------------

// Returns the type of the containers to store the dimensions and the factors
// in order to move properly in the matrix.
template <typename TDPMatrix>
struct _SizeArr {};

template <typename TValue, typename TMatrixSpec>
struct _SizeArr<DPMatrix_<TValue, TMatrixSpec> >
{
    typedef DPMatrix_<TValue, TMatrixSpec> TDPMatrix_;
    typedef typename _DataHost<TDPMatrix_>::Type TDataHost_;
    typedef typename SizeArr_<TDataHost_>::Type Type;
};

template <typename TValue, typename TMatrixSpec>
struct _SizeArr<DPMatrix_<TValue, TMatrixSpec> const>
{
    typedef DPMatrix_<TValue, TMatrixSpec> TDPMatrix_;
    typedef typename _DataHost<TDPMatrix_>::Type TDataHost_;
    typedef typename SizeArr_<TDataHost_>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
struct Spec<DPMatrix_<TValue, TMatrixSpec> >
{
    typedef TMatrixSpec Type;
};

template <typename TValue, typename TMatrixSpec>
struct Spec<DPMatrix_<TValue, TMatrixSpec> const>:
    Spec<DPMatrix_<TValue, TMatrixSpec> >{};


// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
struct Value<DPMatrix_<TValue, TMatrixSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TMatrixSpec>
struct Value<DPMatrix_<TValue, TMatrixSpec> const>
{
    typedef TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
struct Reference<DPMatrix_<TValue, TMatrixSpec> >
{
    typedef TValue & Type;
};

template <typename TValue, typename TMatrixSpec>
struct Reference<DPMatrix_<TValue, TMatrixSpec> const>
{
    typedef TValue const & Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
struct GetValue<DPMatrix_<TValue, TMatrixSpec> >:
    Reference<DPMatrix_<TValue, TMatrixSpec> >{};

template <typename TValue, typename TMatrixSpec>
struct GetValue<DPMatrix_<TValue, TMatrixSpec> const>:
    Reference<DPMatrix_<TValue, TMatrixSpec> const>{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
struct Position<DPMatrix_<TValue, TMatrixSpec> >
{
    typedef typename DPMatrix_<TValue, TMatrixSpec>::THost TDataMatrix_;
    typedef typename Position<TDataMatrix_>::Type Type;
};

template <typename TValue, typename TMatrixSpec>
struct Position<DPMatrix_<TValue, TMatrixSpec> const>:
    Position<DPMatrix_<TValue, TMatrixSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
struct Size<DPMatrix_<TValue, TMatrixSpec> >
{
    typedef typename DPMatrix_<TValue, TMatrixSpec>::THost TDataMatrix_;
    typedef typename Size<TDataMatrix_>::Type Type;
};

template <typename TValue, typename TMatrixSpec>
struct Size<DPMatrix_<TValue, TMatrixSpec> const>:
    Size<DPMatrix_<TValue, TMatrixSpec> >{};

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
struct Host<DPMatrix_<TValue, TMatrixSpec> >
{
    typedef DPMatrix_<TValue, TMatrixSpec> TDPMatrix_;
    typedef typename _DataHost<TDPMatrix_>::Type TDataMatrix_;
    typedef typename Host<TDataMatrix_>::Type Type;
};

template <typename TValue, typename TMatrixSpec>
struct Host<DPMatrix_<TValue, TMatrixSpec> const>
{
    typedef DPMatrix_<TValue, TMatrixSpec> TDPMatrix_;
    typedef typename _DataHost<TDPMatrix_>::Type TDataMatrix_;
    typedef typename Host<TDataMatrix_>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

// There are two iterator types. The standard iterator returns a
// non-rooted iterator to the underlying vector of the hosted two-dimensional
// matrix. The rooted iterator returns the iterator defined by the
// hosted matrix object which is a position iterator.
template <typename TValue, typename TMatrixSpec>
struct Iterator<DPMatrix_<TValue, TMatrixSpec>, Standard const>
{
    typedef DPMatrix_<TValue, TMatrixSpec> TDPMatrix_;
    typedef typename  Host<TDPMatrix_>::Type THost_;
    typedef typename Iterator<THost_, Standard>::Type Type;
};

template <typename TValue, typename TMatrixSpec>
struct Iterator<DPMatrix_<TValue, TMatrixSpec> const, Standard const>
{
    typedef DPMatrix_<TValue, TMatrixSpec> const TDPMatrix_;
    typedef typename  Host<TDPMatrix_>::Type THost_;
    typedef typename Iterator<THost_ const, Standard>::Type Type;
};

template <typename TValue, typename TMatrixSpec>
struct Iterator<DPMatrix_<TValue, TMatrixSpec>, Rooted const>
{
    typedef DPMatrix_<TValue, TMatrixSpec> TDPMatrix_;
    typedef typename  _DataHost<TDPMatrix_>::Type TDataMatrix_;
    typedef typename Iterator<TDataMatrix_, Rooted>::Type Type;
};

template <typename TValue, typename TMatrixSpec>
struct Iterator<DPMatrix_<TValue, TMatrixSpec> const, Rooted const>
{
    typedef DPMatrix_<TValue, TMatrixSpec> TDPMatrix_;
    typedef typename  _DataHost<TDPMatrix_>::Type TDataMatrix_;
    typedef typename Iterator<TDataMatrix_ const, Rooted>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _checkCorrectDimension()
// ----------------------------------------------------------------------------

// Checks whether a given value applies to the correct dimension.
inline bool _checkCorrectDimension(DPMatrixDimension_::TValue dim)
{
    return dim < DPMatrixDimension_::DIMENSION;
}

// ----------------------------------------------------------------------------
// Function _dataHost()
// ----------------------------------------------------------------------------

// Returns a reference to the hosted matrix.
template <typename TValue, typename TMatrixSpec>
inline typename _DataHost<DPMatrix_<TValue, TMatrixSpec> >::Type &
_dataHost(DPMatrix_<TValue, TMatrixSpec>&dpMatrix)
{
    return value(dpMatrix._dataHost);
}

template <typename TValue, typename TMatrixSpec>
inline typename _DataHost<DPMatrix_<TValue, TMatrixSpec> const>::Type &
_dataHost(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix)
{
    return value(dpMatrix._dataHost);
}

// ----------------------------------------------------------------------------
// Function _dataLengths()
// ----------------------------------------------------------------------------

// Returns a reference to the _dataLengths container of the hosted matrix.
template <typename TValue, typename TMatrixSpec>
inline typename _SizeArr<DPMatrix_<TValue, TMatrixSpec> >::Type &
_dataLengths(DPMatrix_<TValue, TMatrixSpec>&dpMatrix)
{
    return _dataLengths(_dataHost(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline typename _SizeArr<DPMatrix_<TValue, TMatrixSpec> const>::Type &
_dataLengths(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix)
{
    return _dataLengths(_dataHost(dpMatrix));
}

// ----------------------------------------------------------------------------
// Function _dataFactors()
// ----------------------------------------------------------------------------

// Returns a reference to the _dataFactors container of the hosted matrix.
template <typename TValue, typename TMatrixSpec>
inline typename _SizeArr<DPMatrix_<TValue, TMatrixSpec> >::Type &
_dataFactors(DPMatrix_<TValue, TMatrixSpec>&dpMatrix)
{
    return _dataFactors(_dataHost(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline typename _SizeArr<DPMatrix_<TValue, TMatrixSpec> const>::Type &
_dataFactors(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix)
{
    return _dataFactors(_dataHost(dpMatrix));
}

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

// Returns a reference to the underlying vector of the hosted matrix.
template <typename TValue, typename TMatrixSpec>
inline typename Host<DPMatrix_<TValue, TMatrixSpec> >::Type &
host(DPMatrix_<TValue, TMatrixSpec>&dpMatrix)
{
    return host(_dataHost(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline typename Host<DPMatrix_<TValue, TMatrixSpec> const>::Type &
host(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix)
{
    return host(_dataHost(dpMatrix));
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

// Sets a new value to the underlying vector of the hosted matrix.
template <typename TValue, typename TMatrixSpec, typename THost>
inline void
setHost(DPMatrix_<TValue, TMatrixSpec> & dpMatrix,
        THost & newHost)
{
    setHost(_dataHost(dpMatrix), newHost);
}

template <typename TValue, typename TMatrixSpec, typename THost>
inline void
setHost(DPMatrix_<TValue, TMatrixSpec> & dpMatrix,
        THost const & newHost)
{
    setHost(_dataHost(dpMatrix), newHost);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

// Returns the value of the matrix at the given host position.
template <typename TValue, typename TMatrixSpec, typename TPosition>
inline typename Reference<DPMatrix_<TValue, TMatrixSpec> >::Type
value(DPMatrix_<TValue, TMatrixSpec> & dpMatrix,
      TPosition const & pos)
{
    return value(_dataHost(dpMatrix), pos);
}

template <typename TValue, typename TMatrixSpec, typename TPosition>
inline typename Reference<DPMatrix_<TValue, TMatrixSpec> const>::Type
value(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix,
      TPosition const & pos)
{
    return value(_dataHost(dpMatrix), pos);
}

// Returns the value of the matrix at the two given coordinates.
template <typename TValue, typename TMatrixSpec, typename TPositionV, typename TPositionH>
inline typename Reference<DPMatrix_<TValue, TMatrixSpec> >::Type
value(DPMatrix_<TValue, TMatrixSpec> & dpMatrix,
      TPositionV const & posDimV,
      TPositionH const & posDimH)
{
    return value(_dataHost(dpMatrix), posDimV, posDimH);
}

template <typename TValue, typename TMatrixSpec, typename TPositionV, typename TPositionH>
inline typename Reference<DPMatrix_<TValue, TMatrixSpec> const>::Type
value(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix,
      TPositionV const & posDimV,
      TPositionH const & posDimH)
{
    return value(_dataHost(dpMatrix), posDimV, posDimH);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

// Returns the length of a given dimension.
template <typename TValue, typename TMatrixSpec>
inline typename Size<DPMatrix_<TValue, TMatrixSpec> const>::Type
length(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix,
       unsigned int dimension)
{
    SEQAN_ASSERT(_checkCorrectDimension(dimension));

    return length(_dataHost(dpMatrix), dimension);
}

// Returns the overall length of the underlying vector of the hosted matrix.
template <typename TValue, typename TMatrixSpec>
inline typename Size<DPMatrix_<TValue, TMatrixSpec> const>::Type
length(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix)
{
    return length(_dataHost(dpMatrix));  // Note that even if the dimensional lengths are set but the matrix was not resized
    // this function returns 0 or the previous length of the host before the resize.
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
inline void
clear(DPMatrix_<TValue, TMatrixSpec> & dpMatrix)
{
    clear(_dataLengths(dpMatrix));
    resize(_dataLengths(dpMatrix), 2, 0);
    clear(_dataFactors(dpMatrix));
    resize(_dataFactors(dpMatrix), 2, 0);
    _dataFactors(dpMatrix)[DPMatrixDimension_::VERTICAL] = 1u;
    clear(host(dpMatrix));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
inline bool
empty(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix)
{
    return empty(host(dpMatrix));
}

// ----------------------------------------------------------------------------
// Function setLength()
// ----------------------------------------------------------------------------

// Sets the new length of a given dimension.
template <typename TValue, typename TMatrixSpec, typename TSize>
inline void
setLength(DPMatrix_<TValue, TMatrixSpec> & dpMatrix,
          unsigned int dimension,
          TSize const & newLength)
{
    SEQAN_ASSERT(_checkCorrectDimension(dimension));
    setLength(_dataHost(dpMatrix), dimension, newLength);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

// Resizes the matrix. Note, the lengths of the dimensions have to be set before.
template <typename TValue, typename TMatrixSpec>
inline void
resize(DPMatrix_<TValue, TMatrixSpec> & dpMatrix)
{
    resize(_dataHost(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline void
resize(DPMatrix_<TValue, TMatrixSpec> & dpMatrix,
       TValue const & fillValue)
{
    resize(_dataHost(dpMatrix), fillValue);
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec>, Standard const>::Type
begin(DPMatrix_<TValue, TMatrixSpec> & dpMatrix, Standard const)
{
    return begin(host(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec> const, Standard const>::Type
begin(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix, Standard const)
{
    return begin(host(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec>, Rooted const>::Type
begin(DPMatrix_<TValue, TMatrixSpec> & dpMatrix, Rooted const)
{
    return begin(_dataHost(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec> const, Rooted const>::Type
begin(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix, Rooted const)
{
    return begin(_dataHost(dpMatrix));
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TValue, typename TMatrixSpec>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec>, Standard const>::Type
end(DPMatrix_<TValue, TMatrixSpec> & dpMatrix, Standard const)
{
    return end(host(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec> const, Standard const>::Type
end(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix, Standard const)
{
    return end(host(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec>, Rooted const>::Type
end(DPMatrix_<TValue, TMatrixSpec> & dpMatrix, Rooted const)
{
    return end(_dataHost(dpMatrix));
}

template <typename TValue, typename TMatrixSpec>
inline typename Iterator<DPMatrix_<TValue, TMatrixSpec> const, Rooted const>::Type
end(DPMatrix_<TValue, TMatrixSpec> const & dpMatrix, Rooted const)
{
    return end(_dataHost(dpMatrix));
}

// ----------------------------------------------------------------------------
// Function coordinate()
// ----------------------------------------------------------------------------

// Returns the coordinate of a host positio in a given dimension.
template <typename TValue, typename TPosition>
inline typename Position<DPMatrix_<TValue, FullDPMatrix> >::Type
coordinate(DPMatrix_<TValue, FullDPMatrix> const & dpMatrix,
           TPosition hostPos,
           typename DPMatrixDimension_::TValue dimension)
{
    return coordinate(_dataHost(dpMatrix), hostPos, dimension);
}

} // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_DP_MATRIX_H_
