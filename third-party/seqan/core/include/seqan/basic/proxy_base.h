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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Proxy base class definition.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_PROXY_BASE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_PROXY_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.Proxy:
..cat:Basic
..summary:Emulates object of another class.
..signature:Proxy<TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
..remarks.text:Use @Metafunction.Value@ to get the emulated type.
An instance of $Proxy$ behaves like an object of its value type.
$Proxy$ can be used as reference type (see @Metafunction.Reference@).
..remarks.text:Note that functions that are both general and specialized for
the value type should be specialized for $Proxy<TSpec>$ too,
since otherwise the general version will be called.
..include:seqan/basic.h
*/

template <typename TSpec>
class Proxy;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

///.Metafunction.Spec.param.T.type:Class.Proxy
///.Metafunction.Spec.class:Class.Proxy

template <typename TSpec>
struct Spec<Proxy<TSpec> >
{
    typedef TSpec Type;
};

template <typename TSpec>
struct Spec<Proxy<TSpec> const>
{
    typedef TSpec Type;
};

// ----------------------------------------------------------------------------
// Metafunction CompareType
// ----------------------------------------------------------------------------

template <typename TSpec, typename T>
struct CompareType<Proxy<TSpec>, T>
{
    typedef typename Value<Proxy<TSpec> >::Type TValue;
    typedef typename RemoveConst_<TValue>::Type TValue_NoConst;
    typedef typename CompareType<TValue_NoConst, T>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function convertImpl()
// ----------------------------------------------------------------------------

// TODO(holtgrew): First variant even necessary?

template <typename TTarget, typename T, typename TSpec>
inline typename Convert<TTarget, Proxy<TSpec> >::Type
convertImpl(Convert<TTarget, T> const,
            Proxy<TSpec> & source)
{
    return convert<TTarget>(getValue(source));
}

template <typename TTarget, typename T, typename TSpec>
inline typename Convert<TTarget, Proxy<TSpec> const>::Type
convertImpl(Convert<TTarget, T> const,
            Proxy<TSpec> const & source)
{
    return convert<TTarget>(getValue(source));
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TRight>
inline bool
operator==(Proxy<TSpec> const & left_,
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator==(TLeft const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator==(Proxy<TLeftSpec> const & left_,
           Proxy<TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TLeftSpec> TLeft;
    typedef Proxy<TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) == convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator==(Proxy<TSpec> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
    return convert<TAccessor>(left_) == convert<TAccessor>(right_);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TRight>
inline bool
operator!=(Proxy<TSpec> const & left_,
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator!=(TLeft const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator!=(Proxy<TLeftSpec> const & left_,
           Proxy<TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TLeftSpec> TLeft;
    typedef Proxy<TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) != convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator!=(Proxy<TSpec> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
    return convert<TAccessor>(left_) != convert<TAccessor>(right_);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TRight>
inline bool
operator<(Proxy<TSpec> const & left_,
          TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator<(TLeft const & left_,
          Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator<(Proxy<TLeftSpec> const & left_,
          Proxy<TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TLeftSpec> TLeft;
    typedef Proxy<TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) < convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator<(Proxy<TSpec> const & left_,
          Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
    return convert<TAccessor>(left_) < convert<TAccessor>(right_);
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TRight>
inline bool
operator<=(Proxy<TSpec> const & left_,
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator<=(TLeft const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator<=(Proxy<TLeftSpec> const & left_,
           Proxy<TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TLeftSpec> TLeft;
    typedef Proxy<TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) <= convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator<=(Proxy<TSpec> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
    return convert<TAccessor>(left_) <= convert<TAccessor>(right_);
}


// ----------------------------------------------------------------------------
// Function operator>()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TRight>
inline bool
operator>(Proxy<TSpec> const & left_,
          TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator>(TLeft const & left_,
          Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator>(Proxy<TLeftSpec> const & left_,
          Proxy<TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TLeftSpec> TLeft;
    typedef Proxy<TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) > convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator>(Proxy<TSpec> const & left_,
          Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
    return convert<TAccessor>(left_) > convert<TAccessor>(right_);
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TRight>
inline bool
operator>=(Proxy<TSpec> const & left_,
           TRight const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TLeft;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeft, typename TSpec>
inline bool
operator>=(TLeft const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TSpec> TRight;
    typedef typename CompareType<TRight, TLeft>::Type TCompareType;
    return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TLeftSpec, typename TRightSpec>
inline bool
operator>=(Proxy<TLeftSpec> const & left_,
           Proxy<TRightSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef Proxy<TLeftSpec> TLeft;
    typedef Proxy<TRightSpec> TRight;
    typedef typename CompareType<TLeft, TRight>::Type TCompareType;
    return convert<TCompareType>(left_) >= convert<TCompareType>(right_);
}

template <typename TSpec>
inline bool
operator>=(Proxy<TSpec> const & left_,
           Proxy<TSpec> const & right_)
{
    SEQAN_CHECKPOINT;
    typedef typename GetValue<Proxy<TSpec> >::Type TAccessor;
    return convert<TAccessor>(left_) >= convert<TAccessor>(right_);
}

// ----------------------------------------------------------------------------
// Function operator>>();  Reading from streams.
// ----------------------------------------------------------------------------

template <typename TStream, typename TSpec>
inline TStream &
operator>>(TStream & strm,
           Proxy<TSpec> & proxy)
{
    typedef Proxy<TSpec> TProxy;
    typedef typename Value<TProxy>::Type TValue;
    TValue temp;
    strm >> temp;
    assignValue(iter(proxy), temp);
    return strm;
}

template <typename TStream, typename TSpec>
inline TStream &
operator>>(TStream & strm,
           Proxy<TSpec> const& proxy)
{
    typedef Proxy<TSpec> TProxy;
    typedef typename Value<TProxy>::Type TValue;
    TValue temp;
    strm >> temp;
    assignValue(iter(proxy), temp);
    return strm;
}

// ----------------------------------------------------------------------------
// Function operator<<();  Writing to streams.
// ----------------------------------------------------------------------------

// TODO(holtgrew): Is the first variant even necessary?

template <typename TStream, typename TSpec>
inline TStream &
operator<<(TStream & strm,
           Proxy<TSpec> & proxy)
{
    strm << getValue(proxy);
    return strm;
}

template <typename TStream, typename TSpec>
inline TStream &
operator<<(TStream & strm,
           Proxy<TSpec> const & proxy)
{
    strm << getValue(proxy);
    return strm;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_PROXY_BASE_H_
