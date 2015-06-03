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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Counting iterator implementation.
// ==========================================================================

#ifndef SEQAN_BASIC_ITERATOR_COUNTING_H_
#define SEQAN_BASIC_ITERATOR_COUNTING_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

template <typename TSpec = void>
struct CountingIteratorImpl_;

typedef CountingIteratorImpl_<void> CountingIterator;

// ============================================================================
// Classes
// ============================================================================

template <typename TSpec, typename TIncrementable>
class Iter<TIncrementable, CountingIteratorImpl_<TSpec> >
{
public:
    TIncrementable data_position;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    Iter(TIncrementable position = 0) :
        data_position(position)
    {}

    template <typename TOther>
    Iter(TOther position) :
        data_position(position)
    {}

    Iter(Iter const & other) :
        data_position(other.data_position)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TSpec, typename TIncrementable>
struct Size<Iter<TIncrementable, CountingIteratorImpl_<TSpec> > >
{
    typedef TIncrementable Type;
};

template <typename TSpec, typename TIncrementable>
struct Position<Iter<TIncrementable, CountingIteratorImpl_<TSpec> > >
{
    typedef TIncrementable Type;
};

template <typename TSpec, typename TIncrementable>
struct Reference<Iter<TIncrementable, CountingIteratorImpl_<TSpec> > >
{
    typedef TIncrementable Type;
};

template <typename TSpec, typename TIncrementable>
struct Difference<Iter<TIncrementable, CountingIteratorImpl_<TSpec> > >
{
    typedef typename MakeSigned<TIncrementable>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIncrementable>
inline typename Position<Iter<TIncrementable, CountingIteratorImpl_<TSpec> > >::Type &
position(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > & me)
{
    return me.data_position;
}

template <typename TSpec, typename TIncrementable>
inline typename Position<Iter<TIncrementable, CountingIteratorImpl_<TSpec> > >::Type const &
position(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & me)
{
    return me.data_position;
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TSpec, typename TPosition>
inline void
setPosition(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > & me, TPosition position_)
{
    me.data_position = position_;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIncrementable>
inline typename Reference<Iter<TIncrementable, CountingIteratorImpl_<TSpec> > >::Type
value(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > & me)
{
    return position(me);
}

template <typename TSpec, typename TIncrementable>
inline typename Reference<Iter<TIncrementable, CountingIteratorImpl_<TSpec> > >::Type
value(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & me)
{
    return position(me);
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TSpec, typename TValue>
inline void
assignValue(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > & me, TValue _value)
{
    setPosition(me, _value);
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIncrementable>
inline bool
operator==(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & left,
           Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & right)
{
    return position(left) == position(right);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIncrementable>
inline bool
operator!=(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & left,
           Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & right)
{
    return position(left) != position(right);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIncrementable>
inline bool
operator<(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & left,
          Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & right)
{
    return position(left) < position(right);
}

template <typename TSpec, typename TIncrementable>
inline bool
operator>(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & left,
          Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & right)
{
    return position(left) > position(right);
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIncrementable>
inline bool
operator<=(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & left,
           Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & right)
{
    return position(left) <= position(right);
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIncrementable>
inline bool
operator>=(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & left,
           Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & right)
{
    return position(left) >= position(right);
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIncrementable>
inline void
goNext(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > & me)
{
    setPosition(me, position(me) + 1);
}

// ----------------------------------------------------------------------------
// Function goPrevious()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TIncrementable>
inline void
goPrevious(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > & me)
{
    setPosition(me, position(me) - 1);
}

// ----------------------------------------------------------------------------
// Function operator+()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TSpec, typename TIntegral>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TIntegral> >, Iter<TIncrementable, CountingIteratorImpl_<TSpec> >)
operator+(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & left, TIntegral right)
{
    Iter<TIncrementable, CountingIteratorImpl_<TSpec> > tmp(left);
    setPosition(tmp, position(left) + right);
    return tmp;
}

template <typename TIncrementable, typename TSpec, typename TIntegral>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TIntegral> >, Iter<TIncrementable, CountingIteratorImpl_<TSpec> >)
operator+(TIntegral left, Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & right)
{
    Iter<TIncrementable, CountingIteratorImpl_<TSpec> > tmp(right);
    setPosition(tmp, position(right) + left);
    return tmp;
}

// ----------------------------------------------------------------------------
// Function operator+=()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TSpec, typename TIntegral>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TIntegral> >, Iter<TIncrementable, CountingIteratorImpl_<TSpec> > &)
operator+=(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > & left, TIntegral right)
{
    setPosition(left, position(left) + right);
    return left;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TSpec, typename TIntegral>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TIntegral> >, Iter<TIncrementable, CountingIteratorImpl_<TSpec> >)
operator-(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & left, TIntegral right)
{
    Iter<TIncrementable, CountingIteratorImpl_<TSpec> > tmp(left);
    setPosition(tmp, position(left) - right);
    return tmp;
}

template <typename TSpec, typename TIncrementable>
inline typename Difference<Iter<TIncrementable, CountingIteratorImpl_<TSpec> > >::Type
operator-(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & left,
          Iter<TIncrementable, CountingIteratorImpl_<TSpec> > const & right)
{
    return position(left) - position(right);
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TIncrementable, typename TSpec, typename TIntegral>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TIntegral> >, Iter<TIncrementable, CountingIteratorImpl_<TSpec> > &)
operator-=(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > & left, TIntegral right)
{
    setPosition(left, position(left) - right);
    return left;
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

// Conversion assignment.
template <typename TIncrementable, typename TSpec, typename TSource>
inline void
assign(Iter<TIncrementable, CountingIteratorImpl_<TSpec> > & target, TSource const & source)
{
    setPosition(target, position(source));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_ITERATOR_COUNTING_H_
