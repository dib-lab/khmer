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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Implementation of the adaptor iterator that makes a non-rooted iterator
// rooted by adding a pointer to the original container.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_ADAPTOR_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_ADAPTOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TValue>
inline size_t length(TValue const * me);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TIterator, typename TSpec = Default>
struct AdaptorIterator;

//////////////////////////////////////////////////////////////////////////////
// Adaptor Iterator
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AdaptorIterator
 * @extends Iter
 * @headerfile <seqan/basic.h>
 * @brief Adapts iterators to RootedIteratorConcept.
 *
 * @signature template <typename TContainer, typename TIterator[, typename TSpec]>
 *            class Iter;
 *
 * @tparam TContainer The container to iterate over.
 * @tparam TIterator  The iterator type to wrap.
 * @tparam TSpec      Specialization tag for AdaptorIterator.  Defaults to <tt>Default</tt>.
 */

/*!
 * @fn AdaptorIterator::Iter
 * @brief Constructor
 *
 * @signature Iter::Iter();
 * @signature Iter::Iter(iter);
 * @signature Iter::Iter(container[, iterator]);
 *
 * @param[in] iter      Other AdaptorIterator to copy from.
 * @param[in] container The container to create an container of.
 * @param[in] iterator  The <tt>Standard</tt> iterator to wrap, defaults to beginning of <tt>container</tt>.
 */

template <typename TContainer, typename TIterator, typename TSpec>
class Iter<TContainer, AdaptorIterator<TIterator, TSpec> >
{
public:
    typename Pointer_<TContainer>::Type data_container;
    TIterator data_iterator;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    SEQAN_HOST_DEVICE
    Iter() : data_container()
    {
        data_iterator = TIterator();
    }

    /*
    // TODO(doering): welches "begin" zur initialisierung von "data_iterator" aufrufen?
    // TODO(holtgrew): Probably call begin(, Standard())...
    Iter(typename Parameter_<TContainer>::Type container_)
        : data_container(_toPointer(container_)),
          data_iterator(begin(container_))
    {
        SEQAN_CHECKPOINT;
    }
    */

    SEQAN_HOST_DEVICE
    Iter(typename Parameter_<TContainer>::Type container_, TIterator it_)
            : data_container(_toPointer(container_)),
              data_iterator(it_)
    {
        SEQAN_CHECKPOINT;
    }

    SEQAN_HOST_DEVICE
    Iter(Iter const & other_)
            : data_container(other_.data_container),
              data_iterator(other_.data_iterator)
    {
        SEQAN_CHECKPOINT;
    }

    // TODO(holtgrew): Use this technique to the other Iter specializations.
    template <typename TContainer_, typename TIterator2>
    SEQAN_HOST_DEVICE
    Iter(Iter<TContainer_, AdaptorIterator<TIterator2, TSpec> > const & other,
         SEQAN_CTOR_ENABLE_IF(IsSameType<TContainer, TContainer_ const>)) :
            data_container(other.data_container), data_iterator(other.data_iterator)
    {
        ignoreUnusedVariableWarning(dummy);\
    }

    // ------------------------------------------------------------------------
    // Assignment Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    SEQAN_HOST_DEVICE
    Iter &
    operator=(Iter const & other_)
    {
        data_container = other_.data_container;
        data_iterator = other_.data_iterator;
        return *this;
    }

    template <typename TContainer_>
    SEQAN_FUNC_ENABLE_IF(IsSameType<TContainer, TContainer_ const>, Iter &)
    SEQAN_HOST_DEVICE
    operator=(Iter<TContainer_, AdaptorIterator<TIterator, TSpec> > const & other_)
    {
        data_container = other_.data_container;
        data_iterator = other_.data_iterator;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Pointer Operators;  Have to be defined within class.
    // ------------------------------------------------------------------------

    // For chaining behaviour of operator->(), see http://stackoverflow.com/a/8782794/84349

    SEQAN_HOST_DEVICE
    TIterator &
    operator->()
    {
        return data_iterator;
    }

    SEQAN_HOST_DEVICE
    TIterator const &
    operator->() const
    {
        return data_iterator;
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    SEQAN_HOST_DEVICE
    operator TIterator() const
    {
        SEQAN_CHECKPOINT;
        return data_iterator;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec>
struct GetValue<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > > :
    GetValue<TIterator> {};

//template <typename TContainer, typename TIterator, typename TSpec>
//struct GetValue<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const> :
//    GetValue<TIterator const> {};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

//template <typename TContainer, typename TIterator, typename TSpec>
//struct Value<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > > :
//    GetValue<TIterator> {};
//
//template <typename TContainer, typename TIterator, typename TSpec>
//struct Value<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const> :
//    GetValue<TIterator const> {};

// ----------------------------------------------------------------------------
// Metafunction IteratorDefaultImp_
// ----------------------------------------------------------------------------

template <typename T>
struct IteratorDefaultImp_<T, Rooted>
{
    typedef typename Iterator<T, Standard>::Type TStandardIterator;
    typedef Iter<T, AdaptorIterator<TStandardIterator> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction Chunk
// ----------------------------------------------------------------------------

// Chunk interface for rooted iterators
template <typename TContainer, typename TValue, typename TSpec>
struct Chunk<Iter<TContainer, AdaptorIterator<TValue*, TSpec> > >
{
    typedef Range<TValue*> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

//template <typename TContainer, typename TIterator, typename TSpec>
//inline SEQAN_HOST_DEVICE typename Parameter_<TContainer>::Type
//container(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
//{
//    SEQAN_CHECKPOINT;
//    return _toParameter<TContainer>(me.data_container);
//}
//
//template <typename TContainer, typename TIterator, typename TSpec>
//inline SEQAN_HOST_DEVICE typename Parameter_<TContainer>::Type
//container(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
//{
//    SEQAN_CHECKPOINT;
//    return _toParameter<TContainer>(me.data_container);
//}

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE TContainer &
container(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
    return _referenceCast<TContainer &>(me.data_container);
}

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE TContainer &
container(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
    return container(const_cast<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &>(me));
}

// ----------------------------------------------------------------------------
// Function setContainer()
// ----------------------------------------------------------------------------

/*!
 * @fn AdaptorIterator#setContainer
 * @brief Set container of an AdaptorIterator.
 *
 * @signature void setContainer(iter, container);
 *
 * @param[in,out] iter      The AdaptorIterator to set the container to.
 * @param[in]     container The container to set for the AdaptorIterator.
 */

// TODO(holtgrew): Also defined in index module, change documentation?
template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE void
setContainer(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
             typename Parameter_<TContainer>::Type container_)
{
    me.data_container = _toPointer(container_);
    setPosition(me, 0);
}

// ----------------------------------------------------------------------------
// Function hostIterator()
// ----------------------------------------------------------------------------

/*!
 * @fn AdaptorIterator#hostIterator
 * @brief Return the host iterator of an AdaptorIterator.
 *
 * @signature TIter hostIterator(iter);
 *
 * @param[in] iter The AdaptorIterator to get the iterator for.
 *
 * @return TIter The wrapped iterator.
 */

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE TIterator &
hostIterator(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
    SEQAN_CHECKPOINT;
    return me.data_iterator;
}

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE TIterator const &
hostIterator(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
    SEQAN_CHECKPOINT;
    return me.data_iterator;
}

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>::Type
position(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
    SEQAN_CHECKPOINT;
    return hostIterator(me) - begin(container(me), Standard());
}

template <typename TContainer, typename TIterator, typename TSpec, typename TContainer2>
inline SEQAN_HOST_DEVICE typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>::Type
position(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me,
         TContainer2 const &)
{
    SEQAN_CHECKPOINT;
    return hostIterator(me) - begin(container(me), Standard());
}

// ----------------------------------------------------------------------------
// Function setPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn AdaptorIterator#setPosition
 * @brief Set position of AdaptorIterator.
 *
 * @signature void setPosition(iter, pos);
 *
 * @param[in,out] iter The AdaptorIterator to set the position for.
 * @param[in]     pos  The position to set.
 */

template <typename TContainer, typename TIterator, typename TSpec, typename TPosition>
inline SEQAN_HOST_DEVICE void
setPosition(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
            TPosition pos_)
{
    SEQAN_CHECKPOINT;
    hostIterator(me) = begin(container(me), Standard()) + pos_;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE typename Reference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type
value(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
    SEQAN_CHECKPOINT;
    return value(hostIterator(me));
}

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE typename Reference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>::Type
value(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
    SEQAN_CHECKPOINT;
    return value(hostIterator(me));
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE typename GetValue<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type
getValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
    SEQAN_CHECKPOINT;
    return getValue(hostIterator(me));
}

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE typename GetValue<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const>::Type
getValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
    SEQAN_CHECKPOINT;
    return getValue(hostIterator(me));
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline SEQAN_HOST_DEVICE void
assignValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    assignValue(hostIterator(me), _value);
}

template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline SEQAN_HOST_DEVICE void
assignValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    assignValue(hostIterator(me), _value);
}

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline SEQAN_HOST_DEVICE void
moveValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    moveValue(hostIterator(me), _value);
}
template <typename TContainer, typename TIterator, typename TSpec, typename TValue>
inline SEQAN_HOST_DEVICE void
moveValue(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    moveValue(hostIterator(me), _value);
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE bool
operator==(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
           Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) == hostIterator(right);
}

// TODO(weese:) Why would we need IterComplementConst here? Disabled it.

//template <typename TContainer, typename TIterator, typename TSpec>
//inline SEQAN_HOST_DEVICE bool
//operator==(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
//           typename IterComplementConst<TIterator>::Type const & right)
//{
//    SEQAN_CHECKPOINT;
//    return hostIterator(left) == right;
//}

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE bool
operator==(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
           TIterator const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) == right;
}

// TODO(weese:) Why would we need IterComplementConst here? Disabled it.

//template <typename TContainer, typename TIterator, typename TSpec>
//inline SEQAN_HOST_DEVICE bool
//operator==(typename IterComplementConst<TIterator>::Type const & left,
//             Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
//{
//    SEQAN_CHECKPOINT;
//    return left == hostIterator(right);
//}

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE bool
operator==(TIterator const & left,
           Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
    SEQAN_CHECKPOINT;
    return left == hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE bool
operator!=(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
           Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) != hostIterator(right);
}

// TODO(weese:) Why would we need IterComplementConst here? Disabled it.

//template <typename TContainer, typename TIterator, typename TSpec>
//inline SEQAN_HOST_DEVICE bool
//operator!=(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
//           typename IterComplementConst<TIterator>::Type const & right)
//{
//    SEQAN_CHECKPOINT;
//    return hostIterator(left) != right;
//}

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE bool
operator!=(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
           TIterator const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) != right;
}

// TODO(weese:) Why would we need IterComplementConst here? Disabled it.

//template <typename TContainer, typename TIterator, typename TSpec>
//inline SEQAN_HOST_DEVICE bool
//operator!=(typename IterComplementConst<TIterator>::Type const & left,
//           Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
//{
//    SEQAN_CHECKPOINT;
//    return left != hostIterator(right);
//}

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE bool
operator!=(TIterator const & left,
           Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
    SEQAN_CHECKPOINT;
    return left != hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE void
goNext(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
    SEQAN_CHECKPOINT;
    goNext(hostIterator(me));
}

// ----------------------------------------------------------------------------
// Function goPrevious()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE void
goPrevious(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
{
    SEQAN_CHECKPOINT;
    goPrevious(hostIterator(me));
}

// ----------------------------------------------------------------------------
// Function operator+()
// ----------------------------------------------------------------------------

// TODO(holtgrew): It should be possible to have default implementations of these operators for all but a core set.

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >
operator+(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
          TIntegral right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) + right);
}

// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >
operator+(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
          int right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) + right);
}

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >
operator+(TIntegral left,
          Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(right), hostIterator(right) + left);
}

// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >
operator+(int left,
          Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(right), hostIterator(right) + left);
}

// ----------------------------------------------------------------------------
// Function operator+=()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator+=(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
           TIntegral right)
{
    SEQAN_CHECKPOINT;
    hostIterator(left) += right;
    return left;
}

// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator+=(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
           int right)
{
    SEQAN_CHECKPOINT;
    hostIterator(left) += right;
    return left;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >
operator-(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
          TIntegral right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) - right);
}

// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> >
operator-(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
          int right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, AdaptorIterator<TIterator, TSpec> >(container(left), hostIterator(left) - right);
}

// TODO(holtgrew): What about Function.difference()?
template <typename TContainer, typename TIterator, typename TSpec, typename TContainer2, typename TIterator2, typename TSpec2>
inline typename Difference<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type
operator-(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & left,
          Iter<TContainer2, AdaptorIterator<TIterator2, TSpec2> > const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) - hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIterator, typename TSpec, typename TIntegral>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator-=(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
           TIntegral right)
{
    SEQAN_CHECKPOINT;
    hostIterator(left) -= right;
    return left;
}

// for <anonymous enum> types
template <typename TContainer, typename TIterator, typename TSpec>
inline Iter<TContainer, AdaptorIterator<TIterator, TSpec> > &
operator-=(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & left,
           int right)
{
    SEQAN_CHECKPOINT;
    hostIterator(left) -= right;
    return left;
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

//template <typename TContainer, typename TIterator, typename TSpec>
//inline SEQAN_HOST_DEVICE bool
//atEnd(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > & me)
//{
//    SEQAN_CHECKPOINT;
//    return atEnd(me, container(me));
//}
//
template <typename TContainer, typename TIterator, typename TSpec>
inline SEQAN_HOST_DEVICE bool
atEnd(Iter<TContainer, AdaptorIterator<TIterator, TSpec> > const & me)
{
    SEQAN_CHECKPOINT;
    return atEnd(me, container(me));
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

// TODO(holtgrew): What about set()?

// Conversion assignment.
template <typename TTargetContainer, typename TIterator, typename TSpec, typename TSource>
inline SEQAN_HOST_DEVICE void
assign(Iter<TTargetContainer, AdaptorIterator<TIterator, TSpec> > & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    target.data_container = container(source);
    target.data_iterator = begin(container(source)) + position(source);
}

// ----------------------------------------------------------------------------
// Function reserveChunk()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TSize>
inline void reserveChunk(Iter<TContainer, TSpec> &iter, TSize size, Output)
{
    typedef Iter<TContainer, TSpec> TIter;

    TContainer &cont = container(iter);
    typename Position<TIter>::Type pos = position(iter);
    typename Size<TIter>::Type len = length(cont);

    if (pos < len)
        return;

    len += size;
    if (len <= capacity(cont))
        return;

    reserve(cont, len);
    setPosition(iter, pos);
}

template <typename TContainer, typename TSpec, typename TSize>
inline void reserveChunk(Iter<TContainer, TSpec> &, TSize, Input)
{}

// ----------------------------------------------------------------------------
// Function getChunk()
// ----------------------------------------------------------------------------

// AdaptorIterator
template <typename TChunk, typename TContainer, typename TValue, typename TSpec>
inline void
getChunk(TChunk &result, Iter<TContainer, AdaptorIterator<TValue*, TSpec> > &rootedIter, Input)
{
    return assignRange(result, hostIterator(rootedIter), end(container(rootedIter), Standard()));
}

template <typename TChunk, typename TContainer, typename TValue, typename TSpec>
inline void
getChunk(TChunk &result, Iter<TContainer, AdaptorIterator<TValue*, TSpec> > &rootedIter, Output)
{
    TContainer &cont = container(rootedIter);
    return assignRange(result, hostIterator(rootedIter), begin(cont, Standard()) + capacity(cont));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_ADAPTOR_H_
