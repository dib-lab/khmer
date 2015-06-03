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
// Positional iterator implementation.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_POSITION_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_POSITION_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

typedef CountingIteratorImpl_<Nothing> PositionIterator;

/*!
 * @class PositionIterator
 * @extends Iter
 * @headerfile <seqan/basic.h>
 * @brief Adapts a position iterator to a rooted iterator.
 *
 * @signature template <typename TContainer>
 *            class Iter<TContainer, PositionIterator>;
 *
 * @tparam TContainer The container to create an iterator for.
 *
 * @section Remarks
 *
 * PositionIterators provide the concept RootedIteratorConcept.
 *
 *
 * @fn PositionIterator::Iter
 * @brief Constructor
 *
 * @signature Iter::Iter();
 * @signature Iter::Iter(other);
 * @signature Iter::Iter(container[, position]);
 *
 * @param[in] other     Other PositionIterator to copy from.
 * @param[in] container A TContainer to get an iterator to.
 * @param[in] position  The position to create the iterator at, defauls to 0.
 */

template <typename TContainer>
class Iter<TContainer, PositionIterator> :
    public Iter<typename Position<TContainer>::Type, CountingIterator>
{
public:
    typedef Iter<typename Position<TContainer>::Type, CountingIterator> TBase;
    typedef typename Position<TContainer>::Type                         TPosition;
    typedef typename Pointer_<TContainer>::Type                         TContainerPointer_;

    TContainerPointer_ data_container;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    Iter() :
        TBase(0),
        data_container(TContainerPointer_(0))
    {}

    Iter(typename Parameter_<TContainer>::Type container_, TPosition position_ = 0) :
        TBase(position_),
        data_container(_toPointer(container_))
    {}

    Iter(Iter const & other_) :
        TBase(static_cast<TBase const &>(other_)),
        data_container(other_.data_container)
    {}

    template <typename TContainer2, typename TSpec2>
    Iter(Iter<TContainer2, TSpec2> const & other_) :
        TBase(position(other_)),
        data_container(_toPointer(container(other_)))
    {}

    // ------------------------------------------------------------------------
    // Pointer Operators;  Have to be defined within class.
    // ------------------------------------------------------------------------

    typename Value<Iter>::Type *
    operator->()
    {
        return &data_container[position(*this)];
    }

    typename Value<Iter>::Type const *
    operator->() const
    {
        return &data_container[position(*this)];
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    // TODO(holtgrew): Should this define a conversion operator to the underlying iterator?
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TContainer>
struct Size<Iter<TContainer, PositionIterator> > : Size<TContainer> {};

template <typename TContainer>
struct Position<Iter<TContainer, PositionIterator> > : Position<TContainer> {};

template <typename TContainer>
struct Reference<Iter<TContainer, PositionIterator> > : Reference<TContainer> {};

template <typename TContainer>
struct Difference<Iter<TContainer, PositionIterator> > : Difference<TContainer> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline typename Parameter_<TContainer>::Type
container(Iter<TContainer, PositionIterator> & me)
{
    return _toParameter<TContainer>(me.data_container);
}

template <typename TContainer>
inline typename Parameter_<TContainer>::Type
container(Iter<TContainer, PositionIterator> const & me)
{
    return _toParameter<TContainer>(me.data_container);
}

// ----------------------------------------------------------------------------
// Function setContainer()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline void
setContainer(Iter<TContainer, PositionIterator> & me, typename Parameter_<TContainer>::Type container_)
{
    typedef Iter<TContainer, PositionIterator> TIter;
    typename Position<TIter>::Type pos = position(me);
    me.data_container = _toPointer(container_);
    setPosition(me, pos);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline typename Reference<Iter<TContainer, PositionIterator> >::Type
value(Iter<TContainer, PositionIterator> & me)
{
    return value(container(me), position(me));
}

template <typename TContainer>
inline typename Reference<Iter<TContainer, PositionIterator> >::Type
value(Iter<TContainer, PositionIterator> const & me)
{
    return value(container(me), position(me));
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TValue>
inline void
assignValue(Iter<TContainer, PositionIterator> & me, TValue const & _value)
{
    assignValue(container(me), position(me), _value);
}

template <typename TContainer, typename TValue>
inline void
assignValue(Iter<TContainer, PositionIterator> const & me, TValue const & _value)
{
    assignValue(container(me), position(me), _value);
}

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): What are these manual forwards for? Can include order fix this?
// manual forwards
template <typename T, typename TValue, typename TPos>
inline void moveValue(T & me, TPos pos, TValue const & _value);
template <typename T, typename TValue, typename TPos>
inline void moveValue(T const & me, TPos pos, TValue const & _value);

template <typename TContainer, typename TValue>
inline void
moveValue(Iter<TContainer, PositionIterator> & me, TValue const & _value)
{
    moveValue(container(me), position(me), _value);
}

template <typename TContainer, typename TValue>
inline void
moveValue(Iter<TContainer, PositionIterator> const & me, TValue const & _value)
{
    moveValue(container(me), position(me), _value);
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

// Conversion assignment.
template <typename TTargetContainer, typename TSource>
inline void
assign(Iter<TTargetContainer, PositionIterator> & target, TSource const & source)
{
    setContainer(target, container(source));
    setPosition(target, position(source));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_POSITION_H_
