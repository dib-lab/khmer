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
// Iterator interface with default implementations.
// ==========================================================================

// TODO(holtgrew): Split into iterator_interface.h and iterator_adapt_pointer.h.

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_INTERFACE_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_INTERFACE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @defgroup ContainerIteratorTags Container Iterator Tags
 * @brief Tags for container iterators.
 *
 * The tags <tt>Standard</tt> and <tt>Rooted</tt> can be used for selecting specific iterator types with the
 * @link ContainerConcept#Iterator @endlink metafunction.  Rooted iterators also carry a pointer to the container
 * they are iterating whereas standard iterators do not carry this information.
 *
 * @tag ContainerIteratorTags#Standard
 * @headerfile <seqan/basic.h>
 * @brief Tag for selecting standard iterators.
 * @signature struct Standard_;
 *            typedef Tag<Standard_> Standard;
 *
 * @tag ContainerIteratorTags#Rooted
 * @headerfile <seqan/basic.h>
 * @brief Tag for selecting rooted iterators.
 * @signature struct Rooted_;
 *            typedef Tag<Rooted_> Rooted;
 */

struct Rooted_;
typedef Tag<Rooted_> const Rooted;

struct Standard_;
typedef Tag<Standard_> const Standard;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultIteratorSpec
// ----------------------------------------------------------------------------

/*!
 * @mfn ContainerConcept#DefaultIteratorSpec
 * @brief Returns the default iterator specialization.
 *
 * @signature DefaultIteratorSpec<TContainer>::Type
 *
 * @tparam TContainer The Container type to query.
 * @return Type       The iterator specialization tag type.
 *
 * Used by @link ContainerConcept#Iterator @endlink to select the default value for <tt>TSpec</tt>.
 *
 * @see ContainerConcept#Iterator
 */

template <typename T>
struct DefaultIteratorSpec
{
    typedef Standard Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultGetIteratorSpec
// ----------------------------------------------------------------------------

/*!
 * @mfn ContainerConcept#DefaultGetIteratorSpec
 * @brief Returns the default iterator specialization for functions.
 *
 * @signature DefaultGetIteratorSpec<TContainer>::Type
 *
 * @tparam TContainer The Container type to query.
 * @return Type       The iterator specialization tag type.
 *
 * Used by functions such as @link ContainerConcept#begin @endlink and @link ContainerConcept#end @endlink for the <tt>TSpec</tt>
 * parameter.
 *
 * @see ContainerConcept#Iterator
 */

template <typename T>
struct DefaultGetIteratorSpec
{
    typedef Rooted Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename T, typename TSpec>
struct IteratorDefaultImp_;

// We use plain pointers as standard iterators.
template <typename T>
struct IteratorDefaultImp_<T, Standard>
{
    typedef typename Value<T>::Type * Type;
};

// (weese): This definition is important and defines default const-iterators. Don't remove.
//          However, there are different places where const-correctness is broken that must be fixed before we can uncomment this

template <typename T>
struct IteratorDefaultImp_<T const, Standard>
{
    typedef typename Value<T>::Type const * Type;
};

//IteratorDefaultImp_<T, Rooted> is implemented in basic_iterator_adaptor.h

// TODO(weese): Mmh. What was the reason to introduce the helper struct IteratorDefaultImp_ instead of directly defining it here.
//              Aah. I guess in to allow to specialize Iterator only in the first template argument. However, right now it is always
//              specialized for both the first and second argument everywhere in the code.
template <typename T, typename TSpec = typename DefaultIteratorSpec<T>::Type>
struct Iterator : IteratorDefaultImp_<T, TSpec>
{
};

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

// TODO(holtgrew): Remove the default implementation; anti-auto-sequence. Also, using plain pointers for strings does not work any more. Will probably only work for rooted/adaptor/positional iterators. Same below.

template <typename T>
struct Container
{
    typedef T Type;
};

// ============================================================================
// Functions
// ============================================================================

// ---------------------------------------------------------------------------
// Function value()
// ---------------------------------------------------------------------------

template <typename T>
inline typename Reference<T>::Type
value(T & me)
{
    SEQAN_CHECKPOINT;
    return *me;
}

template <typename T>
inline typename Reference<T const>::Type
value(T const & me)
{
    SEQAN_CHECKPOINT;
    return *me;
}

// ---------------------------------------------------------------------------
// Function getValue()
// ---------------------------------------------------------------------------

template <typename T>
inline typename GetValue<T>::Type
getValue(T & me)
{
    SEQAN_CHECKPOINT;
    return value(me);
}

template <typename T>
inline typename GetValue<T const>::Type
getValue(T const & me)
{
    SEQAN_CHECKPOINT;
    return value(me);
}

// ---------------------------------------------------------------------------
// Function toGetValue()
// ---------------------------------------------------------------------------

//Nimmt eine Reference und macht daraus einen GetValue
// TODO(doering):toGetValue()

// ---------------------------------------------------------------------------
// Function assignValue()
// ---------------------------------------------------------------------------

template <typename T, typename TValue>
inline void
assignValue(T & me,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    assign(value(me), _value);
}

//const version for iterators as targets
template <typename T, typename TValue>
inline void
assignValue(T const & me,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    assign(value(me), _value);
}

// ---------------------------------------------------------------------------
// Function moveValue()
// ---------------------------------------------------------------------------

/*!
 * @fn OutputIteratorConcept#moveValue
 * @headerfile <seqan/sequence.h>
 * @brief Move a value of a container to a given position.
 *
 * @signature void moveValue(container, pos, value);
 *
 * @param[in,out] container The container to manipulate.
 * @param[in]     pos       The position of the item in the container to manipulate.
 * @param[in,out] value     The value to move to <tt>container[pos]</tt>.
 */

template <typename T, typename TValue>
inline void
moveValue(T & me,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    move(value(me), _value);
}

//const version for iterators as targets
template <typename T, typename TValue>
inline void
moveValue(T const & me,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    move(value(me), _value);
}

// ---------------------------------------------------------------------------
// Function setValue()
// ---------------------------------------------------------------------------

template <typename T, typename TValue>
inline void
setValue(T * & ptr,
         TValue & _value)
{
    SEQAN_CHECKPOINT;
    ptr = &_value;
}

//const version for iterators as targets
template <typename T, typename TValue>
inline void
setValue(T const * & ptr,
         TValue const & _value)
{
    SEQAN_CHECKPOINT;
    ptr = &_value;
}

// ---------------------------------------------------------------------------
// Function container()
// ---------------------------------------------------------------------------

template <typename T>
inline
typename Container<T>::Type &
container(T me)
{
    // TODO(holtgrew): Default implementation with auto-sequences, remove?
    SEQAN_CHECKPOINT;
    return me;
}

// ---------------------------------------------------------------------------
// Function position()
// ---------------------------------------------------------------------------

template <typename T>
inline typename Position<T>::Type
position(T * /*me*/)
{
    // TODO(holtgrew): Default implementation with auto-sequences, remove?
    SEQAN_CHECKPOINT;
    return 0;
}

template <typename TContainer, typename TIterator>
inline SEQAN_HOST_DEVICE typename Position<TContainer>::Type
position(TIterator const & it,
         TContainer const & me)
{
    SEQAN_CHECKPOINT;
    return it - begin(me, Standard());
}

// ---------------------------------------------------------------------------
// Function atBegin()
// ---------------------------------------------------------------------------

// TODO(doering): Was, wenn der Container leer ist?

template <typename T, typename TContainer>
inline SEQAN_HOST_DEVICE bool
atBegin(T const & it, TContainer const & cont)
{
    SEQAN_CHECKPOINT;
    return it == begin(cont, Standard());
}

template <typename T, typename TContainer>
inline SEQAN_HOST_DEVICE bool
atBegin(T const & it, TContainer & cont)
{
    SEQAN_CHECKPOINT;
    return it == begin(cont, Standard());
}

template <typename T, typename TContainer>
inline SEQAN_HOST_DEVICE bool
atBegin(T & it, TContainer & cont)
{
    SEQAN_CHECKPOINT;
    return it == begin(cont, Standard());
}

template <typename T, typename TContainer>
inline SEQAN_HOST_DEVICE bool
atBegin(T & it, TContainer const & cont)
{
    SEQAN_CHECKPOINT;
    return it == begin(cont, Standard());
}

template <typename T>
inline SEQAN_HOST_DEVICE bool
atBegin(T & it)
{
    return atBegin(it, container(it));
}

template <typename T>
inline SEQAN_HOST_DEVICE bool
atBegin(T const & it)
{
    return atBegin(it, container(it));
}

// ---------------------------------------------------------------------------
// Function atEnd()
// ---------------------------------------------------------------------------

template <typename T, typename TContainer>
inline SEQAN_HOST_DEVICE bool
atEnd(T & it,
      TContainer const & cont)
{
    SEQAN_CHECKPOINT;
    return it == end(cont, Standard());
}

template <typename T, typename TContainer>
inline SEQAN_HOST_DEVICE bool
atEnd(T const & it,
      TContainer const & cont)
{
    SEQAN_CHECKPOINT;
    return it == end(cont, Standard());
}

template <typename T, typename TContainer>
inline SEQAN_HOST_DEVICE bool
atEnd(T & it,
      TContainer & cont)
{
    SEQAN_CHECKPOINT;
    return it == end(cont, Standard());
}

template <typename T, typename TContainer>
inline SEQAN_HOST_DEVICE bool
atEnd(T const & it,
      TContainer & cont)
{
    SEQAN_CHECKPOINT;
    return it == end(cont, Standard());
}

//template <typename T>
//inline SEQAN_HOST_DEVICE bool
//atEnd(T & it)
//{
//    SEQAN_CHECKPOINT;
//    return atEnd(it, container(it));
//}

template <typename T>
inline SEQAN_HOST_DEVICE bool
atEnd(T const & it)
{
    SEQAN_CHECKPOINT;
    return atEnd(it, container(it));
}

// ---------------------------------------------------------------------------
// Function goBegin()
// ---------------------------------------------------------------------------

template <typename TIterator, typename TContainer>
inline SEQAN_HOST_DEVICE void
goBegin(TIterator & it,
        TContainer & container)
{
    it = begin(container);
}

// template <typename TIterator, typename TContainer>
// inline void
// goBegin(TIterator & it,
//         TContainer const & container)
// {
//     it = begin(container);
// }

template <typename TIterator>
inline SEQAN_HOST_DEVICE void
goBegin(TIterator & it)
{
    typename Parameter_<typename Container<TIterator>::Type>::Type tmpContainer = container(it);
    goBegin(it, tmpContainer);
}

// ---------------------------------------------------------------------------
// Function goEnd()
// ---------------------------------------------------------------------------

template <typename TIterator, typename TContainer>
inline SEQAN_HOST_DEVICE void
goEnd(TIterator & it,
      TContainer & container)
{
    SEQAN_CHECKPOINT;
    it = end(container);
}

template <typename TIterator, typename TContainer>
inline SEQAN_HOST_DEVICE void
goEnd(TIterator & it,
      TContainer const & container)
{
    SEQAN_CHECKPOINT;
    it = end(container);
}

template <typename TIterator>
inline SEQAN_HOST_DEVICE void
goEnd(TIterator & it)
{
    SEQAN_CHECKPOINT;
    goEnd(it, container(it));
}

// ---------------------------------------------------------------------------
// Function goNext()
// ---------------------------------------------------------------------------

template <typename TIterator>
inline SEQAN_HOST_DEVICE void
goNext(TIterator & it)
{
    SEQAN_CHECKPOINT;
    ++it;
}

// ---------------------------------------------------------------------------
// Function goFurther()
// ---------------------------------------------------------------------------

template <typename TIterator, typename TDiff>
inline SEQAN_HOST_DEVICE void
goFurther(TIterator & it,
          TDiff steps)
{   // return distance type from arbitrary argument
    it += steps;
}

// ---------------------------------------------------------------------------
// Function goPrevious()
// ---------------------------------------------------------------------------

template <typename TIterator>
inline SEQAN_HOST_DEVICE void
goPrevious(TIterator & it)
{
    SEQAN_CHECKPOINT;
    --it;
}

// ---------------------------------------------------------------------------
// Function difference()
// ---------------------------------------------------------------------------

template <typename TIterator>
inline SEQAN_HOST_DEVICE
typename Difference<TIterator>::Type
difference(TIterator const & begin,
           TIterator const & end)
{
    SEQAN_CHECKPOINT;
    return end - begin;
}

// ---------------------------------------------------------------------------
// Function goNil()
// ---------------------------------------------------------------------------

template <typename TIterator>
inline void
goNil(TIterator & me)
{
    SEQAN_CHECKPOINT;
    me = TIterator();
}

template <typename TIterator>
inline void
goNil(TIterator * & me)
{
    SEQAN_CHECKPOINT;
    me = 0;
}

// ---------------------------------------------------------------------------
// Function atNil()
// ---------------------------------------------------------------------------

template <typename TIterator>
inline bool
atNil(TIterator & me)
{
    SEQAN_CHECKPOINT;
    return me == TIterator();
}

template <typename TIterator>
inline bool
atNil(TIterator * me)
{
    SEQAN_CHECKPOINT;
    return me == 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_INTERFACE_H_
