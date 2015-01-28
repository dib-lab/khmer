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
// Iterator interface with default implementations.
// ==========================================================================

// TODO(holtgrew): Split into iterator_interface.h and iterator_adapt_pointer.h.

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ITERATOR_INTERFACE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_ITERATOR_INTERFACE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @defgroup ContainerIteratorTags Container Iterator Tags
 *
 * The tags <tt>Standard</tt> and <tt>Rooted</tt> can be used for selecting specific iterator types with the
 * @link Container#Iterator @endlink metafunction.  Rooted iterators also carry a pointer to the container
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

/**
.Tag.Iterator Spec:
..cat:Iteration
..summary:Specifies the kind of an iterator.
..tag.Rooted:Rooted iterator.
...remarks
....text:This iterator implements some more advanced functions like
@Function.container@ and @Function.position@.
....concept:Concept.RootedIteratorConcept
..tag.Standard:Standard conform iterator.
...remarks
....text:Note that standard iterators need not to implement all functions
that are available for rooted iterators.
..remarks.text:The default iterator spec is given by @Metafunction.DefaultIteratorSpec@.
..see:Metafunction.DefaultIteratorSpec
..see:Concept.IteratorAssociatedTypesConcept
..include:seqan/basic.h
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

/**
.Metafunction.DefaultIteratorSpec:
..hidefromindex
..concept:Concept.ContainerConcept
..summary:Specifies default kind of iterator.
..signature:DefaultIteratorSpec<T>::Type
..param.T:Container type for which the default iterator spec is determined.
...type:Concept.ContainerConcept
..returns.param.Type:Iterator spec of $T$.
..see:Metafunction.Iterator
..include:seqan/basic.h
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

/**
.Metafunction.DefaultGetIteratorSpec:
..hidefromindex
..concept:Concept.ContainerConcept
..summary:Specifies default kind of iterator returned by functions.
..signature:DefaultGetIteratorSpec<T>::Type
..param.T:Container type for which the spec is determined.
...concept:Concept.ContainerConcept
..returns.param.Type:Iterator spec of $T$.
..remarks:This metafunction returns the iterator spec of iterators that are returned by functions like
@Function.begin@, @Function.end@, or @Function.iter@.
..see:Metafunction.Iterator
..see:Metafunction.DefaultIteratorSpec
..include:seqan/basic.h
*/

template <typename T>
struct DefaultGetIteratorSpec
{
    typedef Rooted Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

/**
.Metafunction.Iterator:
..concept:Concept.ContainerConcept
..cat:Iteration
..summary:Type of iterator objects that are used to traverse the container.
..signature:Iterator<T, TSpec>::Type
..param.T:Type for which the iterator type is determined.
...concept:Concept.ContainerConcept
...type:Class.Iter
..param.TSpec:Specifies an @Tag.Iterator Spec.iterator spec@.
...default:The default iterator spec is given by @Metafunction.DefaultIteratorSpec@.
..returns.param.Type:Iterator type of $T$.
..remarks.text:Iterators behave like pointers in some respects.
 For example, you can use $*it$ to access the value object the iterator $it$ points to.
 But note that $Iterator<T>::Type$ can differ from $T *$, depending on $T$.
..see:Metafunction.Position
..include:seqan/basic.h
*/

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

/**
.Metafunction.Container:
..class:Class.Iter
..cat:Iteration
..summary:Type of the container given an iterator.
..signature:Container<T>::Type
..param.T:Iterator type.
...type:Class.Iter
...concept:Concept.RootedIteratorConcept
..returns.param.Type:The container type to $T$.
..include:seqan/basic.h
*/

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

/**
.Function.value
..class:Class.Iter
..signature:Reference value(object)
..param.object:An object that holds a value or an iterator that points to a value.
...type:Class.Iter
...concept:Concept.IteratorAssociatedTypesConcept
..include:seqan/basic.h
*/

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

/**
.Function.getValue
..class:Class.Iter
..cat:Iteration
..signature:GetValue getValue(object)
..param.object:An object that holds a value or points to a value.
...type:Class.Iter
...concept:Concept.IteratorAssociatedTypesConcept
..see:Metafunction.GetValue
..include:seqan/basic.h
*/

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

/**
.Function.assignValue
..class:Class.Iter
..cat:Iteration
..summary:Assigns value to item.
..signature:assignValue(object, value)
..param.object:An object that holds a value or points to a value.
...type:Class.Iter
...concept:Concept.BasicOutputIteratorConcept
..param.value:A value that is assigned to the item $object$ holds or points to.
..remarks.text:This function is similar to @Function.assign@.
The difference is, that $assignValue$ just changes a value stored in $object$ or the value $object$ points to,
while @Function.assign@ changes the whole object.
..see:Function.assign
..include:seqan/basic.h
*/

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

/**
.Function.moveValue
..class:Class.Iter
..cat:Iteration
..summary:Assigns value to item.
..signature:moveValue(object, value)
..param.object:An object that holds a value or points to a value.
...type:Class.Iter
...concept:Concept.BasicOutputIteratorConcept
..param.value:A value that is handed over to the item $object$ holds or points to.
..remarks.text:This function is similar to @Function.move@.
The difference is, that $moveValue$ just changes a value stored in $object$ or the value $object$ points to,
while @Function.move@ changes the whole object.
..see:Function.move
..see:Function.assignValue
..include:seqan/basic.h
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

///.Function.setValue.param.object.type:Concept.BasicOutputIteratorConcept
///.Function.setValue.concept:Concept.BasicOutputIteratorConcept

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

/**
.Function.container
..concept:Concept.RootedIteratorConcept
..cat:Iteration
..summary:Container of an iterator.
..signature:Container container(iterator)
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.RootedIteratorConcept
..returns:The container that $iterator$ traverses.
...metafunction:Metafunction.Container
..include:seqan/basic.h
*/

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

/**
.Function.position
..class:Class.Iter
..concept:Concept.ContainerConcept
..summary:Position of an iterator.
..cat:Iteration
..signature:Position position(iterator [, container])
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.RootedRandomAccessIteratorConcept
..param.container:A container.
...concept:Concept.ContainerConcept
...remarks:If $iterator$ implements @Concept.RootedIteratorConcept@, then $container$ is optional.
...remarks:If $container$ is specified, $iterator$ must be a container of $container$.
..returns:The position of the value in the container $iterator$ points to.
...metafunction:Metafunction.Position
..include:seqan/basic.h
*/

template <typename T>
inline typename Position<T>::Type
position(T * /*me*/)
{
    // TODO(holtgrew): Default implementation with auto-sequences, remove?
    SEQAN_CHECKPOINT;
    return 0;
}

template <typename TContainer, typename TIterator>
inline typename Position<TContainer>::Type
position(TIterator const & it,
         TContainer const & me)
{
    SEQAN_CHECKPOINT;
    return it - begin(me, Standard());
}

// ---------------------------------------------------------------------------
// Function atBegin()
// ---------------------------------------------------------------------------

/**
.Function.atBegin
..class:Class.Iter
..concept:Concept.ContainerConcept
..concept:Concept.RootedIteratorConcept
..cat:Iteration
..summary:Determines whether an iterator is at the beginning position.
..signature:bool atBegin(iterator [, container])
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.RootedIteratorConcept
..param.container:Container of $iterator$. (optional)
...remarks.text:If $iterator$ implements @Concept.RootedIteratorConcept@ then $container$ is optional otherwise $container$ is required.
..returns:$true$ if $iterator$ points to the fist item of the container, otherwise $false$.
..see:Function.begin
..include:seqan/basic.h
*/

// TODO(doering): Was, wenn der Container leer ist?

template <typename T, typename TContainer>
inline bool
atBegin(T const & it, TContainer const & cont)
{
    SEQAN_CHECKPOINT;
    return it == begin(cont, Standard());
}

template <typename T, typename TContainer>
inline bool
atBegin(T const & it, TContainer & cont)
{
    SEQAN_CHECKPOINT;
    return it == begin(cont, Standard());
}

template <typename T, typename TContainer>
inline bool
atBegin(T & it, TContainer & cont)
{
    SEQAN_CHECKPOINT;
    return it == begin(cont, Standard());
}

template <typename T, typename TContainer>
inline bool
atBegin(T & it, TContainer const & cont)
{
    SEQAN_CHECKPOINT;
    return it == begin(cont, Standard());
}

template <typename T>
inline bool
atBegin(T const & it)
{
    SEQAN_CHECKPOINT;
    return atBegin(it, container(it));
}

// ---------------------------------------------------------------------------
// Function atEnd()
// ---------------------------------------------------------------------------

/**
.Function.atEnd
..class:Class.Iter
..concept:Concept.ContainerConcept
..concept:Concept.RootedIteratorConcept
..cat:Iteration
..summary:Determines whether an iterator is at the end position.
..signature:bool atEnd(iterator [, container])
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.RootedIteratorConcept
..param.container:Container of $iterator$.
...remarks.text:If $iterator$ implements @Concept.RootedIteratorConcept@ then $container$ is optional.
....text:$container$ is also optional for iterators to @Adaption.char array.char arrays@.
....text:Otherwise, $container$ is required.
..returns:$true$ if $iterator$ points behind the last item of the container, otherwise $false$.
..see:Function.atBegin
..see:Function.end
..include:seqan/basic.h
*/

template <typename T, typename TContainer>
inline bool
atEnd(T & it,
      TContainer const & cont)
{
    SEQAN_CHECKPOINT;
    return it == end(cont, Standard());
}

template <typename T, typename TContainer>
inline bool
atEnd(T const & it,
      TContainer const & cont)
{
    SEQAN_CHECKPOINT;
    return it == end(cont, Standard());
}

template <typename T, typename TContainer>
inline bool
atEnd(T & it,
      TContainer & cont)
{
    SEQAN_CHECKPOINT;
    return it == end(cont, Standard());
}

template <typename T, typename TContainer>
inline bool
atEnd(T const & it,
      TContainer & cont)
{
    SEQAN_CHECKPOINT;
    return it == end(cont, Standard());
}

template <typename T>
inline bool
atEnd(T & it)
{
    SEQAN_CHECKPOINT;
    return atEnd(it, container(it));
}

template <typename T>
inline bool
atEnd(T const & it)
{
    SEQAN_CHECKPOINT;
    return atEnd(it, container(it));
}

// ---------------------------------------------------------------------------
// Function goBegin()
// ---------------------------------------------------------------------------

/**
.Function.goBegin
..class:Class.Iter
..concept:Concept.RootedIteratorConcept
..cat:Iteration
..summary:Iterates to the first position of a container.
..signature:goBegin(iterator [, container])
..param.iterator:Object that iterates through $container$.
...type:Class.Iter
...concept:Concept.RootedIteratorConcept
...text:$iterator$ is set to the position of the first item in $container$.
..param.container:Container of $iterator$.
...type:Concept.ContainerConcept
...remarks.text:If $iterator$ implements @Concept.RootedIteratorConcept@ then $container$ is optional,
otherwise $container$ is required.
..remarks:This function is equivalent to $iterator = begin(container)$.
..see:Function.begin
..see:Function.atBegin
..see:Function.goEnd
..include:seqan/basic.h
*/

template <typename TIterator, typename TContainer>
inline void
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
inline void
goBegin(TIterator & it)
{
    typename Parameter_<typename Container<TIterator>::Type>::Type tmpContainer = container(it);
    goBegin(it, tmpContainer);
}

// ---------------------------------------------------------------------------
// Function goEnd()
// ---------------------------------------------------------------------------

/**
.Function.goEnd
..class:Class.Iter
..concept:Concept.RootedIteratorConcept
..cat:Iteration
..summary:Iterates to the last position of a container.
..signature:goEnd(iterator [, container])
..param.iterator:Object that iterates through $container$.
...type:Class.Iter
...concept:Concept.RootedIteratorConcept
...text:$iterator$ is set to the position behin the last item in $container$.
..param.container:Container of $iterator$.
...type:Concept.ContainerConcept
...remarks.text:If $iterator$ implements @Concept.RootedIteratorConcept@ then $container$ is optional,
otherwise $container$ is required.
..remarks:This function is equivalent to $iterator = end(container)$.
..see:Function.end
..see:Function.atEnd
..see:Function.goBegin
..see:Function.goEnd
..include:seqan/basic.h
*/

template <typename TIterator, typename TContainer>
inline void
goEnd(TIterator & it,
      TContainer & container)
{
    SEQAN_CHECKPOINT;
    it = end(container);
}

template <typename TIterator, typename TContainer>
inline void
goEnd(TIterator & it,
      TContainer const & container)
{
    SEQAN_CHECKPOINT;
    it = end(container);
}

template <typename TIterator>
inline void
goEnd(TIterator & it)
{
    SEQAN_CHECKPOINT;
    goEnd(it, container(it));
}

// ---------------------------------------------------------------------------
// Function goNext()
// ---------------------------------------------------------------------------

/**
.Function.goNext
..concept:Concept.ForwardIteratorConcept
..cat:Iteration
..summary:Iterates to next position.
..signature:goNext(iterator)
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.ForwardIteratorConcept
...text:$iterator$ is set to the next position of an iteration through its container.
..remarks:This function is equivalent to $++iterator$.
..see:Function.goBegin
..see:Function.goEnd
..include:seqan/basic.h
*/

template <typename TIterator>
inline void
goNext(TIterator & it)
{
    SEQAN_CHECKPOINT;
    ++it;
}

// ---------------------------------------------------------------------------
// Function goFurther()
// ---------------------------------------------------------------------------

/**
.Function.goFurther
..concept:Concept.RandomAccessIteratorConcept
..cat:Iteration
..summary:Iterates some steps further.
..signature:goFurther(iterator, steps)
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.RandomAccessIteratorConcept
...text:$iterator$ is set $steps$ positions further in the iteration through the container.
..param.steps:Number of steps $iterator$ should be moved further.
...remarks:If $iterator$ supports bidirectional iteration, $steps$ could also be negativ.
..remarks:This function is equivalent to $iterator += steps$ for random access iterators.
..see:Function.goNext
..see:Function.goPrevious
..include:seqan/basic.h
*/

template <typename TIterator, typename TDiff>
inline void
goFurther(TIterator & it,
          TDiff steps)
{   // return distance type from arbitrary argument
    it += steps;
}

// ---------------------------------------------------------------------------
// Function goPrevious()
// ---------------------------------------------------------------------------

/**
.Function.goPrevious
..concept:Concept.BidirectionalIteratorConcept
..cat:Iteration
..summary:Iterates to pevious position.
..signature:goPrevious(iterator)
..param.iterator:An iterator.
...type:Class.Iter
...concept:Concept.BidirectionalIteratorConcept
...text:$iterator$ is set to the pevious position of an iteration through its container.
..remarks:This function is equivalent to $--iterator$.
..see:Function.goBegin
..see:Function.goEnd
..see:Function.goNext
..include:seqan/basic.h
*/
template <typename TIterator>
inline void
goPrevious(TIterator & it)
{
    SEQAN_CHECKPOINT;
    --it;
}

// ---------------------------------------------------------------------------
// Function difference()
// ---------------------------------------------------------------------------

/**
.Function.difference
..concept:Concept.RandomAccessIteratorConcept
..cat:Iteration
..summary:The difference between two iterators.
..signature:difference(begin, end)
..param.begin:Iterator to the first position of a range.
...type:Class.Iter
...concept:Concept.RandomAccessIteratorConcept
..param.end:Iterator behind the last position of a range.
...type:Class.Iter
...concept:Concept.RandomAccessIteratorConcept
..returns:Length of the range between $begin$ and $end$.
..remarks:This function is equivalent to $end - begin$.
...text:Usually, $begin$ and $end$ have the same type.
..see:Function.begin
..see:Function.end
..see:Function.length
..include:seqan/basic.h
*/

template <typename TIterator>
inline
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

/**
.Function.goNil
..concept:Concept.RandomAccessIteratorConcept
..cat:Iteration
..summary:Moves iterator to nil position.
..signature:goNil(iterator)
..param.iterator:The iterator that will be moved.
...type:Concept.RandomAccessIteratorConcept
..remarks:$iterator$ is set to an invalid position, e.g. $NULL$ for pointer types.
..see:Function.clear
..include:seqan/basic.h
*/

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

/**
.Function.atNil
..concept:Concept.RandomAccessIteratorConcept
..cat:Iteration
..summary:Tests whether iterator is at nil position.
..signature:bool atNil(iterator)
..param.iterator:An iterator.
...type:Concept.RandomAccessIteratorConcept
..returns:$true$ if $iterator$ points to an ivalid position, e.g. $iterator$ is a $NULL$ pointer.
$false$ otherwise.
..see:Function.goNil
..include:seqan/basic.h
*/

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

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ITERATOR_INTERFACE_H_
