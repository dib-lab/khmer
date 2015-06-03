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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Default implementations for sequences.

// TODO(holtgrew): There appears to be some overlap with string_base.h. Maybe it is a good idea to move everything related to strings to string_base.h and remove all default-container behaviour.

// TODO(holtgrew): These functions have (documentation wise) mostly gone into Container, Sequence Concepts and String class.  This is where they belong.

#ifndef SEQAN_SEQUENCE_SEQUENCE_INTERFACE_H_
#define SEQAN_SEQUENCE_SEQUENCE_INTERFACE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @defgroup OverflowStrategyTags Overflow Strategy Tags
 * @brief The strategy for resizing containers.
 *
 * Changing the capacity of a container can invalidate the iterators of this container.
 *
 * If no overflow tag is specified, most operations use the default overflow strategy given by @link
 * DefaultOverflowImplicit @endlink or @link DefaultOverflowExplicit @endlink, depending on the kind of operation.
 *
 * @see StringConcept#computeGenerousCapacity
 * @see DefaultOverflowImplicit
 * @see DefaultOverflowExplicit
 *
 * @tag OverflowStrategyTags#Limit
 * @headerfile <seqan/sequence.h>
 * @brief Limit the contents to current capacity.
 *
 * @signature typedef Tag<TagLimit_> Limit;
 *
 * All entries that exceed the capacity are lost.
 *
 * @tag OverflowStrategyTags#Generous
 * @headerfile <seqan/sequence.h>
 * @brief Expand if needed, get precautionary extra space.
 *
 * @signature typedef Tag<TagGenerous_> Generous;
 *
 * Whenever the capacity has to be increased, the new capacity is choosen somewhat large than actually needed.  This
 * strategy limits the number of capacity changes, so that resizing takes armotized constant time.  Use this strategy if
 * the total amount of storage is unkown at first.
 *
 * The new capacity is computed by @link StringConcept#computeGenerousCapacity @endlink. By default, it is
 * guaranteed not to exceed about three halfs of the space that is used to store the data.  The user can overload
 * @link StringConcept#computeGenerousCapacity @endlink in order to change this behavior.
 *
 * @tag OverflowStrategyTags#Exact
 * @headerfile <seqan/sequence.h>
 * @brief Expand as far as needed.
 *
 * @signature typedef Tag<TagExact_> Exact;
 *
 * The capacity is only changed if the current capacity is not large enough.  If the capacity can only be expanded up to
 * a certain ammount, it will be increased as far as possible  and the contents are limited to the new capacity.
 *
 * Note that the capacity will never be shrinked.  Use @link ContainerConcept#shrinkToFit @endlink to resize the
 * capacity down to the current length.
 *
 * @tag OverflowStrategyTags#Insist
 * @headerfile <seqan/sequence.h>
 * @brief No capacity check.
 *
 * @signature typedef Tag<TagInsist_> Insist;
 *
 * @note The user has to ensure that the container's capacity is large enough.
 */

struct TagInsist_;
typedef Tag<TagInsist_> Insist;
typedef Tag<TagInsist_> Tight;  // TODO(holtgrew): Necessary?

struct TagLimit_;
typedef Tag<TagLimit_> Limit;

struct TagGenerous_;
typedef Tag<TagGenerous_> Generous;

struct TagExact_;
typedef Tag<TagExact_> Exact;

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// --------------------------------------------------------------------------

/*!
 * @mfn DefaultOverflowImplicit
 * @headerfile <seqan/sequence.h>
 * @brief The default overflow strategy for implicit resize.
 *
 * @signature DefaultOverflowImplicit<T>::Type;
 *
 * @tparam T The type to get the default overflow tag for.
 *
 * @return Type The default overflow tag.  The default implementation returns <tt>Generous</tt>.
 *
 * This function is used for functions that cause an implicit change of a container's size, like e.g. through
 * @link AssignableConcept#assign @endlink, @link ContainerConcept#append @endlink and
 * @link StringConcept#replace @endlink.
 */

template <typename T>
struct DefaultOverflowImplicit
{
    typedef Generous Type;
};

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowExplicit
// --------------------------------------------------------------------------

/*!
 * @mfn DefaultOverflowExplicit
 * @headerfile <seqan/sequence.h>
 * @brief The default overflow strategy for explicit resize.
 *
 * @signature DefaultOverflowExplicit<T>::Type;
 *
 * @tparam T The type to determine overflow strategy.
 *
 * @return Type The resulting expantion tag for <tt>T</tt>.
 *
 * This function is used for functions that change a container's size explicitly, like e.g.
 * @link StringConcept#resize @endlink.
 */

template <typename T>
struct DefaultOverflowExplicit
{
    typedef Generous Type;
};

// --------------------------------------------------------------------------
// Metafunction IsSequence
// --------------------------------------------------------------------------

// TODO(holtgrew): Deprecate in favour of Is<StringConcept>?

/*!
 * @mfn IsSequence
 * @headerfile <seqan/sequence.h>
 * @brief Determines whether a type is a sequence.
 *
 * @signature IsSequence<T>::Type;
 * @signature IsSequence<T>::VALUE;
 *
 * @tparam T The type to query.
 *
 * @return Type  <tt>True</tt> if <tt>T</tt> is a sequence and <tt>False</tt> otherwise.
 * @return VALUE <tt>true</tt> if <tt>T</tt> is a sequence and <tt>false</tt> otherwise.
 *
 * For example, String and Segment as <tt>T</tt> return true.
 */

template <typename T>
struct IsSequence
{
    typedef False Type;
    enum { VALUE = false };
};

template <typename T>
struct IsSequence<T const> : IsSequence<T> {};

// --------------------------------------------------------------------------
// Metafunction AllowsFastRandomAccess
// --------------------------------------------------------------------------

/*!
 * @mfn AllowsFastRandomAccess
 * @headerfile <seqan/sequence.h>
 * @brief Determines whether a sequence efficiently supports random access.
 *
 * @signature AllowsFastRandomAccess<T>::Type;
 * @signature AllowsFastRandomAccess<T>::VALUE;
 *
 * @tparam T The type to query.
 *
 * @return Type  <tt>True</tt> if <tt>T</tt> allows for fast random access and <tt>False</tt> otherwise.
 * @return VALUE <tt>true</tt> if <tt>T</tt> allows for fast random access and <tt>false</tt> otherwise.
 *
 * For example, String and std::vector allow for fast random access, while std::list does not.
 */

template <typename T>
struct AllowsFastRandomAccess
{
    typedef True Type;
    enum { VALUE = true };
};

template <typename T>
struct AllowsFastRandomAccess<T const>
    : public AllowsFastRandomAccess<T> {};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function getObjectId()
// --------------------------------------------------------------------------

/*!
 * @fn ContainerConcept#getObjectId
 * @headerfile <seqan/sequence.h>
 * @brief A value that identifies the underlying sequence.
 *
 * @signature TVoidPtr getObjectId(cont);
 *
 * @param[in] cont The object for which to determine the id.
 *
 * @return TVoidPtr a <tt>void const *</tt> value identying the object.
 *
 * Two sequences should have the same id, if they share the same resource, e.g. the same memory buffer.
 *
 * The exact semantic of the returned id can vary for different classes.  Typically, the id of a string is a <tt>void
 * const *</tt> to the end of the string.
 *
 * @section Examples
 *
 * @code{.cpp}
 * String<char> str = "hallo seqan";
 * bool b1 = (getObjectId(str) == getObjectId(infix(str, 3, 7));   //true
 * bool b2 = (getObjectId(str) == getObjectId(String<char>(str))); //false
 * bool b3 = (getObjectId(str) == getObjectId(toCString(str)));
 * @endcode
 *
 * In this example, <tt>b1</tt> is <tt>true</tt., since the segment object returned by <tt>infix()</tt> is just a filter
 * and uses the buffer of it's host object str.
 *
 * <tt>String&lt;char&gt;(str)</tt> constructs a temporary copy of <tt>str</tt>, so these two strings have different id values.
 *
 * The result of the last comparison depends on the implementation of <tt>toCString</tt> and cannot be predicted at
 * compile time.
 */

template <typename T>
inline void const *
getObjectId(T const & me)
{
    SEQAN_CHECKPOINT;
    return end(me, Standard());
}

// --------------------------------------------------------------------------
// Function shareResources()
// --------------------------------------------------------------------------

/*!
 * @fn shareResources
 * @headerfile <seqan/sequence.h>
 * @brief Determines whether two sequences share the same resource.
 *
 * @signature bool shareResources(s1, s2);
 *
 * @param[in] s1 First sequence.
 * @param[in] s2 Second sequence.
 *
 * @return bool <tt>true</tt> if the two sequences share resources and <tt>false</tt> if not.
 */

template <typename T1, typename T2>
inline bool
shareResources(T1 const & obj1,
               T2 const & obj2)
{
    SEQAN_CHECKPOINT;
    return getObjectId(obj1) == getObjectId(obj2);
}

// --------------------------------------------------------------------------
// Function _beginDefault()
// --------------------------------------------------------------------------

//* ???Anti Default Sequences
// TODO(holtgrew): Evil -- each value is a container of length 1.
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, Standard>::Type
_beginDefault(T & me,
               Standard)
{
    SEQAN_CHECKPOINT;
    return & me;
}
// TODO(holtgrew): Evil -- each value is a container of length 1.
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Standard>::Type
_beginDefault(T const & me,
               Standard)
{
    SEQAN_CHECKPOINT;
    return & me;
}
//*/

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, Rooted>::Type
_beginDefault(T & me,
               Rooted)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<T, Rooted>::Type TIterator;
    return TIterator(me, begin(me, Standard()));
}
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Rooted>::Type
_beginDefault(T const & me,
               Rooted)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<T const, Rooted>::Type TIterator;
    return TIterator(me, begin(me, Standard()));
}

// --------------------------------------------------------------------------
// Function begin()
// --------------------------------------------------------------------------

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type
begin(T & me)
{
    SEQAN_CHECKPOINT;
    return begin(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type
begin(T const & me)
{
    SEQAN_CHECKPOINT;
    return begin(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

//folgende forward Deklaration wurde wegen Phaenomene bei VC++ 2003 hinzugenommen
//implemented in string_pointer.h
template <typename TValue>
SEQAN_HOST_DEVICE inline typename Iterator<TValue const *, Standard>::Type
begin(TValue const * me,
      Standard);

template <typename T, typename TSpec>
SEQAN_HOST_DEVICE inline typename Iterator<T, Tag<TSpec> const>::Type
begin(T & me,
      Tag<TSpec> const tag_)
{
    SEQAN_CHECKPOINT;
    return _beginDefault(me, tag_);
}
template <typename T, typename TSpec>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Tag<TSpec> const>::Type
begin(T const & me,
      Tag<TSpec> const tag_)
{
    SEQAN_CHECKPOINT;
    return _beginDefault(me, tag_);
}

/*
template <typename TValue>
inline typename Iterator<TValue *, Standard>::Type
begin(TValue * me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return me;
}

//folgende Version wurde wegen eines seltsamen Phaenomens bei VC++ hinzugenommen
template <typename TValue>
inline typename Iterator<TValue const *, Standard>::Type
begin(TValue const * me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return me;
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue *, Standard>::Type
begin(TValue * me,
      Tag<TSpec> const tag_)
//    Standard)
{
    SEQAN_CHECKPOINT;
    return me;
}

//folgende Version wurde wegen eines seltsamen Phaenomens bei VC++ hinzugenommen
template <typename TValue, typename TSpec>
inline typename Iterator<TValue const *, Standard>::Type
begin(TValue const * me,
      Tag<TSpec> const tag_)
//    Standard)
{
    SEQAN_CHECKPOINT;
    return me;
}
*/

// --------------------------------------------------------------------------
// Function beginPosition()
// --------------------------------------------------------------------------

template <typename T>
inline typename Position<T>::Type
beginPosition(T &)
{
    SEQAN_CHECKPOINT;
    return 0;
}

template <typename T>
inline typename Position<T>::Type
beginPosition(T const &)
{
    SEQAN_CHECKPOINT;
    return 0;
}

// --------------------------------------------------------------------------
// Function _endDefault()
// --------------------------------------------------------------------------

//* ???Anti Default Sequences
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, Standard>::Type
_endDefault(T & me,
             Standard)
{
    SEQAN_CHECKPOINT;
    return (& me) + 1;
}
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Standard>::Type
_endDefault(T const & me,
             Standard)
{
    SEQAN_CHECKPOINT;
    return (& me) + 1;
}
//*/

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, Rooted>::Type
_endDefault(T & me,
             Rooted)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<T, Rooted>::Type TIterator;
    return TIterator(me, end(me, Standard()));
}
template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Rooted>::Type
_endDefault(T const & me,
             Rooted)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<T const, Rooted>::Type TIterator;
    return TIterator(me, end(me, Standard()));
}

// --------------------------------------------------------------------------
// Function end()
// --------------------------------------------------------------------------

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type
end(T & me)
{
    SEQAN_CHECKPOINT;
    return end(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

template <typename T>
SEQAN_HOST_DEVICE inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type
end(T const & me)
{
    SEQAN_CHECKPOINT;
    return end(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

template <typename T, typename TSpec>
SEQAN_HOST_DEVICE inline typename Iterator<T, Tag<TSpec> const>::Type
end(T & me,
    Tag<TSpec> const tag_)
{
    SEQAN_CHECKPOINT;
    return _endDefault(me, tag_);
}

template <typename T, typename TSpec>
SEQAN_HOST_DEVICE inline typename Iterator<T const, Tag<TSpec> const>::Type
end(T const & me,
    Tag<TSpec> const tag_)
{
    SEQAN_CHECKPOINT;
    return _endDefault(me, tag_);
}

// --------------------------------------------------------------------------
// Function endPosition()
// --------------------------------------------------------------------------

template <typename T>
inline typename Position<T>::Type
endPosition(T & me)
{
    SEQAN_CHECKPOINT;
    return length(me);
}

template <typename T>
inline typename Position<T>::Type
endPosition(T const & me)
{
    SEQAN_CHECKPOINT;
    return length(me);
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

//* ???Anti Default Sequences
template <typename T, typename TPos>
SEQAN_HOST_DEVICE inline typename Reference<T>::Type
value(T & me,
      TPos /*pos*/)
{
    SEQAN_CHECKPOINT;
    return me;
}

template <typename T, typename TPos>
SEQAN_HOST_DEVICE inline typename Reference<T const>::Type
value(T const & me,
      TPos /*pos*/)
{
    SEQAN_CHECKPOINT;
    return me;
}
//*/

// --------------------------------------------------------------------------
// Function getValue()
// --------------------------------------------------------------------------

template <typename T, typename TPos>
inline typename GetValue<T>::Type
getValue(T & me,
         TPos pos)
{
    SEQAN_CHECKPOINT;
    return (typename GetValue<T>::Type) value(me, pos);
}

template <typename T, typename TPos>
inline typename GetValue<T const>::Type
getValue(T const & me,
         TPos pos)
{
    SEQAN_CHECKPOINT;
    return value(me, pos);
}

// --------------------------------------------------------------------------
// Function front()
// --------------------------------------------------------------------------

template <typename T>
inline typename Reference<T>::Type
front(T & me)
{
    SEQAN_CHECKPOINT;
    return value(me, 0);
}
template <typename T>
inline typename Reference<T const>::Type
front(T const & me)
{
    SEQAN_CHECKPOINT;
    return value(me, 0);
}

// --------------------------------------------------------------------------
// Function back()
// --------------------------------------------------------------------------

template <typename T>
SEQAN_HOST_DEVICE inline typename Reference<T const>::Type
back(T const & me)
{
    SEQAN_CHECKPOINT;
    return value(me, length(me) - 1);
}

template <typename T>
SEQAN_HOST_DEVICE inline typename Reference<T>::Type
back(T & me)
{
    SEQAN_CHECKPOINT;
    return value(me, length(me) - 1);
}

template <typename T>
SEQAN_HOST_DEVICE inline typename Reference<T const>::Type
backPrev(T const & me)
{
    return value(me, length(me) - 2);
}

template <typename T>
SEQAN_HOST_DEVICE inline typename Reference<T>::Type
backPrev(T & me)
{
    return value(me, length(me) - 2);
}

// --------------------------------------------------------------------------
// Function iter()
// --------------------------------------------------------------------------

template <typename T, typename TPos>
inline typename Iterator<T, typename DefaultGetIteratorSpec<T>::Type>::Type
iter(T & me,
     TPos pos)
{
    SEQAN_CHECKPOINT;
    return iter(me, pos, typename DefaultGetIteratorSpec<T>::Type());
}

template <typename T, typename TPos>
inline typename Iterator<T const, typename DefaultGetIteratorSpec<T>::Type>::Type
iter(T const & me,
     TPos pos)
{
    SEQAN_CHECKPOINT;
    return iter(me, pos, typename DefaultGetIteratorSpec<T>::Type());
}

template <typename T, typename TPos, typename TTag>
inline typename Iterator<T, Tag<TTag> const>::Type
iter(T & me,
     TPos pos,
     Tag<TTag> const tag_)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LEQ_MSG(pos, static_cast<TPos>(length(me)), "Trying to get an iterator behind a container through iter().");
    return begin(me, tag_) + pos;
}

template <typename T, typename TPos, typename TTag>
inline typename Iterator<T const, Tag<TTag> const>::Type
iter(T const & me,
     TPos pos,
     Tag<TTag> const tag_)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LEQ_MSG(pos, static_cast<TPos>(length(me)), "Trying to get an iterator behind a container through iter().");
    return begin(me, tag_) + pos;
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

template <typename T, typename TValue, typename TPos>
inline void
assignValue(T & me,
            TPos pos,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    assign(value(me, pos), _value);
}

// --------------------------------------------------------------------------
// Function moveValue()
// --------------------------------------------------------------------------

/*!
 * @fn ContainerConcept#moveValue
 * @headerfile <seqan/sequence.h>
 * @brief Move a value into a container at a given position.
 *
 * @signature void moveValue(container, pos, value);
 *
 * @param[in,out] container The container to manipulate.
 * @param[in]     pos       The position of the item in the container to manipulate.
 * @param[in,out] value     The value to move to <tt>container[pos]</tt>.
 */

template <typename T, typename TValue, typename TPos>
inline void
moveValue(T & me,
          TPos pos,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    move(value(me, pos), _value);
}

template <typename T, typename TValue, typename TPos>
inline void
moveValue(T const & me,
          TPos pos,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    move(value(me, pos), _value);
}

// --------------------------------------------------------------------------
// Function length()
// --------------------------------------------------------------------------

//* ???Anti Default Sequences
template <typename T>
inline typename Size<T>::Type
length(T const & /*me*/)
{
    SEQAN_CHECKPOINT;
    return 1;
}
//*/

// --------------------------------------------------------------------------
// Function capacity()
// --------------------------------------------------------------------------

/*!
 * @fn StringConcept#capacity
 * @headerfile <seqan/sequence.h>
 * @brief Returns the capacity of a sequence.
 *
 * @signature TSize capacity(seq);
 *
 * @param[in] seq The sequence to query for its capacity.
 *
 * @return TSize Returns the capacity of the sequence.  <tt>TSize</tt> is the result of
 *               <tt>Size&lt;TSequence&gt;::type</tt> where <tt>TSequence</tt> is the type of <tt>seq</tt>.
 *
 * The size of a sequence can never exceed its capacity but some container support resizing of the capacity.  Some
 * functions do that implicitely if they are called with a suitable @link OverflowStrategyTags tag @endlink.  The
 * function reserve can be used to change the capacity explicitely.
 */

template <typename T>
inline SEQAN_HOST_DEVICE typename Size<T const>::Type
capacity(T const & me)
{
    SEQAN_CHECKPOINT;
    return length(me);
}

// --------------------------------------------------------------------------
// Function empty()
// --------------------------------------------------------------------------

template <typename T>
SEQAN_HOST_DEVICE inline bool
empty(T const & me)
{
    SEQAN_CHECKPOINT;
    return (length(me) == 0);
}

// --------------------------------------------------------------------------
// Function _computeSizeForCapacity()
// --------------------------------------------------------------------------

// note: for value types of size 1 or 2,
// an extra position for the termination character is allocated.
// This speeds up a conversion to a c style string (see Spec.CStyle String)
// note that this extra position is necessary not only for char and wchar_t,
// but also for all other value types of size 1 and 2 to make the application
// of the funciton move for in-place alphabet conversion.


template <typename T, typename TSize>
inline TSize
_computeSizeForCapacity(T const & /*me*/,
                      TSize capacity)
{
    SEQAN_CHECKPOINT;
    typedef typename Value<T>::Type TValue;
    if (sizeof(TValue) <= 2) return capacity + 1;
    else return capacity;
}

// --------------------------------------------------------------------------
// Function computeGenerousCapacity()
// --------------------------------------------------------------------------

/*!
 * @fn StringConcept#computeGenerousCapacity
 * @headerfile <seqan/sequence.h>
 * @brief Capacity for generous expansion.
 *
 * @signature TSize computeGenerousCapacity(seq, capacity);
 *
 * @param[in,out] seq       The sequence to compute the generous capacity for.
 * @param[in]     capacity  The minimal required capacity.
 *
 * @return TSize A value larger than <tt>capacity</tt> that should be used as the new capacity for <tt>container</tt>
 *               when it is expanded using the <tt>Generous</tt> overflow strategy.
 */

template <typename T, typename TSize>
inline TSize
computeGenerousCapacity(T const & /*me*/,
                         TSize capacity)
{
    SEQAN_CHECKPOINT;
    if (capacity < 32) return 32;       // returned value is implicitly >= capacity + 1
    return capacity + (capacity >> 1);
}

// --------------------------------------------------------------------------
// Function _storageUpdated()
// --------------------------------------------------------------------------

/*
template <typename T>
inline void
_storageUpdated(T & me,
                void const *)
{
}

template <typename T>
inline void
_storageUpdated(T & me)
{
    _storageUpdated_(me, (T *) 0);
}

template <typename T>
inline void
_storageUpdated(T const & me)
{
    _storageUpdated_(me, (T const *) 0);
}
*/

// --------------------------------------------------------------------------
// Function assign()
// --------------------------------------------------------------------------

template<typename TTarget, typename TSource>
inline void
assign(TTarget & target,
       TSource & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
assign(TTarget const & target,
       TSource & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TSource>
inline void
assign(TTarget & target,
       TSource const & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
assign(TTarget const & target,
       TSource const & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

// --------------------------------------------------------------------------
// Function append()
// --------------------------------------------------------------------------

/*!
 * @fn ContainerConcept#append
 * @headerfile <seqan/sequence.h>
 * @brief Concatenate a container to another.
 *
 * @signature void append(target, source);
 *
 * @param[in,out] target The @link ContainerConcept container @endlink to append <tt>source</tt> to.
 * @param[in]     source This @link ContainerConcept container @endlink will be appended to <tt>source</tt>.
 */

template<typename TTarget, typename TSource>
inline void
append(TTarget & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    append(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget const & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    append(target, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    append(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget const & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    append(target, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget & target,
       TSource & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget const & target,
       TSource & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget & target,
       TSource const & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TSource>
inline void
append(TTarget const & target,
       TSource const & source,
       typename Size<TTarget>::Type limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

// --------------------------------------------------------------------------
// Function appendValue()
// --------------------------------------------------------------------------

/*!
 * @fn ContainerConcept#appendValue
 * @headerfile <seqan/sequence.h>
 * @brief Append a value to a container.
 *
 * @signature void appendValue(target, val[, tag]);
 *
 * @param[in,out] target The container to append <tt>val</tt> to.
 * @param[in]     val    The value to append to <tt>target</tt>.
 * @param[in]     tag    The resize tag to use.  Defaults to What DefaultOverflowImplicit returns for the type of
 *                       <tt>target</tt>.
 */

template <typename T, typename TValue>
inline void
appendValue(T SEQAN_FORWARD_ARG me,
            TValue SEQAN_FORWARD_CARG _value)
{
    appendValue(SEQAN_FORWARD(T, me), SEQAN_FORWARD(TValue, _value), typename DefaultOverflowImplicit<T>::Type());
}

#ifndef SEQAN_CXX11_STANDARD

template <typename T, typename TValue>
inline void
appendValue(T const & me,
            TValue const & _value)
{
    appendValue(me, _value, typename DefaultOverflowImplicit<T const>::Type());
}

#endif  // #ifndef SEQAN_CXX11_STANDARD

// --------------------------------------------------------------------------
// Function insert()
// --------------------------------------------------------------------------

/*!
 * @fn StringConcept#insert
 * @headerfile <seqan/sequence.h>
 * @brief Inserts a sequence into another sequence.
 *
 * @signature void insert(seq, pos, src[, tag]);
 *
 * @param[in,out] seq The sequence to insert element sinto.
 * @param[in]     pos The position to start inserting at.
 * @param[in]     src The sequence to insert at pos.
 * @param[in]     tag The resize tag, defaults to what <tt>OverflowStrategyImplicit</tt> returns.
 */

template <typename T, typename TPosition, typename TSeq, typename TExpand>
inline void
insert(T & me,
       TPosition pos,
       TSeq const & insertSeq,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    replace(me, pos, pos, insertSeq, Tag<TExpand>());
}

template <typename T, typename TPosition, typename TSeq, typename TExpand>
inline void
insert(T const & me,
       TPosition pos,
       TSeq const & insertSeq,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    replace(me, pos, pos, insertSeq, Tag<TExpand>());
}

template <typename T, typename TPosition, typename TSeq>
inline void
insert(T & me,
       TPosition pos,
       TSeq const & insertSeq)
{
    SEQAN_CHECKPOINT;
    replace(me, pos, pos, insertSeq, typename DefaultOverflowImplicit<T>::Type());
}

template <typename T, typename TPosition, typename TSeq>
inline void
insert(T const & me,
       TPosition pos,
       TSeq const & insertSeq)
{
    SEQAN_CHECKPOINT;
    replace(me, pos, pos, insertSeq, typename DefaultOverflowImplicit<T const>::Type());
}

// --------------------------------------------------------------------------
// Function insertValue()
// --------------------------------------------------------------------------

/*!
 * @fn StringConcept#insertValue
 * @headerfile <seqan/sequence.h>
 * @brief Inserts an element into a sequence.
 *
 * @signature void insertValue(seq, pos, val[, tag]);
 *
 * @param[in,out] seq  The @link StringConcept sequence @endlink to insert element into.
 * @param[in]     pos  The position to insert at.
 * @param[in]     val  The value to insert at <tt>pos</tt> into <tt>seq<tt/>.
 * @param[in]     tag  The resize tag, defaults to what <tt>OverflowStrategyImplicit</tt> returns.
 */

template <typename T, typename TPosition, typename TValue>
inline void
insertValue(T & me,
            TPosition pos,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    insertValue(me, pos, _value, typename DefaultOverflowImplicit<T>::Type());
}

template <typename T, typename TPosition, typename TValue>
inline void
insertValue(T const & me,
            TPosition pos,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    insertValue(me, pos, _value, typename DefaultOverflowImplicit<T const>::Type());
}

// --------------------------------------------------------------------------
// Function replace()
// --------------------------------------------------------------------------

/*!
 * @fn StringConcept#replace
 * @headerfile <seqan/sequence.h>
 * @brief Replaces a part of a sequence with another sequence.
 *
 * @signature void replace(target, posBegin, posEnd, source[, limit][, resizeTag]);
 *
 * @param[in,out] target    The sequence to modify.
 * @param[in]     posBegin  The begin position of the range to replace.
 * @param[in]     posEnd    The end position of the range to replace.
 * @param[in]     source    The source sequence to replace <tt>[posBegin, posEnd)</tt> with.
 * @param[in]     limit     Largest size of <tt>target</tt> after the operation.
 * @param[in]     resizeTag Specify the resizing behaviour.  Defaults to what <tt>DefaultOverflowImplicit</tt>
 *                          returns.
 */

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource & source)
{
    replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget const & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource & source)
{
    replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source)
{
    replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget const & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source)
{
    replace(target, pos_begin, pos_end, source, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource & source,
        typename Size<TTarget>::Type limit)
{
    replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget const & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource & source,
        typename Size<TTarget>::Type limit)
{
    replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source,
        typename Size<TTarget>::Type limit)
{
    replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTarget, typename TPositionBegin, typename TPositionEnd, typename TSource>
inline void
replace(TTarget const & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source,
        typename Size<TTarget>::Type limit)
{
    replace(target, pos_begin, pos_end, source, limit, typename DefaultOverflowImplicit<TTarget const>::Type());
}

// --------------------------------------------------------------------------
// Function _capacityReturned()
// --------------------------------------------------------------------------

// TODO(holtgrew): Is this internal or a helper?

template <typename T, typename TSize, typename TExpand>
inline typename Size<T>::Type
_capacityReturned(T & me,
                  TSize,
                  Tag<TExpand>)
{
    return capacity(me);
}

template <typename T, typename TSize>
inline typename Size<T>::Type
_capacityReturned(T &,
                  TSize new_capacity,
                  Insist const & )
{
    return new_capacity;
}

// --------------------------------------------------------------------------
// Function reserve()
// --------------------------------------------------------------------------

/*!
 * @fn String#reserve
 * @brief Increases the capacity.
 *
 * @signature TSize reserve(str, new_capacity[, tag]);
 *
 * @param[in,out] str         The String to reserve space in.
 * @param[in]     newCapacity The new capacity <tt>str</tt> will get.
 * @param[in]     tag         Specifies the strategy that is applied for changing the capacity.
 *
 * @return TSize The amount of the requested capacity that was available.  That is the function returns the minimum of
 *               <tt>newCapacity</tt> and <tt>capacity(me)</tt>.
 *
 * This function allows to increase the capacity but not the length of a container.
 *
 * Use @link StringConcept#resize @endlink if you want to change the size of a container.
 *
 * @section Remarks
 *
 * At the end of the operation, <tt>capacity(me)</tt> can be larger than <tt>new_capacity</tt>.  If
 * <tt>new_capacity</tt> is smaller than <tt>capacity(me)</tt> at the beginning of the operation, the operation need not
 * to change the capacity at all.
 *
 * This operation does not changes the content of <tt>object</tt>.
 *
 * This operation may invalidate iterators of <tt>object</tt>.
 */

template <typename T, typename TSize, typename TExpand>
inline typename Size<T>::Type
reserve(T & me,
        TSize const & new_capacity,
        Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    return _capacityReturned(me, new_capacity, tag);
}

template <typename T, typename TSize>
inline typename Size<T>::Type
reserve(T & me,
        TSize const & new_capacity)
{
    SEQAN_CHECKPOINT;
    return reserve(me, new_capacity, typename DefaultOverflowExplicit<T>::Type());
}

// --------------------------------------------------------------------------
// Function resize()
// --------------------------------------------------------------------------

template <typename T, typename TSize>
inline typename Size<T>::Type
resize(T & me,
       TSize new_length)
{
    SEQAN_CHECKPOINT;
    return resize(me, new_length, typename DefaultOverflowExplicit<T>::Type());
}

template <typename T, typename TSize, typename TValue>
inline typename Size<T>::Type
resize(T & me,
       TSize new_length,
       TValue const & val)
{
    SEQAN_CHECKPOINT;
    return resize(me, new_length, val, typename DefaultOverflowExplicit<T>::Type());
}

// --------------------------------------------------------------------------
// Function resizeSpace()
// --------------------------------------------------------------------------

// TODO(holtgrew): Deprecated!

/*!
 * @fn String#resizeSpace
 * @headerfile <seqan/sequence.h>
 * @brief Makes free space in container
 *
 * @signature TSize resizeSpace(str, size, posBegin, posEnd [, limit][, resizeTag]);
 *
 * @param[in,out] str       The String to modify.
 * @param[in]     size      Number of characters that should be freed.
 * @param[in]     posEnd    Position behind the last item in <tt>object</tt> that is to be destroyed.  If
 *                          <tt>posEnd == posBegin</tt>, no item in <tt>object</tt> will be destroyed.
 * @param[in]     posBegin  Position of the first item in <tt>object</tt> that is to be destroyed.
 * @param[in]     limit     Maximal length <tt>object</tt> can get after this operation. (optional)
 * @param[in]     resizeTag Strategy that is applied if <tt>object</tt> has not enough capacity to store the
 *                          complete content. (optional)
 *
 * @return TSize The number of free characters.Depeding on resizeTag, this could be <tt>size</tt> or less than
 *               <tt>size</tt> if <tt>object</tt> has not enough <tt>capacity</tt>.
 */

template<typename T, typename TSize, typename TBeginPosition, typename TEndPosition>
inline TSize
resizeSpace(T & me,
            TSize size,
            TBeginPosition pos_begin,
            TEndPosition pos_end)
{
    SEQAN_CHECKPOINT;
    return resizeSpace(me, size, pos_begin, pos_end, typename DefaultOverflowExplicit<T>::Type());
}

template<typename T, typename TSize, typename TBeginPosition, typename TEndPosition, typename TLimit>
inline TSize
resizeSpace(T & me,
            TSize size,
            TBeginPosition pos_begin,
            TEndPosition pos_end,
            TLimit limit)
{
    SEQAN_CHECKPOINT;
    return resizeSpace(me, size, pos_begin, pos_end, limit, typename DefaultOverflowExplicit<T>::Type());
}

// --------------------------------------------------------------------------
// Function erase()
// --------------------------------------------------------------------------

template<typename T, typename TBeginPosition, typename TEndPosition>
inline void
erase(T & me,
      TBeginPosition pos,
      TEndPosition pos_end)
{
    SEQAN_CHECKPOINT;
    resizeSpace(me, 0, pos, pos_end);
}

template<typename T, typename TPosition>
inline void
erase(T & me,
      TPosition pos)
{
    SEQAN_CHECKPOINT;
    resizeSpace(me, 0, pos, pos + 1);
}

// For segments, we also have to define the version for const-containers.

template<typename T, typename TBeginPosition, typename TEndPosition>
inline void
erase(T const & me,
      TBeginPosition pos,
      TEndPosition pos_end)
{
    SEQAN_CHECKPOINT;
    resizeSpace(me, 0, pos, pos_end);
}

template<typename T, typename TPosition>
inline void
erase(T const & me,
      TPosition pos)
{
    SEQAN_CHECKPOINT;
    resizeSpace(me, 0, pos, pos + 1);
}

// --------------------------------------------------------------------------
// Function eraseBack()
// --------------------------------------------------------------------------

template <typename T>
inline void eraseBack(T & me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GT_MSG(length(me), 0u, "String must have more than 0 characters in eraseBack()!");
    resize(me, length(me) - 1);
}

// --------------------------------------------------------------------------
// Function shrinkToFit()
// --------------------------------------------------------------------------

/*!
 * @fn ContainerConcept#shrinkToFit
 * @headerfile <seqan/sequence.h>
 * @brief Resizes container to minimum capacity.
 *
 * @signature void shrinkToFit(cont);
 *
 * @param[in] cont The container to shrink.
 */

template <typename T>
inline void
shrinkToFit(T & me)
{
    SEQAN_CHECKPOINT;

//  following line has no effect as in SeqAn it is not yet possible
//  to reduce the memory consumption of a string with resize/reserve
//
//  reserve(me, length(me), Exact());

    T tmp;
    assign(tmp, me, Exact());
    swap(me, tmp);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_SEQUENCE_INTERFACE_H_
