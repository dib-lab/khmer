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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Concept definitions for containers.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_CONTAINER_CONCEPT_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_CONTAINER_CONCEPT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T> struct Infix;
template <typename T> struct Prefix;
template <typename T> struct Suffix;

// ============================================================================
// Concepts
// ============================================================================

// Forwards.
struct Standard_;
typedef Tag<Standard_> const Standard;
template <typename TContainer, typename TSpec> struct Iterator;

/*!
 * @concept ContainerConcept
 * @extends AssignableConcept
 * @extends DestructibleConcept
 * @brief A container is an object that stores other objects (<i>elements</i>).
 * @headerfile <seqan/basic.h>
 *
 * @signature ContainerConcept<T>
 *
 * Containers store multiple entries of the same type (the <i>element type</i>) and provide means to access these
 * items.  More specific, each container has an iterator type that is used for accessing its elements.
 *
 * There is no guarantee for the elements to be in a particular order (the order can vary between two iterations) and
 * no guarantee for the time complexity of element access.  Furthermore, there is no guarantee that there can be more
 * than one iterator in the container.  Modification of a container through an iterator invalidates all other iterators.
 *
 * Refinements of the Container concept or specific implementations can provide these guarantees, however.
 *
 * A container owns its elements and the elements are destructed when their owning container is destructed.  The
 * elements must fulfill the concepts @link AssignableConcept @endlink and @link DestructibleConcept @endlink.
 */

/*!
 * @mfn ContainerConcept#Value
 * @brief Returns the value type of the container.
 *
 * @signature Value<TContainer>::Type
 *
 * @tparam TContainer The Container to query.
 *
 * @return Type The element type of the container.
 *
 * The value type is the type that can be used for storing copies of the elements in the container.
 *
 * @section Valid Expressions
 *
 * The variable v has the value type of the container TContainer whereas it is an iterator into the container. Thus,
 * copies of values from TContainer (*it) ca be stored in v.
 *
 * @code{.cpp}
 * Value<TContainer>::Type v = *it;
 * @endcode
 */

/*!
 * @mfn ContainerConcept#GetValue
 * @brief Returns the get-value type of the container.
 *
 * @signature GetValue<TContainer>::Type
 *
 * @tparam TContainer The Container to query.
 *
 * @return Type The get-value type of the container.
 *
 * The get-value type of the container is a type for efficient read-only access to the elements in the container.
 * For small types (e.g. <tt>int</tt>), this can be a copy (thus <tt>int</tt>), for larger types, this can be
 * a const reference to the value type.
 *
 * @section Valid Expressions
 *
 * The variable v has the get-value type of the container TContainer whereas it is an iterator into the container.
 * Thus, we can store a get-value in v.
 *
 * @code{.cpp}
 * GetValue<TContainer>::Type v = *it;
 * @endcode
 */

/*!
 * @mfn ContainerConcept#Reference
 * @brief Returns the reference type of the container.
 *
 * @signature Reference<TContainer>::Type
 *
 * @tparam TContainer The Container to query.
 *
 * @return Type The reference type of the container.
 *
 * Different from STL containers, the const-ness of <tt>TContainer</tt> decides whether the returned type is a
 * const reference or a reference for modifying elements.
 *
 * Note that the reference type is not guaranteed to be <tt>TValue &amp;</tt> if the value type of the container
 * is <tt>TValue</tt>.  The reference can be implemented as a proxy, similar to <tt>std::vector&lt;bool&gt;</tt>.
 *
 * @section Valid Expressions
 *
 * The variable r has the reference type of the container TContainer whereas it is an iterator into the container.
 * Thus, we can store a reference to a value in the container in r.  Then, we can assign the value of v, a value
 * of the container.
 *
 * @code{.cpp}
 * Reference<TContainer>::Type r = *it;  // reference into container
 * r = v;  // updates value in container, thus also *it
 * @endcode
 */

/*!
 * @mfn ContainerConcept#Iterator
 * @brief Returns the iterator type of the container.
 *
 * @signature Iterator<TContainer[, TSpec]>::Type
 *
 * @tparam TContainer The Container type to query.
 * @tparam TSpec      Optionally, a tag for selecting the kind of iterator.  If not given, then
 *                    @link ContainerConcept#DefaultIteratorSpec @endlink of TContainer is used.  When given, one of
 *                    <tt>Standard</tt> and <tt>Rooted</tt>.
 *
 * @return Type       The iterator type.
 *
 * Different from the STL the <tt>const</tt> attribute of <tt>TContainer</tt> determines whether the resturned
 * iterator is a const iterator or a non-const iterator.
 *
 * @see ContainerIteratorTags
 */

/*!
 * @mfn ContainerConcept#Size
 * @brief Returns the size type of a container.
 *
 * @signature Size<TContainer>::Type
 *
 * @tparam TContainer The Container type to query.
 * @return Type       The type to use for storing container sizes.
 */

// TODO(holtgrew): This should become RandomAccessContainer#Position.
/*!
 * @mfn ContainerConcept#Position
 * @brief Returns the position type of a container.
 *
 * @signature Positin<TContainer>::Type
 *
 * @tparam TContainer The Container type to query.
 * @return Type       The type to use for storing container positions.
 */

/*!
 * @mfn ContainerConcept#Difference
 * @brief Returns the type for distances between two iterators.
 *
 * @signature Size<TContainer>::Type
 *
 * @tparam TContainer The Container type to query.
 * @return Type       The type to use for storing iterator distances sizes.
 *
 * This must be the same type as the distance type of the containers iterators.
 */

/*!
 * @fn ContainerConcept#begin
 * @brief Returns an iterator to the beginning of the container.
 *
 * @signature TIterator begin(c[, tag]);
 *
 * @param[in] c   The container to get the begin iterator for (type <tt>TContainer</tt>).
 * @param[in] tag An optional tag for selecting the type of the iterator.  One of <tt>Standard</tt> and <tt>Rooted</tt>.
 *                When left out, @link ContainerConcept#DefaultGetIteratorSpec @endlink of <tt>TContainer</tt> is used.
 *
 * @return TIterator Iterator to the beginning of the container, the type is selected by @link ContainerConcept#Iterator @endlink with
 *                   the given (or default) tag.
 *
 * When empty, <tt>begin(c) == end(c)</tt>.
 */

/*!
 * @fn ContainerConcept#end
 * @brief Returns an iterator to the end of the container.
 *
 * @signature TIterator end(c[, tag]);
 *
 * @param[in] c   The container to get the end iterator for (type <tt>TContainer</tt>).
 * @param[in] tag An optional tag for selecting the type of the iterator.  One of <tt>Standard</tt> and <tt>Rooted</tt>.
 *                When left out, @link ContainerConcept#DefaultGetIteratorSpec @endlink of <tt>TContainer</tt> is used.
 *
 * @return TIterator Iterator to the end of the container, the type is selected by @link ContainerConcept#Iterator @endlink with
 *                   the given (or default) tag.
 *
 * When empty, <tt>begin(c) == end(c)</tt>.
 */

/*!
 * @fn ContainerConcept#length
 * @brief Returns the size of the container.
 *
 * @signature TSize length(c);
 *
 * @param[in] c The container to query for its size.
 *
 * @return TSize The number of elements in the container.
 */

/*!
 * @fn ContainerConcept#empty
 * @brief Returns whether the container is empty.
 *
 * @signature bool empty(c);
 *
 * @param[in] c The container to query.
 *
 * @return bool Whether or not the container is empty.
 */

/*!
 * @fn ContainerConcept#swap
 * @brief Swap the contents of two containers.
 *
 * @signature void swap(c1, c2);
 *
 * @param[in,out] c1 The first container.
 * @param[in,out] c2 The second container.
 *
 * Swaps the contents of <tt>c1</tt> and <tt>c2</tt>.  The <tt>swap</tt> function must be defined in the same
 * namespace as the container for Koenig lookup to work.  In the heart of sorting algorithms, for example,
 * the swap function is properly used as follows.  This way, both the generic <tt>std::swap</tt> and the specialized
 * <tt>swap</tt> function of the container are available:
 *
 * @code{.cpp}
 * TContainer c1, c2; // ...
 * using std::swap;
 * swap(c1, c2);
 * @endcode
 */

// mutable container concept
template <typename TContainer>
struct ContainerConcept :
    Assignable<TContainer>,
    DefaultConstructible<TContainer>,
    CopyConstructible<TContainer>
{
    typedef typename Value<TContainer>::Type                TValue;
    typedef typename GetValue<TContainer>::Type             TGetValue;
    typedef typename Reference<TContainer>::Type            TReference;
    typedef typename Size<TContainer>::Type                 TSize;
    typedef typename Position<TContainer>::Type             TPosition;
    typedef typename Difference<TContainer>::Type           TDifference;
    typedef typename Iterator<TContainer, Standard>::Type   TIterator;

    TContainer  c, c2;
    TValue      val;
    TSize       size;
    TPosition   pos;
    TDifference diff;
    TIterator   iter;

    SEQAN_CONCEPT_ASSERT((AlphabetConcept<TValue>));
    SEQAN_CONCEPT_ASSERT((SignedIntegerConcept<TDifference>));
    SEQAN_CONCEPT_ASSERT((UnsignedIntegerConcept<TSize>));

    SEQAN_CONCEPT_USAGE(ContainerConcept)
    {
        // test return of const values
        sameType(getValue(c, 0), val);

        // TODO(holtgrew): Index based access requires random-access and *linear* container.
        // test whether returned references/proxies
        // can be assigned to val and vice versa
        val = value(c, 0);
        value(c, 0) = val;
        moveValue(c, pos, val);

        // sameType(value(str, 0), val); would not work
        // for Strings returning proxies, e.g. String<.., Packed>

        // TODO(holtgrew): Too strong assumption about random access iterators IMO.
        // test iterators
        sameType(iter, begin(c, Standard()));
        sameType(iter, end(c, Standard()));
        sameType(diff, end(c, Standard()) - begin(c, Standard()));

        // length and empty
        sameType(size, length(c));
        sameType(true, empty(c));

        // TODO: infix/suffix/prefix

        // swap containers
//        swap(c, c2);          // swap is not yet supported by every string
    }
};

SEQAN_CONCEPT_REFINE(StringConcept, (TString), (ContainerConcept)(PropertyMapConcept))
{
    typedef typename Value<TString>::Type                 TValue;
    typedef typename Size<TString>::Type                  TSize;
    typedef typename Position<TString>::Type              TPos;
    typedef typename Difference<TString>::Type            TDifference;
    typedef typename Iterator<TString, Standard>::Type    TIterator;

    TValue      val;
    TSize       size;
    TPos        pos;

    TString     str, str2;

    SEQAN_CONCEPT_USAGE(StringConcept)
    {
        pos = 0u;

        // clear
        clear(str);

        // append
        append(str, str2);
        appendValue(str, val);

        // capacity
        sameType(size, capacity(str));
    }
};

// --------------------------------------------------------------------------
// Metafunction IsContiguous
// --------------------------------------------------------------------------

/*!
 * @mfn IsContiguous
 * @headerfile <seqan/sequence.h>
 * @brief Determines whether a container stores its elements contiguously in memory.
 *
 * @signature IsContiguous<T>::Type;
 * @signature IsContiguous<T>::VALUE;
 *
 * @tparam T The type that is tested for being a string.
 *
 * @return Type  Either <tt>True</tt> or <tt>False</tt>, depending on whether <tt>T</tt> is stored contiguously.
 * @return VALUE Either <tt>true</tt> or <tt>false</tt>, depending on whether <tt>T</tt> is stored contiguously.
 *
 * A sequence container is "contiguous", if its elements are stored in a single contiguous array.  Examples for
 * contiguous sequences are AllocString or char arrays.
 *
 * If an object <tt>obj</tt> is a contiguous sequence, then <tt>begin(obj)</tt> can be converted to a pointer to the
 * first element of the content array.
 */

template <typename T>
struct IsContiguous
{
    typedef False Type;
    enum { VALUE = false };
};

template <typename T>
struct IsContiguous<T const>
    : public IsContiguous<T> {};


//void testStringConcepts()
//{
//    SEQAN_CONCEPT_ASSERT((StringConcept<String<char, Alloc<> > >));
//    SEQAN_CONCEPT_ASSERT((StringConcept<String<Pair<int, double>, Alloc<> > >));
////    SEQAN_CONCEPT_ASSERT((StringConcept<String<bool, Packed<> > >));  // doesn't compile yet
////    SEQAN_CONCEPT_ASSERT((StringConcept<String<Dna5, Packed<> > >));
//    SEQAN_CONCEPT_ASSERT((StringConcept<String<int, Array<50> > >));
//}

/*!
 * @concept ForwardContainerConcept
 * @extends ContainerConcept
 * @headerfile <seqan/basic.h>
 * @brief A forward container is a Container whose elements follow a strict linear ordering.
 *
 * @signature ForwardContainerConcept<T>
 *
 * The order of the elements will not change spontaneously from iteration to iteration.  The linear iteration of
 * elements allows the ForwardContainer to define a lexical ordering if its elements have the
 * @link LessThanComparableConcept @endlink concept.
 */

/*!
 * @concept ReversibleContainerConcept
 * @extends ForwardContainerConcept
 * @headerfile <seqan/basic.h>
 * @brief A reversible container is a forward container that can also be iterated in reverse order.
 *
 * @signature ReversibleContainerConcept<T>
 */

/*!
 * @concept RandomAccessContainerConcept
 * @extends ReversibleContainerConcept
 * @headerfile <seqan/basic.h>
 * @brief A random access container is a reversible container whose iterator is a random access iterator.
 *
 * @signature RandomAccessContainerConcept<T>
 *
 * Random access containers support amortized constant time access to its elements.
 */

/*!
 * @fn RandomAccessContainerConcept::operator[]
 * @brief Returns a reference to an arbitrary element in the sequence.
 *
 * @signature TReference T::operator[](pos);
 *
 * @param[in] pos Position in the sequence (convertible to @link ContainerConcept#Position position @endlink type).
 *
 * @return TReference A reference to into the container with position <tt>pos</tt>.
 */

// TODO(holtgrew): Really deprecated?
/*!
 * @fn RandomAccessContainerConcept#value
 * @brief Global function variant of subscript operator.
 * @deprecated Use the subscript operator (<tt>operator[]</tt>) instead.
 *
 * @signature TReference value(seq, pos);
 *
 * @param[in] seq The sequence to get value in.
 * @param[in] pos Position in the sequence (convertible to @link ContainerConcept#Position position @endlink type).
 *
 * @return TReference A reference to into the container with position <tt>pos</tt>.
 */

/*!
 * @fn RandomAccessContainerConcept#assignValue
 * @brief Assign value in RandomAccessContainer.
 *
 * @signature void assignValue(cont, pos, val);
 *
 * @param[in,out] seq The RandomAccessContainer to modify.
 * @param[in]     pos The position to modify value at.
 * @param[in]     val The value to assign to the given position.
 */

// TODO(holtgrew): Really deprecated?
/*!
 * @fn RandomAccessContainerConcept#getValue
 * @brief Get-value retrieval from container.
 * @deprecated Use the subscript operator (<tt>operator[]</tt>) instead.
 *
 * @signature TGetValue getValue(seq, pos);
 *
 * @param[in] seq The sequence to get value in.
 * @param[in] pos Position in the sequence (convertible to @link ContainerConcept#Position position @endlink type).
 *
 * @return TGetValue The get-value (type is @link ContainerConcept#GetValue @endlink of the sequence type).
 */

/*!
 * @concept StringConcept
 * @brief Sequences are dense linear containers that have positions.
 * @extends RandomAccessContainerConcept
 * @headerfile <seqan/basic.h>
 *
 * @signature StringConcept<T>
 */

/*!
 * @fn StringConcept#iter
 * @headerfile <seqan/sequence.h>
 * @brief Iterator to the item at the given position in a container.
 *
 * @signature TIterator iter(seq, pos[, tag]);
 *
 * @param[in] seq    The sequence to get the iterator for.
 * @param[in] pos    The position to get the iterator for.
 * @param[in] tag    The tag to pick the type of the iterator.
 *
 * @return TIterator The resulting iterator.  If <tt>TTag</tt> is the type of <tt>tag</tt> and <tt>TSequence</tt> the
 *                   type of <tt>seq</tt> then TIterator is of the type <tt>Iterator&lt;TSequence,
 *                   TTag&gt;::Type</tt>.
 *
 * @section Remarks
 *
 * If <tt>pos</tt> is out of range then the iterator is invalid.
 */

/*!
 * @fn StringConcept#append
 * @brief Append a sequence to another one.
 *
 * @signature void append(seq, other);
 *
 * @param[in,out] seq   The sequence to append the other sequence to.
 * @param[in]     other The other sequence to append to <tt>seq</tt>.  Of same type as <tt>seq</tt>.
 */

/*!
 * @fn StringConcept#appendValue
 * @brief Append a value to a sequence.
 *
 * @signature void appendValue(seq, val);
 *
 * @param[in,out] seq The sequence to append a value to (type <tt>TSequence</tt>).
 * @param[in]     val A value to append to the sequence.  Convertible to <tt>Value&lt;TSequence&gt;::Type</tt>.
 */

/*!
 * @fn StringConcept#front
 * @brief Return reference to the first element.
 *
 * @signature TReference front(seq);
 *
 * @param[in] seq The sequence to get the first element of.
 *
 * @return TReference A reference to the first element of <tt>seq</tt>.
 */

/*!
 * @fn StringConcept#back
 * @brief Return reference to the last element.
 *
 * @signature TReference back(seq);
 *
 * @param[in] seq The sequence to get the last element of.
 *
 * @return TReference A reference to the last element of <tt>seq</tt>.
 */

/*!
 * @fn StringConcept#resize
 * @brief Resize a sequence.
 *
 * @signature void resize(seq, len[, val]);
 *
 * @param[in,out] seq Sequence to resize.
 * @param[in]     len Length to resize <tt>seq</tt> to.
 * @param[in]     val When increasing the size, <tt>val</tt> is used to fill new entries.  When omitted,
 *                    <tt>TValue()</tt> is used where <tt>TValue</tt> is the @link ContainerConcept#Value @endlink
 *                    type of the sequence.
 */

/*!
 * @fn StringConcept#clear
 * @brief Remove all elements from the sequences.
 *
 * @signature void clear(seq);
 *
 * @param[in,out] seq Sequence to clear.
 */

/*!
 * @fn StringConcept#erase
 * @brief Erase an element or a range of elements from a sequence.
 *
 * @signature void erase(seq, pos[, posEnd)
 *
 * @param[in,out] seq    Sequence to remove range from.
 * @param[in]     pos    Begin position of the range to remove.
 * @param[in]     posEnd Optional end position of the range to remove.  If omitted, <tt>pos + 1</tt> is used.
 */

/*!
 * @fn StringConcept#eraseFront
 * @brief Erase first element in a sequence.
 *
 * @signature void eraseFront(seq);
 *
 * @param[in,out] seq The sequence to remove the first element from.
 */

/*!
 * @fn StringConcept#eraseBack
 * @brief Erase last element in a sequence.
 *
 * @signature void eraseBack(seq);
 *
 * @param[in,out] seq The sequence to remove the last element from.
 */

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_CONTAINER_CONCEPT_H_
