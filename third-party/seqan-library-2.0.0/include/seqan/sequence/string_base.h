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
// Tags, declarations and generic code for the String class and its
// specializations.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_ARRAY_BASE_H_
#define SEQAN_SEQUENCE_STRING_ARRAY_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSpec = void>
struct Alloc {};

// TODO(holtgrew): This requires some work: explain it, maybe rather put this into a group since the text object appears in no function's signatures.

/*!
 * @concept TextConcept
 * @brief Concept for a type that can be as text of an index.
 * @headerfile <seqan/sequence.h>
 *
 * @signature concept TextConcept;
 *
 * Certain algorithms and data structures can work for both strings and string sets but need to treat these two
 * types slightly different.  Examples are index data structures and algorithms that build the indices and use
 * the indices for lookup.
 *
 * To facilitate writing of generic algorithms, the TextConcept concept gives a common interface to both for this
 * kind of algorithms.
 *
 * @see String
 * @see StringSet
 */

/*!
 * @mfn TextConcept#StringSetLimits
 * @brief Return type of string set limits for TextConcept types.
 *
 * @signature StringSetLimits<TText>::Type;
 *
 * @tparam TText The type of the text.
 *
 * @return Type The type of string set limits objects.
 */

/*!
 * @mfn TextConcept#SAValue
 * @brief The default alphabet type of a suffix array, i.e. the type to store a
 *        position of a string or string set.
 *
 * @signature SAValue<TText>::Type
 *
 * @tparam TText The text type to query.
 *
 * @return TReturn A type to store a position in a <tt>TText</tt>.  This could be an integer for strings or a
 *                 pair of integers for string sets.
 *
 * @section Usage
 *
 * This type should be removed for functions returning positions in texts such as online or index-based search.
 * Thus, always use this metafunction for declaring position variables.
 *
 * Use the functions @link TextConcept#posLocalize @endlink, @link TextConcept#posGlobalize @endlink, @link
 * TextConcept#getSeqNo @endlink, and @link TextConcept#getSeqOffset @endlink for conversion between local
 * and global positions in texts.
 *
 * @section Examples
 *
 * The following shows the original definition of the SAValue metafunction in SeqAn.
 *
 * @code{.cpp}
 * template <typename TString, typename TSpec>
 * struct SAValue<StringSet<TString, TSpec> >
 * {
 *     typedef Pair<
 *             typename Size<StringSet<TString, TSpec> >::Type,
 *             typename SAValue<TString>::Type,
 *             Pack
 *         > Type;
 * };
 * @endcode
 */

/*!
 * @fn TextConcept#stringSetLimits
 * @brief Return string delimiter positions for TextConcept types.
 *
 * @signature TStringSetLimits stringSetLimits(text);
 *
 * @param[in] text The text to query for its string set limits.
 *
 * @return TStringSetLimits The string set limits (of type @link TextConcept#StringSetLimits @endlink).
 */

/*!
 * @fn TextConcept#posLocalToX
 * @brief Converts a local to a local/global position.
 *
 * @signature void posLocalToX(dst, localPos, limits);
 *
 * @param[in] dst      The local or global position (pair or integer value) is written here.
 * @param[in] localPos The local position.
 * @param[in] limits   The string limits as returned by @link TextConcept#stringSetLimits @endlink.
 */

/*!
 * @class String
 * @implements StringConcept
 * @implements TextConcept
 * @implements SegmentableConcept
 * @headerfile <seqan/sequence.h>
 * @brief @link StringConcept Sequence @endlink container class.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class String<TValue, TSpec>;
 *
 * @tparam TValue The element type of the string.
 * @tparam TSpec  The tag for selecting the string specialization.
 *
 * The String class is for storing sequences and thus at the core of the sequence analysis library SeqAn.  They
 * are models for the @link StringConcept sequence concept @endlink but extend the sequence concept by allowing
 * implicit conversion of other sequence into strings as long as the element conversion works:
 *
 * @snippet demos/sequence/string.cpp initializing strings
 *
 * Aside from that, the usual operations (appending, insertion, removing, element access) are available as well.
 *
 * @snippet demos/sequence/string.cpp usual operations
 *
 * Strings have a size (the actual number of elements) and a capacity (the number of elements that memory has been
 * allocated for).  Note that clearing a string does not free the memory (as the STL, SeqAn assumes that strings will
 * later require a similar amount of memory as before).  Using @link ContainerConcept#shrinkToFit @endlink, the user can
 * force a re-allocation of the memory such that the string afterward uses the minimal amount of memory to accomodate
 * all of its objects.
 *
 * @snippet demos/sequence/string.cpp clear and resize
 *
 * @section Examples
 *
 * This example shows a brute force pattern matching scheme for two character Strings.  Creation of String "text" shows
 * the usage of some available String operating functions.  See class @link StringSet @endlink for an example of a
 * String container with other than simple type values.  See class @link Index @endlink example for efficiently finding
 * the same pattern matches using an index.
 *
 * @include demos/sequence/string2.cpp
 *
 * The output is as follows:
 *
 * @include demos/sequence/string2.cpp.stdout
 *
 * @see StringSet
 */

/*!
 * @fn String::String
 * @brief Constructor.
 *
 * @signature String::String()
 * @signature String::String(other)
 *
 * @param[in] other The source for the copy constructor.  Can be of any @link StringConcept sequence @endlink type
 *                  as long as <tt>other</tt>'s elements are convertible to the value type of this string.
 *
 * Default and copy constructor are implemented.
 */

/*!
 * @fn String::operator=
 * @brief The String assignment operator allows assignment of convertible sequences.
 *
 * @signature TString String::operator=(other)
 *
 * @param[in] other The other string.  Must be a sequence whose elements are convertible into this String's type.
 *
 * @return TString Reference to the String objecta after assignment.
 */

// TODO(holtgrew): The conversion functions rather belong into their own group than to the concept. The original documentation was a bit misleading and needs to be updated.

/*!
 * @fn TextConcept#posLocalize
 * @brief Converts a local/global to a local position.
 * @headerfile <seqan/sequence.h>
 *
 * @signature void posLocalize(result, pos, limits)
 *
 * @param[in] pos    A local or global position (pair or integer value).
 * @param[in] limits The limits string returned by @link TextConcept#stringSetLimits @endlink.
 * @param[in] result Reference to the resulting corresponding local position of
 *                   <tt>pos</tt>.
 */

/*!
 * @fn TextConcept#posGlobalize
 * @brief Converts a local/global to a global position.
 * @headerfile <seqan/sequence.h>
 *
 * @signature TPos posGlobalize(pos, limits)
 *
 * @param[in] pos A local or global position (pair or integer value). Types: Pair
 * @param[in] limits The limits string returned by @link TextConcept#stringSetLimits @endlink.
 *
 * @return TPos The corresponding global position of <tt>pos</tt>. If
 *                 <tt>pos</tt> is an integral type <tt>pos</tt> is returned. If
 *                 not, <tt>limits[getSeqNo(pos, limits)] + getSeqOffset(pos,
 *                 limits)</tt> is returned.
 */

/*!
 * @fn TextConcept#getSeqNo
 * @brief Returns the sequence number of a position.
 * @headerfile <seqan/sequence.h>
 *
 * @signature TSeqNo getSeqNo(pos[, limits]);
 *
 * @param[in] pos A position. Types: Pair
 * @param[in] limits The limits string returned by @link TextConcept#stringSetLimits @endlink.
 *
 * @return TSeqNo A single integer value that identifies the string within the stringset <tt>pos</tt> points at. If
 *                <tt>limits</tt> is omitted or @link Nothing @endlink <tt>getSeqNo</tt> returns 0.If <tt>pos</tt> is a
 *                local position (of class @link Pair @endlink) then <tt>i1</tt> is returned.If <tt>pos</tt> is a global
 *                position (integer type and <tt>limits</tt> is a @link String @endlink) then <tt>pos</tt> is converted
 *                to a local position and <tt>i1</tt> is returned.
 */

/*!
 * @fn TextConcept#getSeqOffset
 * @brief Returns the local sequence offset of a position.
 * @headerfile <seqan/sequence.h>
 *
 * @signature TOffset getSeqOffset(pos[, limits]);
 *
 * @param[in] pos A position. Types: Pair
 * @param[in] limits The limits string returned by @link TextConcept#stringSetLimits @endlink.
 *
 * @return TOffset A single integer value that identifies the position within the string <tt>pos</tt> points at.If
 *                 <tt>limits</tt> is omitted or @link Nothing @endlink <tt>getSeqNo</tt> returns <tt>pos</tt>.  If
 *                 <tt>pos</tt> is a local position (of class @link Pair @endlink) then <tt>i2</tt> is returned.If
 *                 <tt>pos</tt> is a global position (integer type and <tt>limits</tt> is a @link String @endlink) then
 *                 <tt>pos</tt> is converted to a local position and <tt>i2</tt> is returned.
 */

template <typename TValue, typename TSpec = Alloc<> >
class String;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Value<String<TValue, TSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Value<String<TValue, TSpec> const >
        : public Value<String<TValue, TSpec> >
{};

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Spec<String<TValue, TSpec> >
{
    typedef TSpec Type;
};
template <typename TValue, typename TSpec>
struct Spec<String<TValue, TSpec> const>:
    public Spec<String<TValue, TSpec> >
{
};

// ----------------------------------------------------------------------------
// Metafunction StringSpec
// ----------------------------------------------------------------------------

template <typename T>
struct StringSpec
{
    typedef Alloc<> Type;
};

template <typename T>
struct StringSpec<T const> : StringSpec<T> {};

template <typename TValue, typename TSpec>
struct StringSpec<String<TValue, TSpec> >
{
    typedef TSpec   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Chunk
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct Chunk<String<TValue, TSpec> >:
    Chunk<typename Iterator<String<TValue, TSpec>, Rooted>::Type> {};

// ----------------------------------------------------------------------------
// Metafunction IsSequence
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct IsSequence<String<TValue, TSpec> > : True {};

// ----------------------------------------------------------------------------
// Concept StringConcept
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
SEQAN_CONCEPT_IMPL((String<TValue, TSpec>), (StringConcept));         // resizable container

template <typename TValue, typename TSpec>
SEQAN_CONCEPT_IMPL((String<TValue, TSpec> const), (ContainerConcept));  // read-only container

// ----------------------------------------------------------------------------
// Internal Metafunction TempCopy_
// ----------------------------------------------------------------------------

template <typename T>
struct TempCopy_
{
    typedef typename Value<T>::Type TValue_;
    typedef typename RemoveConst_<TValue_>::Type TValueNotConst_;
    typedef String<TValueNotConst_, Alloc<> > Type;
};

// ============================================================================
// Functions
// ============================================================================

// TODO(holtgrew): Where to move this documentation/specification-only stuff?

// ----------------------------------------------------------------------------
// Function swap()
// ----------------------------------------------------------------------------

template <typename TAlphabet, typename TSpec>
inline void
swap(String<TAlphabet, TSpec> & left,
     String<TAlphabet, TSpec> & right)
{

    typedef String<TAlphabet, TSpec> TString;

    TString tmp(left, Move());
    move(left, right);
    move(right, tmp);
}

// ----------------------------------------------------------------------------
// Function shareResources()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline bool
shareResources(String<TValue, TSpec> const & obj1,
               TValue const & obj2)
{
    return (begin(obj1) >= &obj2) && (end(obj1) <= &obj2);
}

template <typename TValue, typename TSpec>
inline bool
shareResources(TValue const & obj1,
               String<TValue, TSpec> const & obj2)
{
    return (begin(obj2) >= &obj1) && (end(obj2) <= &obj1);
}

// TODO(holtgrew): Where to move this documentation/specification-only stuff?
// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TPos>
SEQAN_HOST_DEVICE inline typename Reference< String<TValue, TSpec> >::Type
value(String<TValue, TSpec> & me,
      TPos const & pos)
{
    typedef typename Position< String<TValue, TSpec> >::Type TStringPos SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT_MSG(static_cast<TStringPos>(pos), static_cast<TStringPos>(length(me)), "Trying to access an element behind the last one!");
    return *(begin(me, Standard()) + pos);
}

template <typename TValue, typename TSpec, typename TPos>
SEQAN_HOST_DEVICE inline typename Reference< String<TValue, TSpec> const >::Type
value(String<TValue, TSpec> const & me,
      TPos const & pos)
{
    typedef typename Position< String<TValue, TSpec> const >::Type TStringPos SEQAN_TYPEDEF_FOR_DEBUG;
    SEQAN_ASSERT_LT_MSG(static_cast<TStringPos>(pos), static_cast<TStringPos>(length(me)), "Trying to access an element behind the last one!");
    return *(begin(me, Standard()) + pos);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
SEQAN_HOST_DEVICE inline typename Size< String<TValue, TSpec> const>::Type
length(String<TValue, TSpec> const & me)
{
    return end(me, Standard()) - begin(me, Standard());
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
SEQAN_HOST_DEVICE inline bool
empty(String<TValue, TSpec> const & me)
{
    return end(me, Standard()) == begin(me, Standard());
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void
clear(String<TValue, TSpec> & me)
{
    arrayDestruct(begin(me, Standard()), end(me, Standard()));
    _setLength(me, 0);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TExpand>
struct ClearSpaceStringBase_
{
};

// ----------------------------------------------------------------------------
// Internal Function _clearSpace()
// ----------------------------------------------------------------------------

template <>
struct ClearSpaceStringBase_<Insist>
{
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size)
    {
        arrayDestruct(begin(seq, Standard()), end(seq, Standard()));
        _setLength(seq, size);
        return size;
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type limit)
    {
        arrayDestruct(begin(seq, Standard()), end(seq, Standard()));
        if (limit < size)
        {
            size = limit;
        }
        _setLength(seq, size);
        return size;
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end)
    {
        typename Size<T>::Type new_length = length(seq) + size - (end - start);
        arrayClearSpace(begin(seq, Standard()) + start, length(seq) - start, end - start, size);
        _setLength(seq, new_length);
        return size;
    }

    template <typename T>
    static typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Value<T>::Type * seq_buffer = begin(seq);
        typename Size<T>::Type seq_length = length(seq);

        if (limit > start + size)
        {
            typename Size<T>::Type removed_size = end - start;
            typename Size<T>::Type new_length = seq_length - removed_size + size;
            if (limit < new_length)
            {
                arrayDestruct(seq_buffer + limit, seq_buffer + new_length);
                seq_length -= new_length - limit;
            }
            arrayClearSpace(seq_buffer + start, seq_length - start, end - start, size);
            _setLength(seq, new_length);
            return size;
        }
        else
        {
            arrayDestruct(seq_buffer + start, seq_buffer + seq_length);
            _setLength(seq, limit);
            if (limit > start) return limit - start;
            else return 0;
        }
    }
/*
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, Insist());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, limit, Insist());
    }
*/
};


template <>
struct ClearSpaceStringBase_<Limit>
{

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size)
    {
        return _clearSpace(seq, size, capacity(seq), Insist());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type limit)
    {
        typename Size<T>::Type seq_capacity = capacity(seq);
        if (limit > seq_capacity)
        {
            limit = seq_capacity;
        }
        return _clearSpace(seq, size, limit, Insist());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end)
    {
        return _clearSpace(seq, size, start, end, capacity(seq), Insist());
    }

    template <typename T>
    static typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Size<T>::Type seq_capacity = capacity(seq);
        if (limit > seq_capacity)
        {
            limit = seq_capacity;
        }
        return _clearSpace(seq, size, start, end, limit, Insist());
    }

/*
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, Insist());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, limit, Insist());
    }
*/
};

template <typename TExpand>
struct ClearSpaceExpandStringBase_
{
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size)
    {
        arrayDestruct(begin(seq, Standard()), end(seq, Standard()));
        typename Size<T>::Type old_capacity = capacity(seq);
        typename Value<T>::Type * old_array = _reallocateStorage(seq, size, TExpand());
        if (old_array)
        {
            _deallocateStorage(seq, old_array, old_capacity);
        }
        _setLength(seq, size);
        return size;
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type limit)
    {
        arrayDestruct(begin(seq, Standard()), end(seq, Standard()));
        if (limit < size)
        {
            size = limit;
        }
        typename Size<T>::Type old_capacity = capacity(seq);
        typename Value<T>::Type * old_array = _reallocateStorage(seq, size, limit, TExpand());
        if (old_array)
        {
            _deallocateStorage(seq, old_array, old_capacity);
        }
        _setLength(seq, size);
        return size;
    }

    template <typename T>
    static typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end)
    {
        typename Size<T>::Type old_length = length(seq);
        typename Size<T>::Type removed_size = end - start;
        typename Size<T>::Type new_length = old_length - removed_size + size;

        SEQAN_ASSERT_LEQ(start, end);
        SEQAN_ASSERT_LEQ(start, old_length);

        typename Size<T>::Type old_capacity = capacity(seq);
        typename Value<T>::Type * old_array = _reallocateStorage(seq, new_length, TExpand());
        typename Value<T>::Type * seq_array = begin(seq);

        if (old_array)
        {
            arrayConstructMove(old_array, old_array + start, seq_array);
            arrayConstructMove(old_array + end, old_array + old_length, seq_array + start + size);
            _deallocateStorage(seq, old_array, old_capacity);
            // TODO(weese:) Did we miss to destruct the old_array here, e.g. with arrayDestruct
        }
        else
        {
            arrayClearSpace(seq_array + start, old_length - start, removed_size, size);
        }

        _setLength(seq, new_length);

        return size;
    }

    template <typename T>
    static typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Size<T>::Type old_length = length(seq);
        typename Size<T>::Type removed_size = end - start;
        typename Size<T>::Type need_length = old_length - removed_size + size;

        SEQAN_ASSERT_LEQ(start, end);
        SEQAN_ASSERT_LEQ(start, old_length);

        typename Size<T>::Type new_length = need_length;
        typename Size<T>::Type length_to_copy = old_length;
        if (limit < need_length)
        {
            new_length = limit;
            length_to_copy = new_length - size + removed_size;
        }

        bool keep_second_part = (new_length > start + size);

        typename Size<T>::Type old_capacity = capacity(seq);
        typename Value<T>::Type * old_array = _reallocateStorage(seq, new_length, limit, TExpand());
        typename Value<T>::Type * seq_array = begin(seq);

        if (old_array)
        {//new buffer allocated
            typename Size<T>::Type keep_start_length = (start > new_length) ? new_length : start;
            arrayConstructMove(old_array, old_array + keep_start_length, seq_array);
            if (keep_second_part)
            {
                arrayConstructMove(old_array + end, old_array + length_to_copy, seq_array + start + size);
            }
            _deallocateStorage(seq, old_array, old_capacity);
        }
        else
        {
            if (keep_second_part)
            {
                arrayClearSpace(seq_array + start, length_to_copy - start, end - start, size);
                if (length_to_copy < old_length)
                {
                    arrayDestruct(seq_array + length_to_copy, seq_array + old_length);
                }
            }
            else
            {
                arrayDestruct(seq_array + start, seq_array + old_length);
            }
        }

        _setLength(seq, new_length);

        if (keep_second_part) return size;
        else if (new_length > start) return new_length - start;
        else return 0;
    }

/*
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, TExpand());
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Iterator<T>::Type start,
        typename Iterator<T>::Type end,
        typename Size<T>::Type limit)
    {
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, limit, TExpand());
    }
*/
};

template <>
struct ClearSpaceStringBase_<Exact>:
    ClearSpaceExpandStringBase_<Exact>
{
};

template <>
struct ClearSpaceStringBase_<Generous>:
    ClearSpaceExpandStringBase_<Generous>
{
};

template<typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
_clearSpace(String<TValue, TSpec> & me,
        TSize size,
        Tag<TExpand>)
{
    return ClearSpaceStringBase_<Tag<TExpand> >::_clearSpace_(me, size);
}

template<typename TValue, typename TSpec, typename TSize, typename TCapacity, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
_clearSpace(String<TValue, TSpec> & me,
        TSize size,
        TCapacity limit,
        Tag<TExpand>)
{
    return ClearSpaceStringBase_<Tag<TExpand> >::_clearSpace_(me, size, limit);
}

template<typename TValue, typename TSpec, typename TSize, typename TPosition, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
_clearSpace(String<TValue, TSpec> & me,
            TSize size,
            TPosition pos_begin,
            TPosition pos_end,
            Tag<TExpand>)
{
    return ClearSpaceStringBase_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end);
}

template<typename TValue, typename TSpec, typename TSize, typename TPosition, typename TCapacity, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
_clearSpace(String<TValue, TSpec> & me,
            TSize size,
            TPosition pos_begin,
            TPosition pos_end,
            TCapacity limit,
            Tag<TExpand>)
{
    return ClearSpaceStringBase_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end, limit);
}

// ----------------------------------------------------------------------------
// Function resizeSpace()
// ----------------------------------------------------------------------------

template<typename TValue, typename TSpec, typename TSize, typename TBeginPosition, typename TEndPosition, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
resizeSpace(String<TValue, TSpec> & me,
            TSize size,
            TBeginPosition pos_begin,
            TEndPosition pos_end,
            Tag<TExpand> tag)
{
    typedef typename Size< String<TValue, TSpec> >::Type TSize_;
    typedef typename Position<String<TValue, TSpec> >::Type TPos_;
    TSize_ ret_ = _clearSpace(
        me,
        static_cast<TSize_>(size),
        static_cast<TPos_>(pos_begin),
        static_cast<TPos_>(pos_end),
        tag);
    arrayConstruct(iter(me, pos_begin), iter(me, pos_begin) + ret_);
    return ret_;
}

template<typename TValue, typename TSpec, typename TSize, typename TBeginPosition, typename TEndPosition, typename TCapacity, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
resizeSpace(String<TValue, TSpec> & me,
            TSize size,
            TBeginPosition pos_begin,
            TEndPosition pos_end,
            TCapacity limit,
            Tag<TExpand> tag)
{
    typedef typename Size< String<TValue, TSpec> >::Type TSize_;
    typedef typename Position<String<TValue, TSpec> >::Type TPos_;
    TSize_ ret_ = _clearSpace(
        me,
        static_cast<TSize_>(size),
        static_cast<TPos_>(pos_begin),
        static_cast<TPos_>(pos_end),
        static_cast<TSize_>(limit),
        tag);
    arrayConstruct(iter(me, pos_begin), iter(me, pos_begin) + ret_);
    return ret_;
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

// Facade version without overflow tag.  Forwards to version with overflow
// tag, using Metafunction.DefaultOverflowImplicity.

template<typename TTargetValue, typename TTargetSpec, typename TSource>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
       TSource & source)
{
    typedef String<TTargetValue, TTargetSpec> TTarget;
    assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTargetValue, typename TTargetSpec, typename TSource>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
       TSource const & source)
{
    typedef String<TTargetValue, TTargetSpec> TTarget;
    assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

// Helper struct for assigning strings.

template <typename TExpand>
struct AssignString_
{
    template <typename TTarget, typename TSource>
    static inline void
    assign_(
        TTarget & target,
        TSource & source)
    {
        if (empty(source) && empty(target))
            return;  // Do nothing if both source and target are empty.
        if (!getObjectId(source) || !shareResources(target, source))
        {
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()));
        }
        else
        {
            if ((void *) &target == (void *) &source) return;

            typename TempCopy_<TSource>::Type temp(source, length(source));
            assign(target, temp, TExpand());
        }
    }

    template <typename TTarget, typename TSource>
    static inline void
    assign_(
        TTarget & target,
        TSource & source,
        typename Size<TTarget>::Type limit)
    {
        if (!getObjectId(source) || !shareResources(target, source))
        {
            typename Size<TTarget>::Type part_length = _clearSpace(target, typename Size<TTarget>::Type(length(source)), limit, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()));
        }
        else
        {
            if ((void *) &target == (void *) &source) return;

            typename Size<TTarget>::Type source_length = length(source);
            if (source_length > limit) source_length = limit;

            typename TempCopy_<TSource>::Type temp(source, source_length);
            assign(target, temp, TExpand());
        }
    }
};

// Interface for assign with overflow strategy tag.  Forwards to static
// functions of helper struct.

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TExpand>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
       TSource const & source,
       Tag<TExpand>)
{
    AssignString_<Tag<TExpand> >::assign_(target, source);
}
template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TSize, typename TExpand>
inline void
assign(String<TTargetValue, TTargetSpec> & target,
       TSource const & source,
       TSize limit,
       Tag<TExpand>)
{
    AssignString_<Tag<TExpand> >::assign_(target, source, limit);
}

//// TODO(holtgrew): Still required with dropped VC++ 2003 support?
////this variant is a workaround for the "const array"-bug of VC++
//
//template<typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TExpand>
//inline void
//assign(String<TTargetValue, TTargetSpec> & target,
//       TSourceValue const * source,
//       Tag<TExpand>)
//{
//    AssignString_<Tag<TExpand> >::assign_(target, source);
//}
//template<typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSize, typename TExpand>
//inline void
//assign(String<TTargetValue, TTargetSpec> & target,
//       TSourceValue const * source,
//       TSize limit,
//       Tag<TExpand>)
//{
//    AssignString_<Tag<TExpand> >::assign_(target, source, limit);
//}

// ----------------------------------------------------------------------------
// Function move()
// ----------------------------------------------------------------------------

//implementation of move for contiguous sequences
//note: there is a problem, if sizeof(TSourceValue) and sizeof(TTargetValue) are not a multiple
//  of each other, since in this case the correct size cannot be determined afterwards
//  when calling the deallocate function.
//  ???TODO
template <typename TTarget, typename TSource>
void
_moveContiguous(TTarget & target,
                TSource & source)
{
    typedef typename Value<TSource>::Type TSourceValue;
    typedef typename Value<TTarget>::Type TTargetValue;

    clear(target);

    typename Iterator<TSource, Standard>::Type source_begin = begin(source, Standard());
    typename Iterator<TTarget, Standard>::Type target_begin = (typename Iterator<TTarget, Standard>::Type) begin(source, Standard());

    typename Size<TTarget>::Type size = sizeof(TSourceValue) * capacity(source);
    if (size >= sizeof(TTargetValue))
    {
        if (sizeof(TSourceValue) <= 2) ++size; //regard the "end of string termination" case
        typename Size<TTarget>::Type target_capacity = size / sizeof(TTargetValue);
        if (sizeof(TTargetValue) <= 2) --target_capacity; //regard the "end of string termination" case

        typename Size<TTarget>::Type target_length = length(source);
        if (target_length > target_capacity)
        {
            target_length = target_capacity;
        }

        if (sizeof(TSourceValue) >= sizeof(TTargetValue))
        {
            arrayMoveForward(source_begin, source_begin + target_length, target_begin);
        }
        else
        {
            arrayMoveBackward(source_begin, source_begin + target_length, target_begin);
        }

        _setBegin(target, target_begin);
        _setLength(target, target_length);
        _setCapacity(target, target_capacity);

        typedef typename Iterator<TSource, Standard>::Type TSourceIterator;
        _setBegin(source, TSourceIterator(0));
        _setLength(source, 0);
        _setCapacity(source, 0);
    }
    else
    {
        clear(source);
    }
}

template<typename TTargetValue, typename TTargetSpec, typename TSource>
inline void
move(String<TTargetValue, TTargetSpec> & target,
     TSource & source)
{
    typedef String<TTargetValue, TTargetSpec> TTarget;
    move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTargetValue, typename TTargetSpec, typename TSource>
inline void
move(String<TTargetValue, TTargetSpec> & target,
     TSource const & source)
{
    typedef String<TTargetValue, TTargetSpec> TTarget;
    move(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TTag>
inline void
move(String<TTargetValue, TTargetSpec> & target,
     TSource & source,
     Tag<TTag> const & tag)
{
    assign(target, source, tag);
}

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TTag>
inline void
move(String<TTargetValue, TTargetSpec> & target,
     TSource const & source,
     Tag<TTag> const & tag)
{
    assign(target, source, tag);
}

// TODO(holtgrew): Garbage?

//////////////////////////////////////////////////////////////////////////////
// valueConstructMove:
// it is usually better for strings to default construct and move instead of
// copy construct strings

/*
template <typename TIterator, typename TValue, typename TSpec>
inline void
valueConstructMove(TIterator it,
                   String<TValue, TSpec> const & value)
{
    valueConstruct(it);
    move(*it, value);
}
*/

// ----------------------------------------------------------------------------
// Function _stringCheckForOverlap
// ----------------------------------------------------------------------------

template <typename TIter1, typename TIter2, typename TSize>
inline bool
_stringCheckForPossibleOverlap(TIter1 const &, TIter2 const &, TSize)
{
    return true;
}

template <typename TIter, typename TSize>
inline bool
_stringCheckForPossibleOverlap(TIter const &it1, TIter const &it2, TSize length)
{
    // return false if [it1,it1+length) and [it2,it2+length) are not overlapping
    return !(it2 + length <= it1 || it1 + length <= it2);
}

// ----------------------------------------------------------------------------
// Function append()
// ----------------------------------------------------------------------------

template <typename TExpand>
struct AppendString_
{
    template <typename TTarget, typename TSource>
    static inline void
    append_(TTarget & target,
            TSource & source)
    {
        if (!getObjectId(source) || !shareResources(target, source) ||
            !_stringCheckForPossibleOverlap(begin(source, Standard()), end(target, Standard()), length(source)))
        {
            typename Size<TTarget>::Type target_length = length(target);
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), target_length, target_length, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()) + target_length);
        }
        else
        {
            typename TempCopy_<TSource>::Type temp(source, length(source));
            append(target, temp, TExpand());
        }
    }

    template <typename TTarget, typename TSource>
    static inline void
    append_(TTarget & target,
            TSource & source,
            typename Size<TTarget>::Type limit)
    {
        typename Iterator<TTarget, Standard>::Type target_begin = begin(target, Standard());
        if (!getObjectId(source) || !shareResources(target, source))
        {
            typename Size<TTarget>::Type target_length = length(target);
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), target_length, target_length, limit, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()) + target_length);
        }
        else
        {
            typename Size<TTarget>::Type target_length = length(target);
            if (target_length >= limit)
            {
                arrayDestruct(target_begin + limit, target_begin + target_length);
                _setLength(target, limit);
            }
            else
            {
                limit -= target_length;
                typename Size<TTarget>::Type source_length = length(source) ;
                if (source_length > limit) source_length = limit;

                typename TempCopy_<TSource>::Type temp(source, source_length);
                append(target, temp, TExpand());
            }
        }
    }
};

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TExpand>
inline void
append(String<TTargetValue, TTargetSpec> & target,
       TSource const & source,
       Tag<TExpand>)
{
    AppendString_<Tag<TExpand> >::append_(target, source);
}

template<typename TTargetValue, typename TTargetSpec, typename TSource, typename TExpand>
inline void
append(String<TTargetValue, TTargetSpec> & target,
       TSource const & source,
       typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
       Tag<TExpand>)
{
    AppendString_<Tag<TExpand> >::append_(target, source, limit);
}

// TODO(holtgrew): Still required with dropped VC++ 2003 support?
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TExpand>
inline void
append(String<TTargetValue, TTargetSpec> & target,
       TSourceValue * source,
       Tag<TExpand>)
{
    AppendString_<Tag<TExpand> >::append_(target, source);
}

template<typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TExpand>
inline void
append(String<TTargetValue, TTargetSpec> & target,
       TSourceValue * source,
       typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
       Tag<TExpand>)
{
    AppendString_<Tag<TExpand> >::append_(target, source, limit);
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

template <typename TExpand>
struct AppendValueToString_
{
    template <typename T, typename TValue>
    static inline void
    appendValue_(T & me,
                 TValue SEQAN_FORWARD_CARG _value)
    {
        typedef typename Value<T>::Type TTargetValue;
        typedef typename Size<T>::Type TSize;

        TSize me_length = length(me);
        if (capacity(me) <= me_length)
        {
            TTargetValue temp_copy(SEQAN_FORWARD(TValue, _value)); //temp copy because resize could invalidate _value
            // TODO(holtgrew): The resize() function will default construct the last element. This is slow. Get rid of this.
            TSize new_length = reserve(me, me_length + 1, TExpand());
            if (me_length < new_length)
            {
                // *(begin(me) + me_length) = temp_copy;
                valueConstruct(begin(me, Standard()) + me_length, SEQAN_FORWARD(TTargetValue, temp_copy)); //??? this should be valueMoveConstruct
                _setLength(me, me_length + 1);
            }
        }
        else
        {
            valueConstruct(begin(me, Standard()) + me_length, SEQAN_FORWARD(TValue, _value));
            _setLength(me, me_length + 1);
        }
    }
};

template <typename TTargetValue, typename TTargetSpec, typename TValue, typename TExpand>
inline void
appendValue(String<TTargetValue, TTargetSpec> & me,
            TValue SEQAN_FORWARD_CARG _value,
            Tag<TExpand>)
{
    AppendValueToString_<Tag<TExpand> >::appendValue_(me, SEQAN_FORWARD(TValue, _value));
}

// ----------------------------------------------------------------------------
// Function appendValue(Serial)
// ----------------------------------------------------------------------------

template <typename TTargetValue, typename TTargetSpec, typename TValue, typename TExpand>
inline void
appendValue(String<TTargetValue, TTargetSpec> & me,
            TValue SEQAN_FORWARD_CARG _value,
            Tag<TExpand> const & expandTag,
            Serial)
{
    appendValue(me, SEQAN_FORWARD(TValue, _value), expandTag);
}

//////////////////////////////////////////////////////////////////////////////
// insertValue
//////////////////////////////////////////////////////////////////////////////

template <typename TExpand>
struct InsertValueToString_
{
    template <typename T, typename TPosition, typename TValue>
    static inline void
    insertValue_(T & me,
                TPosition pos,
                TValue & _value)
    {
        typename Value<T>::Type temp_copy = _value; //temp copy because resizeSpace could invalidate _value
        resizeSpace(me, 1, pos, pos, TExpand());
        if ((typename Size<T>::Type) pos < length(me))
            moveValue(me, pos, temp_copy);
    }
};

template <typename TTargetValue, typename TTargetSpec, typename TPosition, typename TValue, typename TExpand>
inline void
insertValue(String<TTargetValue, TTargetSpec> & me,
            TPosition pos,
            TValue const & _value,
            Tag<TExpand>)
{
    InsertValueToString_<Tag<TExpand> >::insertValue_(me, pos, _value);
}

// ----------------------------------------------------------------------------
// Function replace()
// ----------------------------------------------------------------------------

template <typename TExpand>
struct ReplaceString_
{
    template <typename TTarget, typename TSource>
    static inline void
    replace_(TTarget & target,
             typename Size<TTarget>::Type pos_begin,
             typename Size<TTarget>::Type pos_end,
             TSource & source)
    {
        if (!getObjectId(source) || !shareResources(target, source))
        {
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), pos_begin, pos_end, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()) + pos_begin);
        }
        else
        {
            typename TempCopy_<TSource>::Type temp(source, length(source));
            replace(target, pos_begin, pos_end, temp, TExpand());
        }
    }

    template <typename TTarget, typename TSource>
    static inline void
    replace_(TTarget & target,
             typename Size<TTarget>::Type pos_begin,
             typename Size<TTarget>::Type pos_end,
             TSource & source,
             typename Size<TTarget>::Type limit)
    {
        if (!getObjectId(source) || !shareResources(target, source))
        {
            typename Size<TTarget>::Type part_length = _clearSpace(target, length(source), pos_begin, pos_end, limit, TExpand());
            arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()) + pos_begin);
        }
        else
        {
            if (pos_begin >= limit)
            {
                arrayDestruct(begin(target) + limit, end(target));
                _setLength(target, limit);
            }
            else
            {
                limit -= pos_begin;
                typename Size<TTarget>::Type source_length = length(source) ;
                if (source_length > limit) source_length = limit;

                typename TempCopy_<TSource>::Type temp(source, source_length);
                replace(target, pos_begin, pos_end, temp, limit, TExpand());
            }
        }
    }

};

template<typename TTargetValue, typename TTargetSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(String<TTargetValue, TTargetSpec> & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source,
        Tag<TExpand>)
{
    ReplaceString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source);
}

template<typename TTargetValue, typename TTargetSpec, typename TPositionBegin, typename TPositionEnd, typename TSource, typename TExpand>
inline void
replace(String<TTargetValue, TTargetSpec> & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSource const & source,
        typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
        Tag<TExpand>)
{
    ReplaceString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source, limit);
}

// TODO(holtgrew): Still required with dropped VC++ 2003 support?
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TTargetSpec, typename TPositionBegin, typename TPositionEnd, typename TSourceValue, typename TExpand>
inline void
replace(String<TTargetValue, TTargetSpec> & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSourceValue const * source,
        Tag<TExpand>)
{
    ReplaceString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source);
}

template<typename TTargetValue, typename TTargetSpec, typename TPositionBegin, typename TPositionEnd, typename TSourceValue, typename TExpand>
inline void
replace(String<TTargetValue, TTargetSpec> & target,
        TPositionBegin pos_begin,
        TPositionEnd pos_end,
        TSourceValue const * source,
        typename Size< String<TTargetValue, TTargetSpec> >::Type limit,
        Tag<TExpand>)
{
    ReplaceString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source, limit);
}

// ----------------------------------------------------------------------------
// Internal Function _reallocateStorage()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity)
{
    return _allocateStorage(me, new_capacity);
}

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity,
    Exact)
{
    typedef typename Size<String<TValue, TSpec> >::Type TStringSize;
    if (static_cast<TStringSize>(new_capacity) <= capacity(me))
        return 0;

    return _reallocateStorage(me, new_capacity);
}

template <typename TValue, typename TSpec, typename TSize, typename TSize2>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity,
    TSize2 limit,
    Exact)
{
    typedef typename Size<String<TValue, TSpec> >::Type TStringSize;
    if (static_cast<TStringSize>(new_capacity) <= capacity(me))
        return 0;

    if (new_capacity > limit) new_capacity = limit;
    return _reallocateStorage(me, new_capacity);
}

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity,
    Generous)
{
    typedef typename Size<String<TValue, TSpec> >::Type TStringSize;
    if (static_cast<TStringSize>(new_capacity) <= capacity(me))
        return 0;

    new_capacity = computeGenerousCapacity(me, new_capacity);
    return _reallocateStorage(me, new_capacity);
}

template <typename TValue, typename TSpec, typename TSize, typename TSize2>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> & me,
    TSize new_capacity,
    TSize2 limit,
    Generous)
{
    typedef typename Size<String<TValue, TSpec> >::Type TStringSize;
    if (static_cast<TStringSize>(new_capacity) <= capacity(me))
        return 0;

    new_capacity = computeGenerousCapacity(me, new_capacity);
    if (new_capacity > limit) new_capacity = limit;
    return _reallocateStorage(me, new_capacity);
}

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> &,
    TSize,
    Insist)
{
    return 0;
}

template <typename TValue, typename TSpec, typename TSize, typename TSize2>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> &,
    TSize,
    TSize2,
    Insist)
{
    return 0;
}

template <typename TValue, typename TSpec, typename TSize>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> &,
    TSize,
    Limit)
{
    return 0;
}

template <typename TValue, typename TSpec, typename TSize, typename TSize2>
inline typename Value<String<TValue, TSpec> >::Type *
_reallocateStorage(
    String<TValue, TSpec> &,
    TSize,
    TSize2,
    Limit)
{
    return 0;
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize_>
inline void
_reserveStorage(
    String<TValue, TSpec> & /*seq*/,
    TSize_ /*new_capacity*/,
    Insist)
{
    // do nothing
}

template <typename TValue, typename TSpec, typename TSize_>
inline void
_reserveStorage(
    String<TValue, TSpec> & /*seq*/,
    TSize_ /*new_capacity*/,
    Limit)
{
    // do nothing
}

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline void
_reserveStorage(
    String<TValue, TSpec> & seq,
    TSize_ new_capacity,
    Tag<TExpand> tag)
{
    typedef typename Size< String<TValue, TSpec> >::Type TSize;

    TSize old_capacity = capacity(seq);

//  if (old_capacity == (TSize)new_capacity) return;
//  if (!IsSameType<TExpand,TagExact_>::VALUE && old_capacity > (TSize)new_capacity) return;
    if (old_capacity >= (TSize)new_capacity) return;

    TSize seq_length = length(seq);
    typename Value< String<TValue, TSpec> >::Type * old_array = _reallocateStorage(seq, new_capacity, tag);
    if (old_array)
    {//buffer was replaced, destruct old buffer
//        arrayConstructCopy(old_array, old_array + seq_length, begin(seq, Standard()));
        arrayConstructMove(old_array, old_array + seq_length, begin(seq, Standard()));
        arrayDestruct(old_array, old_array + seq_length);
        _deallocateStorage(seq, old_array, old_capacity);
    }
    _setLength(seq, seq_length);
}

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
reserve(
    String<TValue, TSpec> & seq,
    TSize_ new_capacity,
    Tag<TExpand> tag)
{
    _reserveStorage(seq, new_capacity, tag);
    return _capacityReturned(seq, new_capacity, tag);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TExpand>
struct _Resize_String
{
    template <typename T>
    static inline typename Size<T>::Type
    resize_(
        T & me,
        typename Size<T>::Type new_length)
    {
        typedef typename Size<T>::Type TSize;
        TSize me_length = length(me);
        if (new_length < me_length)
        {
            arrayDestruct(begin(me, Standard()) + new_length, begin(me, Standard()) + me_length);
        }
        else
        {
            typename Size<T>::Type me_capacity = capacity(me);
            if (new_length > me_capacity)
            {
                TSize new_capacity = reserve(me, new_length, TExpand());
                if (new_capacity < new_length)
                {
                    new_length = new_capacity;
                }
            }
            if (new_length > me_length)
            {
                arrayConstruct(begin(me, Standard()) + me_length, begin(me, Standard()) + new_length);
            }
        }

        _setLength(me, new_length);
        return new_length;
    }

    template <typename T, typename TValue>
    static inline typename Size<T>::Type
    resize_(
        T & me,
        typename Size<T>::Type new_length,
        TValue const & val)
    {
        typedef typename Size<T>::Type TSize;
        TSize me_length = length(me);
        if (new_length < me_length)
        {
            arrayDestruct(begin(me, Standard()) + new_length, begin(me, Standard()) + me_length);
        }
        else
        {
            TSize me_capacity = capacity(me);
            if (new_length > me_capacity)
            {
                TValue tempCopy = val;  // reserve could invalidate val
                TSize new_capacity = reserve(me, new_length, TExpand());
                if (new_capacity < new_length)
                {
                    new_length = new_capacity;
                }
                arrayConstruct(begin(me, Standard()) + me_length, begin(me, Standard()) + new_length, tempCopy);
            } else
                if (new_length > me_length)
                {
                    arrayConstruct(begin(me, Standard()) + me_length, begin(me, Standard()) + new_length, val);
                }
        }

        _setLength(me, new_length);
        return new_length;
    }
};

template <typename TValue, typename TSpec, typename TSize, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
resize(
    String<TValue, TSpec> & me,
    TSize new_length,
    Tag<TExpand>)
{
    return _Resize_String<Tag<TExpand> >::resize_(me, new_length);
}

template <typename TValue, typename TSpec, typename TSize, typename TValue2, typename TExpand>
inline typename Size< String<TValue, TSpec> >::Type
resize(
    String<TValue, TSpec> & me,
    TSize new_length,
    TValue2 const & val,
    Tag<TExpand>)
{
    return _Resize_String<Tag<TExpand> >::resize_(me, new_length, val);
}

// ----------------------------------------------------------------------------
// Function operator+=(), shortcut to append()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Maybe move to a string_shortcuts.h?

template <typename TLeftValue, typename TLeftSpec, typename TRight >
inline
String<TLeftValue, TLeftSpec> &
operator+=(String<TLeftValue, TLeftSpec> & left,
           TRight const & right)
{
    append(left, right);
    return left;
}

// The comparator functions forward to the code in sequence_lexical.h.

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TLeftValue, typename TLeftSpec, typename TRight >
inline bool
operator==(String<TLeftValue, TLeftSpec> const & left,
           TRight const & right)
{
    typename Comparator<String<TLeftValue, TLeftSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}

template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator==(TLeftValue * left,
           String<TRightValue, TRightSpec> const & right)
{
    typename Comparator<String<TRightValue, TRightSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TLeftValue, typename TLeftSpec, typename TRight >
inline bool
operator!=(String<TLeftValue, TLeftSpec> const & left,
           TRight const & right)
{
    typename Comparator<String<TLeftValue, TLeftSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}

template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator!=(TLeftValue * left,
           String<TRightValue, TRightSpec> const & right)
{
    typename Comparator<String<TRightValue, TRightSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TLeftValue, typename TLeftSpec, typename TRight>
inline bool
operator<(String<TLeftValue, TLeftSpec> const & left,
          TRight const & right)
{
    return isLess(left, right, typename DefaultPrefixOrder<String<TLeftValue, TLeftSpec> >::Type());
}

template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator<(TLeftValue * left,
          String<TRightValue, TRightSpec> const & right)
{
    return isLess(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TLeftValue, typename TLeftSpec, typename TRight>
inline bool
operator<=(String<TLeftValue, TLeftSpec> const & left,
           TRight const & right)
{
    return isLessOrEqual(left, right, typename DefaultPrefixOrder<String<TLeftValue, TLeftSpec> >::Type());
}

template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator<=(TLeftValue * left,
           String<TRightValue, TRightSpec> const & right)
{
    return isLessOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

// ----------------------------------------------------------------------------
// Function operator>()
// ----------------------------------------------------------------------------

template <typename TLeftValue, typename TLeftSpec, typename TRight>
inline bool
operator>(String<TLeftValue, TLeftSpec> const & left,
          TRight const & right)
{
    return isGreater(left, right, typename DefaultPrefixOrder<String<TLeftValue, TLeftSpec> >::Type());
}
template <typename TLeftValue, typename TRightValue, typename TRightSpec >
inline bool
operator>(TLeftValue * left,
          String<TRightValue, TRightSpec> const & right)
{
    return isGreater(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TLeftValue, typename TLeftSpec, typename TRight>
inline bool
operator>=(String<TLeftValue, TLeftSpec> const & left,
           TRight const & right)
{
    return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<String<TLeftValue, TLeftSpec> >::Type());
}
template <typename TLeftValue, typename TRightValue, typename TRightSpec>
inline bool
operator>=(TLeftValue * left,
             String<TRightValue, TRightSpec> const & right)
{
    return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

// ----------------------------------------------------------------------------
// Function beginPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn String#beginPosition
 * @headerfile <seqan/sequence.h>
 * @brief Return 0 for compatibility with @link Segment @endlink.
 *
 * @signature TPos beginPosition(str);
 *
 * @param[in] seg The String to use.
 *
 * @return TPos Always 0.
 */


// ----------------------------------------------------------------------------
// Function endPosition()
// ----------------------------------------------------------------------------

/*!
 * @fn String#endPosition
 * @headerfile <seqan/sequence.h>
 * @brief Return length of string for compatibility with @link Segment @endlink.
 *
 * @signature TPos endPosition(str);
 *
 * @param[in] seg The string to use.
 *
 * @return TPos Length of the string.
 */

// ----------------------------------------------------------------------------
// Function reserveChunk()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize>
inline void reserveChunk(String<TValue, TSpec> &str, TSize size, Output)
{
    reserve(str, length(str) + size);
}

template <typename TValue, typename TSpec, typename TSize>
inline void reserveChunk(String<TValue, TSpec> const &, TSize, Input)
{}

// ----------------------------------------------------------------------------
// Function getChunk()
// ----------------------------------------------------------------------------

template <typename TChunk, typename TValue, typename TSpec>
inline void
getChunk(TChunk &result, String<TValue, TSpec> &cont, Output)
{
    return assignRange(result, end(cont, Standard()), begin(cont, Standard()) + capacity(cont));
}

// ----------------------------------------------------------------------------
// Function advanceChunk()
// ----------------------------------------------------------------------------

// extend target string size
template <typename TValue, typename TSpec, typename TSize>
inline void advanceChunk(String<TValue, TSpec> &str, TSize size)
{
    _setLength(str, length(str) + size);
}

template <typename TValue, typename TStringSpec, typename TSpec, typename TSize>
inline void advanceChunk(Iter<String<TValue, TStringSpec>, TSpec> &iter, TSize size)
{
    typedef String<TValue, TStringSpec> TContainer;
    typedef Iter<TContainer, TSpec> TIter;

    iter += size;

    TContainer &cont = container(iter);
    typename Position<TIter>::Type pos = position(iter);

    if (pos > length(cont))
        _setLength(cont, pos);
}

// ----------------------------------------------------------------------------
// Function operator<<() for streams.
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator<<(TStream & target,
           String<TValue, TSpec> const & source)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, source);
    return target;
}

// ----------------------------------------------------------------------------
// Function operator>>() for streams.
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TSpec>
inline TStream &
operator>>(TStream & source,
           String<TValue, TSpec> & target)
{
    typename DirectionIterator<TStream, Input>::Type it = directionIterator(source, Input());;
    read(it, target);
    return source;
}

// ----------------------------------------------------------------------------
// Function assignValueById
// ----------------------------------------------------------------------------

template<typename TValue, typename TSpec, typename TId, typename TValue2>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TId> >, void)
assignValueById(String<TValue, TSpec> & me,
                TId id,
                TValue2 const & obj)
{
    if (length(me) <= id)
        resize(me, id + 1, TValue());
    assignValue(me, id, obj);
}

// ----------------------------------------------------------------------------
// Function getValueById
// ----------------------------------------------------------------------------

template<typename TValue, typename TSpec, typename TId>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TId> >, TValue)
getValueById(String<TValue, TSpec> & me,
             TId id)
{
    if (id < length(me))
        return getValue(me, id);
    else
        return TValue();
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_ARRAY_BASE_H_
