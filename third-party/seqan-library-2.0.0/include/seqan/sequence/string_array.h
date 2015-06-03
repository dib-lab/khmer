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
//         David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Implementation of constant-sized Array String class.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_ARRAY_H_
#define SEQAN_SEQUENCE_STRING_ARRAY_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class ArrayString Array String
 * @extends String
 * @headerfile <seqan/sequence.h>
 * @brief Fast, static-size string.
 *
 * @signature template <typename TValue, size_t CAPACITY>
 *            class String<TValue, Array<CAPACITY> >;
 *
 * @tparam TValue The value type, that is the type of the items/characters
 *                stored in the string.Use @link Value @endlink to get the value
 *                type for a given class.
 * @tparam CAPACITY A positive integer that specifies the capacity of the
 *                string.Note that the capacity of an Array String is fixed at
 *                compile-time.
 *
 * Having static-sized strings is useful as members of structs for external memory algorithms, for example.
 */

template <size_t CAPACITY>
struct Array;

// to save memory, we compute the smallest size type required to
// represent values 0..CAPACITY, i.e. CAPACITY+1 different values

template <typename TValue, size_t CAPACITY>
struct Size<String<TValue, Array<CAPACITY> > >:
    BitVector_<Log2<CAPACITY+1>::VALUE>
{};

template <typename TValue, size_t CAPACITY>
class String<TValue, Array<CAPACITY> >
{
public:
    typedef typename Size<String>::Type TSize;

    TSize data_length;
    TValue data_begin[CAPACITY];

    String():
        data_length(0u)
    {}

    template <typename TSource>
    String(TSource & source):
        data_length(0u)
    {
        assign(*this, source);
    }

    template <typename TSource>
    String(TSource const & source):
        data_length(0u)
    {
        assign(*this, source);
    }

    template <typename TSource>
    String & operator=(TSource const & source)
    {
        assign(*this, source);
        return *this;
    }

    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def. (auto inlined)
    // ----------------------------------------------------------------------

    template <typename TPos>
    typename Reference<String>::Type
    operator[](TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    typename Reference<String const>::Type
    operator[](TPos pos) const
    {
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// ----------------------------------------------------------------------------

template <typename TValue, size_t CAPACITY>
struct DefaultOverflowImplicit<String<TValue, Array<CAPACITY> > >
{
    typedef Limit Type;
};

template <typename TValue, size_t CAPACITY>
struct DefaultOverflowImplicit<String<TValue, Array<CAPACITY> > const>
{
    typedef Limit Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultOverflowExplicit
// ----------------------------------------------------------------------------

template <typename TValue, size_t CAPACITY>
struct DefaultOverflowExplicit<String<TValue, Array<CAPACITY> > >
{
    typedef Limit Type;
};

template <typename TValue, size_t CAPACITY>
struct DefaultOverflowExplicit<String<TValue, Array<CAPACITY> > const>
{
    typedef Limit Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsContiguous
// ----------------------------------------------------------------------------

template <typename TValue, size_t CAPACITY>
struct IsContiguous<String<TValue, Array<CAPACITY> > >
{
    typedef True Type;
    enum { VALUE = true };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TValue, size_t CAPACITY>
inline typename Iterator<String<TValue, Array<CAPACITY> >, Standard>::Type
begin(String<TValue, Array<CAPACITY> > & me,
      Standard const &)
{
    return me.data_begin;
}
template <typename TValue, size_t CAPACITY>
inline typename Iterator<String<TValue, Array<CAPACITY> > const, Standard>::Type
begin(String<TValue, Array<CAPACITY> > const & me,
      Standard const & )
{
    return me.data_begin;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TValue, size_t CAPACITY>
inline typename Iterator<String<TValue, Array<CAPACITY> >, Standard>::Type
end(String<TValue, Array<CAPACITY> > & me,
    Standard const &)
{
    return me.data_begin + me.data_length;
}
template <typename TValue, size_t CAPACITY>
inline typename Iterator<String<TValue, Array<CAPACITY> > const, Standard>::Type
end(String<TValue, Array<CAPACITY> > const & me,
    Standard const &)
{
    return me.data_begin + me.data_length;
}

// ----------------------------------------------------------------------------
// Function capacity()
// ----------------------------------------------------------------------------

template <typename TValue, size_t CAPACITY>
inline typename Size<String<TValue, Array<CAPACITY> > >::Type
capacity(String<TValue, Array<CAPACITY> > const &)
{
    return CAPACITY;
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TValue, size_t CAPACITY, typename TSize, typename TExpand>
inline typename Size<String<TValue, Array<CAPACITY> > >::Type
reserve(String<TValue, Array<CAPACITY> > & me,
        TSize,
        Tag<TExpand>)
{
    return capacity(me);
}


// ----------------------------------------------------------------------------
// Function _setLength()
// ----------------------------------------------------------------------------

template <typename TValue, size_t CAPACITY, typename TSize>
inline void
_setLength(String<TValue, Array<CAPACITY> > & me,
           TSize new_length)
{
    SEQAN_ASSERT_LEQ_MSG(static_cast<size_t>(new_length), CAPACITY, "New length would exceed Array String's capacity!");
    me.data_length = new_length;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_ARRAY_H_
