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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
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
 * @signature template <typename TValue, size_t LENGTH>
 *            class String<TValue, Array<LENGTH> >;
 * 
 * @tparam TValue The value type, that is the type of the items/characters
 *                stored in the string.Use @link Value @endlink to get the value
 *                type for a given class.
 * @tparam LENGTH A positive integer that specifies the capacity of the
 *                string.Note that the capacity of an Array String is fixed at
 *                compile-time.
 * 
 * Having static-sized strings is useful as members of structs for external memory algorithms, for example.
 */

/**
.Spec.Array String:
..cat:Strings
..general:Class.String
..summary:Fast, static-size string.
..remarks:This is useful as members of structs for external memory algorithms, for example.
..signature:String<TValue, Array<LENGTH> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.LENGTH:A positive integer that specifies the capacity of the string.
...remarks:Note that the capacity of an Array String is fixed at compile-time.
..include:seqan/sequence.h
*/

template <size_t LENGTH>
struct Array;

template <typename TValue, size_t LENGTH>
class String<TValue, Array<LENGTH> >
{
public:
    mutable TValue   data_begin[LENGTH];    // TODO(weese): mutable seems to be a workaround for non-const iterators returned by
    TValue * data_end;                      //              Iterator<> for const Alloc/Array strings. Must be fixed to work without mutable.

    String()
    {
        SEQAN_CHECKPOINT;
        data_end = data_begin;
    }

    template <typename TSource>
    String(TSource & source)
    {
        SEQAN_CHECKPOINT;
        data_end = data_begin;
        assign(*this, source);
    }

    template <typename TSource>
    String(TSource const & source)
    {
        SEQAN_CHECKPOINT;
        data_end = data_begin;
        assign(*this, source);
    }

    String(String const & source)
    {
        SEQAN_CHECKPOINT;
        data_end = data_begin;
        assign(*this, source);
    }

    template <typename TSource>
    inline
    String & operator=(TSource const & source)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source);
        return *this;
    }

    String & operator=(String const & source)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source);
        return *this;
    }

/*  // TODO(weese): Empty destructors are superfluos
    ~String()
    {
    }
*/
    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    template <typename TPos>
    inline typename Reference<String>::Type
    operator[](TPos pos)
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<String const>::Type
    operator[](TPos pos) const
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// ----------------------------------------------------------------------------

template <typename TValue, size_t LENGTH>
struct DefaultOverflowImplicit<String<TValue, Array<LENGTH> > >
{
    typedef Limit Type;
};

template <typename TValue, size_t LENGTH>
struct DefaultOverflowImplicit<String<TValue, Array<LENGTH> > const>
{
    typedef Limit Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultOverflowExplicit
// ----------------------------------------------------------------------------

template <typename TValue, size_t LENGTH>
struct DefaultOverflowExplicit<String<TValue, Array<LENGTH> > >
{
    typedef Limit Type;
};

template <typename TValue, size_t LENGTH>
struct DefaultOverflowExplicit<String<TValue, Array<LENGTH> > const>
{
    typedef Limit Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsContiguous
// ----------------------------------------------------------------------------

template <typename TValue, size_t LENGTH>
struct IsContiguous<String<TValue, Array<LENGTH> > >
{
    typedef True Type;
    enum { VALUE = true };
};

// ----------------------------------------------------------------------------
// Metafunction LENGTH
// ----------------------------------------------------------------------------

///.Metafunction.LENGTH.param.T.type:Spec.Array String
template <typename TValue, size_t LENGTH_>
struct LENGTH<String<TValue, Array<LENGTH_> > >
{
    enum { VALUE = LENGTH_ };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TValue, size_t LENGTH>
inline typename Iterator<String<TValue, Array<LENGTH> >, Standard>::Type
begin(String<TValue, Array<LENGTH> > & me,
      Standard const &)
{
    SEQAN_CHECKPOINT;
    return me.data_begin;
}
template <typename TValue, size_t LENGTH>
inline typename Iterator<String<TValue, Array<LENGTH> > const, Standard>::Type
begin(String<TValue, Array<LENGTH> > const & me,
      Standard const & )
{
    SEQAN_CHECKPOINT;
    return me.data_begin;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TValue, size_t LENGTH>
inline typename Iterator<String<TValue, Array<LENGTH> >, Standard>::Type
end(String<TValue, Array<LENGTH> > & me,
    Standard const &)
{
    SEQAN_CHECKPOINT;
    return me.data_end;
}
template <typename TValue, size_t LENGTH>
inline typename Iterator<String<TValue, Array<LENGTH> > const, Standard>::Type
end(String<TValue, Array<LENGTH> > const & me,
    Standard const &)
{
    SEQAN_CHECKPOINT;
    return me.data_end;
}

// ----------------------------------------------------------------------------
// Function capacity()
// ----------------------------------------------------------------------------

template <typename TValue, size_t LENGTH>
inline typename Size<String<TValue, Array<LENGTH> > >::Type
capacity(String<TValue, Array<LENGTH> > const &)
{
    SEQAN_CHECKPOINT;
    return LENGTH;
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TValue, size_t LENGTH, typename TSize, typename TExpand>
inline typename Size<String<TValue, Array<LENGTH> > >::Type
reserve(String<TValue, Array<LENGTH> > & me,
        TSize,
        Tag<TExpand>)
{
SEQAN_CHECKPOINT
    return capacity(me);
}


// ----------------------------------------------------------------------------
// Function _setLength()
// ----------------------------------------------------------------------------

/**
.Internal._setLength.param.object.type:Spec.Array String
*/
template <typename TValue, size_t LENGTH, typename TSize>
inline void
_setLength(String<TValue, Array<LENGTH> > & me,
           TSize new_length)
{
    SEQAN_CHECKPOINT;
    me.data_end = me.data_begin + new_length;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_ARRAY_H_
