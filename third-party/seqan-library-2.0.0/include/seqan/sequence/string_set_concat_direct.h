// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Implementation of ConcatDirect string set, a string set storing the
// concatenation of all strings within one string.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_SET_CONCAT_DIRECT_H_
#define SEQAN_SEQUENCE_STRING_SET_CONCAT_DIRECT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSpec = void>
struct ConcatDirect;                    // contains 1 string (the concatenation of n strings)

// TODO(holtgrew): Change name of specialization to ConcatDirect Owner StringSet?

/*!
 * @class ConcatDirectStringSet ConcatDirect StringSet
 * @extends OwnerStringSet
 * @headerfile <seqan/sequence.h>
 * @brief Owner StringSet implementation that stores strings in one large underlying string.
 *
 * @signature template <typename TString>
 *            class StringSet<TString, Owner<ConcatDirect> >;
 *
 * @tparam TString The type of the string to store in the string set.
 *
 * Storing multiple strings in one larger one with storing the positions between strings leads to a very compact
 * representation with a predictable memory layout.
 *
 * At the moment, ConcatDirect StringSet objects only support appending data.
 *
 * @var TConcatenator ConcatDirectStringSet::concat;
 * @brief The concatenation string.  Concatenates all sequences of the StringSet without gaps.
 */

template <typename TString, typename TSpec>
class StringSet<TString, Owner<ConcatDirect<TSpec> > >
{
public:
    typedef typename StringSetLimits<StringSet>::Type   TLimits;
    typedef typename Concatenator<StringSet>::Type      TConcatenator;

    TLimits         limits;
    TConcatenator   concat;

    StringSet()
    {
        _initStringSetLimits(*this);
    }

    template <typename TOtherString, typename TOtherSpec>
    StringSet(StringSet<TOtherString, TOtherSpec> &other)
    {
        _initStringSetLimits(*this);
        assign(*this, other);
    }

    template <typename TOtherString, typename TOtherSpec>
    StringSet(StringSet<TOtherString, TOtherSpec> const &other)
    {
        _initStringSetLimits(*this);
        assign(*this, other);
    }

    template <typename TOtherSpec>
    StringSet(String<TString, TOtherSpec> const &other)
    {
        _initStringSetLimits(*this);
        assign(*this, other);
    }

    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    template <typename TPos>
    SEQAN_HOST_DEVICE inline typename Reference<StringSet>::Type
    operator[](TPos pos)
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }

    template <typename TPos>
    SEQAN_HOST_DEVICE inline typename Reference<StringSet const>::Type
    operator[](TPos pos) const
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }

    template <typename TStringSet>
    StringSet & operator= (TStringSet const &other)
    {
        assign(*this, other);
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction Concatenator
// --------------------------------------------------------------------------

template <typename TString, typename TSpec >
struct Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >
{
    typedef TString Type;
};

// --------------------------------------------------------------------------
// Metafunction Value
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Value<StringSet<TString, Owner<ConcatDirect<TSpec> > > >
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type> {};

template <typename TString, typename TSpec>
struct Value<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type const> {};

// --------------------------------------------------------------------------
// Metafunction GetValue
// --------------------------------------------------------------------------

/* // we had a problem with constructing a Finder<GetValue<TConcatStringSet>::Type>
template <typename TString, typename TSpec >
struct GetValue<StringSet<TString, Owner<ConcatDirect<TSpec> > > >:
    Infix<TString> {};

template <typename TString, typename TSpec >
struct GetValue<StringSet<TString, Owner<ConcatDirect<TSpec> > > const >:
    Infix<TString const> {};
*/

template <typename TString, typename TSpec >
struct GetValue<StringSet<TString, Owner<ConcatDirect<TSpec> > > >
{
    typedef typename Infix<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type const Type;
};

template <typename TString, typename TSpec >
struct GetValue<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>
{
    typedef typename Infix<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>::Type const Type;
};

// --------------------------------------------------------------------------
// Metafunction Reference
// --------------------------------------------------------------------------

template <typename TString, typename TSpec >
struct Reference<StringSet<TString, Owner<ConcatDirect<TSpec> > > >
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type> {};

template <typename TString, typename TSpec >
struct Reference<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type const> {};

// --------------------------------------------------------------------------
// Metafunction Prefix
// --------------------------------------------------------------------------

template <typename TString, typename TSpec >
struct Prefix<StringSet<TString, Owner<ConcatDirect<TSpec> > > >
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type> {};

template <typename TString, typename TSpec >
struct Prefix<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type const> {};

// --------------------------------------------------------------------------
// Metafunction Suffix
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
struct Suffix<StringSet<TString, Owner<ConcatDirect<TSpec> > > >
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type> {};

template <typename TString, typename TSpec>
struct Suffix<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type const> {};

// --------------------------------------------------------------------------
// Metafunction Infix
// --------------------------------------------------------------------------

template <typename TString, typename TSpec >
struct Infix<StringSet<TString, Owner<ConcatDirect<TSpec> > > >
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type> {};

template <typename TString, typename TSpec >
struct Infix<StringSet<TString, Owner<ConcatDirect<TSpec> > > const >
    : Infix<typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type const> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec>
typename View<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
view(StringSet<TString, Owner<ConcatDirect<TSpec> > > & stringSet)
{
    typename View<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type stringSetView;

    concat(stringSetView) = view(concat(stringSet));
    stringSetLimits(stringSetView) = view(stringSetLimits(stringSet));

    return stringSetView;
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TString2, typename TSpec2>
void assign(StringSet<TString, Owner<ConcatDirect<TSpec> > > & stringSet,
            StringSet<TString2, Owner<ConcatDirect<TSpec2> > > const & other)
{
    assign(concat(stringSet), concat(other));
    assign(stringSetLimits(stringSet), stringSetLimits(other));
}

template <typename TString, typename TSpec, typename TString2, typename TSpec2>
void assign(StringSet<TString, Owner<ConcatDirect<TSpec> > > & stringSet,
            StringSet<TString2, Owner<ConcatDirect<TSpec2> > > & other)
{
    assign(concat(stringSet), concat(other));
    assign(stringSetLimits(stringSet), stringSetLimits(other));
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPos, typename TSequence >
inline void assignValue(
    StringSet<TString, Owner<ConcatDirect<TSpec> > > & me,
    TPos pos,
    TSequence const & seq)
{
    typedef StringSet<TString, Owner<ConcatDirect<TSpec> > > TStringSet;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename StringSetLimits<TStringSet>::Type TLimits;
    typedef typename Value<TLimits>::Type TLimitValue;
    typedef typename MakeSigned<TLimitValue>::Type TSignedLimitValue;

    TSignedLimitValue oldSize = length(me[pos]);
    replace(me.concat, me.limits[pos], me.limits[pos + 1], seq);
    if (_validStringSetLimits(me))
    {
        TSignedLimitValue delta = (TSignedLimitValue)length(seq) - oldSize;
        TSize size = length(me);
        while (pos < size)
            me.limits[++pos] += delta;
    }
}

// --------------------------------------------------------------------------
// Function _initStringSetLimits
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline void _initStringSetLimits(StringSet<TString, TSpec> & me)
{
    appendValue(me.limits, 0);
}

// --------------------------------------------------------------------------
// Function _validStringSetLimits
// --------------------------------------------------------------------------

template <typename TString, typename TSpec >
inline bool _validStringSetLimits(StringSet<TString, Owner<ConcatDirect<TSpec> > > const &)
{
    return true;
}

// --------------------------------------------------------------------------
// Function _refreshStringSetLimits()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec >
inline void _refreshStringSetLimits(StringSet<TString, Owner<ConcatDirect<TSpec> > > &) {}

// --------------------------------------------------------------------------
// Function appendValue()
// --------------------------------------------------------------------------

template <typename TString, typename TString2, typename TSpec, typename TExpand>
inline void appendValue(
    StringSet<TString, Owner<ConcatDirect<TSpec> > > & me,
    TString2 const & obj,
    Tag<TExpand> tag)
{
    appendValue(me.limits, lengthSum(me) + length(obj), tag);
    append(me.concat, obj, tag);
}

// --------------------------------------------------------------------------
// Function insertValue()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPos, typename TSequence, typename TExpand >
inline void insertValue(
    StringSet<TString, Owner<ConcatDirect<TSpec> > > & me,
    TPos pos,
    TSequence const & seq,
    Tag<TExpand> tag)
{
    typedef StringSet<TString, Owner<ConcatDirect<TSpec> > > TStringSet;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename StringSetLimits<TStringSet>::Type TLimits;
    typedef typename Value<TLimits>::Type TLimitValue;

    replace(me.concat, me.limits[pos], me.limits[pos], seq, tag);
    insertValue(me.limits, pos, me.limits[pos], tag);
    TLimitValue delta = (TLimitValue)length(seq);
    TSize size = length(me);
    while (pos <size)
        me.limits[++pos] += delta;
}

// --------------------------------------------------------------------------
// Function erase()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPos >
inline void erase(
    StringSet<TString, Owner<ConcatDirect<TSpec> > > & me,
    TPos pos,
    TPos pos_end)
{
    typedef StringSet<TString, Owner<ConcatDirect<TSpec> > > TStringSet;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename StringSetLimits<TStringSet>::Type TLimits;
    typedef typename Value<TLimits>::Type TLimitValue;

    erase(me.concat, me.limits[pos], me.limits[pos_end]);

    TLimitValue lengthSum = 0;
    for (TSize i = pos; i < pos_end; ++i)
        lengthSum += me.limits[i];

    erase(me.limits, pos, pos_end);

    TSize size = length(me);
    while (pos <size)
        me.limits[++pos] -= lengthSum;
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec >
inline void clear(StringSet<TString, Owner<ConcatDirect<TSpec> > > & me)
{
    clear(me.concat);
    resize(me.limits, 1, Exact());
}

// --------------------------------------------------------------------------
// Function length()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
SEQAN_HOST_DEVICE inline typename Size<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
length(StringSet<TString, Owner<ConcatDirect<TSpec> > > const & me)
{
    return length(me.limits) - 1;
}

// --------------------------------------------------------------------------
// Function resize()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TSize, typename TExpand >
inline typename Size<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
resize(StringSet<TString, Owner<ConcatDirect<TSpec> > > & me, TSize new_size, Tag<TExpand> tag)
{
    if (new_size < length(me.limits))
    {
        resize(me.concat, me.limits[new_size]);
        return resize(me.limits, new_size + 1, tag) - 1;
    } else
        return resize(me.limits, new_size + 1, back(me.limits), tag) - 1;
}


// --------------------------------------------------------------------------
// Function reserve()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TSize, typename TExpand>
inline typename Size<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
reserve(StringSet<TString, Owner<ConcatDirect<TSpec> > > & me,
        TSize const & new_capacity,
        Tag<TExpand> tag)
{
    return reserve(me.limits, new_capacity + 1, tag) - 1;
}

// --------------------------------------------------------------------------
// Function prefix()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPosition >
inline typename Prefix<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
prefix(StringSet<TString, Owner<ConcatDirect<TSpec> > > & me, TPosition pos)
{
    return infix(me.concat, stringSetLimits(me)[getSeqNo(pos, stringSetLimits(me))], posGlobalize(pos, stringSetLimits(me)));
}

template <typename TString, typename TSpec, typename TPosition >
inline typename Prefix<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>::Type
prefix(StringSet<TString, Owner<ConcatDirect<TSpec> > > const & me, TPosition pos)
{
    return infix(me.concat, stringSetLimits(me)[getSeqNo(pos, stringSetLimits(me))], posGlobalize(pos, stringSetLimits(me)));
}

// --------------------------------------------------------------------------
// Function suffix()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPosition >
inline typename Suffix<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
suffix(StringSet<TString, Owner<ConcatDirect<TSpec> > > & me, TPosition pos)
{
    return infix(me.concat, posGlobalize(pos, stringSetLimits(me)), stringSetLimits(me)[getSeqNo(pos, stringSetLimits(me)) + 1]);
}

template <typename TString, typename TSpec, typename TPosition >
inline typename Suffix<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>::Type
suffix(StringSet<TString, Owner<ConcatDirect<TSpec> > > const & me, TPosition pos)
{
    return infix(me.concat, posGlobalize(pos, stringSetLimits(me)), stringSetLimits(me)[getSeqNo(pos, stringSetLimits(me)) + 1]);
}

// --------------------------------------------------------------------------
// Function infix()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPosBegin, typename TPosEnd >
inline typename Infix<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
infix(StringSet<TString, Owner<ConcatDirect<TSpec> > > & me, TPosBegin posBegin, TPosEnd posEnd)
{
    return infix(me.concat, posGlobalize(posBegin, stringSetLimits(me)), posGlobalize(posEnd, stringSetLimits(me)));
}

template <typename TString, typename TSpec, typename TPosBegin, typename TPosEnd >
inline typename Infix<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>::Type
infix(StringSet<TString, Owner<ConcatDirect<TSpec> > > const & me, TPosBegin posBegin, TPosEnd posEnd)
{
    return infix(me.concat, posGlobalize(posBegin, stringSetLimits(me)), posGlobalize(posEnd, stringSetLimits(me)));
}

// --------------------------------------------------------------------------
// Function infixWithLength()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPosition, typename TSize >
inline typename Infix<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
infixWithLength(StringSet<TString, Owner<ConcatDirect<TSpec> > > & me, TPosition pos, TSize length)
{
    return infixWithLength(me.concat, posGlobalize(pos, stringSetLimits(me)), length);
}

template <typename TString, typename TSpec, typename TPosition, typename TSize >
inline typename Infix<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>::Type
infixWithLength(StringSet<TString, Owner<ConcatDirect<TSpec> > > const & me, TPosition pos, TSize length)
{
    return infixWithLength(me.concat, posGlobalize(pos, stringSetLimits(me)), length);
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPos >
SEQAN_HOST_DEVICE inline typename Value<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type
value(StringSet<TString, Owner<ConcatDirect<TSpec> > > & me, TPos pos)
{
    return infix(me.concat, me.limits[pos], me.limits[pos + 1]);
}

template <typename TString, typename TSpec, typename TPos >
SEQAN_HOST_DEVICE inline typename Value<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>::Type
value(StringSet<TString, Owner<ConcatDirect<TSpec> > > const & me, TPos pos)
{
    return infix(me.concat, me.limits[pos], me.limits[pos + 1]);
}

// --------------------------------------------------------------------------
// Function concat()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
inline typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > >::Type &
concat(StringSet<TString, Owner<ConcatDirect<TSpec> > > & me)
{
    return me.concat;
}

template <typename TString, typename TSpec>
inline typename Concatenator<StringSet<TString, Owner<ConcatDirect<TSpec> > > const>::Type &
concat(StringSet<TString, Owner<ConcatDirect<TSpec> > > const & me)
{
    return me.concat;
}

// --------------------------------------------------------------------------
// Function swap()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
void swap(StringSet<TString, Owner<ConcatDirect<TSpec> > > & lhs,
          StringSet<TString, Owner<ConcatDirect<TSpec> > > & rhs)
{
    swap(lhs.limits, rhs.limits);
    swap(lhs.concat, rhs.concat);
}

// ============================================================================
// Device Functions
// ============================================================================
// NOTE(esiragusa): All functions are equivalent to the originals - overloaded to remove SEQAN_HOST_DEVICE :(

#ifdef PLATFORM_CUDA

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec, typename TPos >
inline typename Value<StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > >::Type
value(StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > & me, TPos pos)
{
    return infix(me.concat, me.limits[pos], me.limits[pos + 1]);
}

template <typename TValue, typename TAlloc, typename TSpec, typename TPos >
inline typename Value<StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > const>::Type
value(StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > const & me, TPos pos)
{
    return infix(me.concat, me.limits[pos], me.limits[pos + 1]);
}

// --------------------------------------------------------------------------
// Function length()
// --------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
inline typename Size<StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > >::Type
length(StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > const & me)
{
    return length(me.limits) - 1;
}

// --------------------------------------------------------------------------
// Function back()
// --------------------------------------------------------------------------

template <typename TValue, typename TAlloc, typename TSpec>
inline typename Reference<StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > const>::Type
back(StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > const & me)
{
    return value(me, length(me) - 1);
}

template <typename TValue, typename TAlloc, typename TSpec>
inline typename Reference<StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > >::Type
back(StringSet<thrust::device_vector<TValue, TAlloc>, Owner<ConcatDirect<TSpec> > > & me)
{
    return value(me, length(me) - 1);
}

#endif // PLATFORM_CUDA

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_SET_CONCAT_DIRECT_H_
