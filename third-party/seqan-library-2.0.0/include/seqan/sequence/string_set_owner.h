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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Implementation of the StringSet specialization Owner.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_SET_OWNER_H_
#define SEQAN_SEQUENCE_STRING_SET_OWNER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// The tag Owner is defined in string_set_base since it is the default
// specialization.

// TODO(holtgrew): Shouldn't Default be void?

// template <typename TSpec = Default>
// struct Owner;

/*!
 * @class OwnerStringSet Owner StringSet
 * @extends StringSet
 * @headerfile <seqan/sequence.h>
 * @brief String set implementation that owns the string.
 *
 * @signature template <typename TString>
 *            class StringSet<TString, Owner<Default> >;
 *
 * @tparam TString The type of the string to store in the string set.
 */

// TODO(holtgrew): Change name of specialization to Owner StringSet.

template <typename TString>
class StringSet<TString, Owner<Default> >
{
public:

    typedef String<TString>                             TStrings;
    typedef typename StringSetLimits<StringSet>::Type   TLimits;
    typedef typename Concatenator<StringSet>::Type      TConcatenator;

    TStrings        strings;
    TLimits         limits;
    bool            limitsValid;        // is true if limits contains the cumulative sum of the sequence lengths
    TConcatenator   concat;

    StringSet() :
        limitsValid(true)
    {
        _initStringSetLimits(*this);
    }

    template <typename TOtherString, typename TOtherSpec>
    StringSet(StringSet<TOtherString, TOtherSpec> const &other) :
        limitsValid(true)
    {
        _initStringSetLimits(*this);
        assign(*this, other);
    }

    template <typename TOtherSpec>
    StringSet(String<TString, TOtherSpec> const &other) :
        limitsValid(true)
    {
        _initStringSetLimits(*this);
        assign(*this, other);
    }

    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    template <typename TPos>
    inline typename Reference<StringSet>::Type
    operator[] (TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<StringSet const>::Type
    operator[] (TPos pos) const
    {
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

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function appendValue()
// --------------------------------------------------------------------------

template <typename TString, typename TString2, typename TExpand >
inline void appendValue(
    StringSet<TString, Owner<Default> > & me,
    TString2 const & obj,
    Tag<TExpand> tag)
{
    // we rather invalidate limits here to allow to do modify appended strings:
    // appendValue(back(stringSet), 'A');

//    if (_validStringSetLimits(me))
//        appendValue(me.limits, lengthSum(me) + length(obj), tag);
    me.limitsValid = false;
    appendValue(me.strings, obj, tag);
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPos, typename TSequence >
inline void assignValue(
    StringSet<TString, Owner<TSpec> > & me,
    TPos pos,
    TSequence const & seq)
{
    typedef StringSet<TString, Owner<TSpec> > TStringSet;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename StringSetLimits<TStringSet>::Type TLimits;
    typedef typename Value<TLimits>::Type TLimitValue;
    typedef typename MakeSigned<TLimitValue>::Type TSignedLimitValue;

    TSignedLimitValue oldSize = length(me[pos]);
    assign(me[pos], seq);
    if (_validStringSetLimits(me))
    {
        TSignedLimitValue delta = (TSignedLimitValue)length(seq) - oldSize;
        TSize size = length(me);
        while (pos < size)
            me.limits[++pos] += delta;
    }
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename TString >
inline void clear(StringSet<TString, Owner<Default> > & me)
{
    SEQAN_CHECKPOINT;
    clear(me.strings);
    resize(me.limits, 1, Exact());
    me.limitsValid = true;
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TString, typename TPos >
inline typename Reference<StringSet<TString, Owner<Default> > >::Type
value(StringSet<TString, Owner<Default> > & me, TPos pos)
{
    return me.strings[pos];
}

template <typename TString, typename TPos >
inline typename Reference<StringSet<TString, Owner<Default> > const >::Type
value(StringSet<TString, Owner<Default> > const & me, TPos pos)
{
    return me.strings[pos];
}

// --------------------------------------------------------------------------
// Function erase()
// --------------------------------------------------------------------------

template <typename TString, typename TPos>
inline typename Size<StringSet<TString, Owner<Default> > >::Type
erase(StringSet<TString, Owner<Default> > & me, TPos pos)
{
    erase(me.strings, pos);
    me.limitsValid = false;
    return length(me);
}

template <typename TString, typename TPos, typename TPosEnd>
inline typename Size<StringSet<TString, Owner<Default> > >::Type
erase(StringSet<TString, Owner<Default> > & me, TPos pos, TPosEnd posEnd)
{
    erase(me.strings, pos, posEnd);
    me.limitsValid = false;
    return length(me);
}

// --------------------------------------------------------------------------
// Function getValueById()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TId>
inline typename Reference<StringSet<TString, Owner<TSpec> > >::Type
getValueById(StringSet<TString, Owner<TSpec> >& me,
            TId const id)
{
    SEQAN_CHECKPOINT;
    if (id < (TId) length(me)) return value(me, id);
    static TString tmp = TString();
    return tmp;
}

template <typename TString, typename TSpec, typename TId>
inline typename Reference<StringSet<TString, Owner<TSpec> > const>::Type
getValueById(StringSet<TString, Owner<TSpec> > const & me,
            TId const id)
{
    if (id < (TId) length(me)) return value(me, id);
    static TString tmp = TString();
    return tmp;
}

// --------------------------------------------------------------------------
// Function assignValueById()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TId>
inline typename Id<StringSet<TString, Owner<TSpec> > >::Type
assignValueById(StringSet<TString, Owner<TSpec> > & me,
                TString& obj,
                TId id)
{
    SEQAN_CHECKPOINT;
    if (id >= (TId) length(me.strings))
    {
        resize(me.strings, id+1, TString());
        resize(me.limits, length(me.limits) + 1, Generous());
    }
    assignValue(me, id, obj);
    me.limitsValid = false;
    return id;
}

// --------------------------------------------------------------------------
// Function removeValueById()
// --------------------------------------------------------------------------

template<typename TString, typename TSpec, typename TId>
inline void
removeValueById(StringSet<TString, Owner<TSpec> > & me, TId const id)
{
    SEQAN_CHECKPOINT;
    erase(me.strings, id);
    resize(me.limits, length(me.limits) - 1, Generous());
    me.limitsValid = empty(me);
}

// --------------------------------------------------------------------------
// Function positionToId()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPos>
inline typename Id<StringSet<TString, Owner<TSpec> > >::Type
positionToId(StringSet<TString, Owner<TSpec> > &,
            TPos const pos)
{
    SEQAN_CHECKPOINT;
    return pos;
}

// --------------------------------------------------------------------------
// Function positionToId()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TPos>
inline typename Id<StringSet<TString, Owner<TSpec> > >::Type
positionToId(StringSet<TString, Owner<TSpec> > const &,
            TPos const pos)
{
    SEQAN_CHECKPOINT;
    return pos;
}

// --------------------------------------------------------------------------
// Function idToPosition()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec, typename TId>
inline typename Id<StringSet<TString, Owner<TSpec> > >::Type
idToPosition(StringSet<TString, Owner<TSpec> > const&,
            TId const id)
{
    SEQAN_CHECKPOINT;
    return id;
}

// --------------------------------------------------------------------------
// Function swap()
// --------------------------------------------------------------------------

template <typename TString, typename TSpec>
void swap(StringSet<TString, Owner<TSpec> > & lhs,
          StringSet<TString, Owner<TSpec> > & rhs)
{
    using std::swap;

    swap(lhs.strings, rhs.strings);
    swap(lhs.limits, rhs.limits);
    swap(lhs.limitsValid, rhs.limitsValid);
    swap(lhs.concat, rhs.concat);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_SET_OWNER_H_
