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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Implementation of the StringSet specialization Dependent<Generous>.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_SET_DEPENDENT_GENEROUS_H_
#define SEQAN_SEQUENCE_STRING_SET_DEPENDENT_GENEROUSH_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class GenerousDependentStringSet Generous Dependent StringSet
 * @extends DependentStringSet
 * @headerfile <seqan/sequence.h>
 * @brief Dependent StringSet implementation with efficient sequence access by sequence id.
 *
 * @signature template <typename TString>
 *            class StringSet<TString, Depedent<Generous> >;
 *
 * @tparam TString The type of the string to store in the string set.
 *
 * See @link TightDependentStringSet @endlink for a Dependent StringSet implementation with a more memory efficient
 * representation at higher costs for by-id element access.
 */

// The Dddoc documentation of the Dependent specialization is in string_set_tight.h
// since Tight is the default specialization of Dependent.

// StringSet with individual sequences in a string of string pointers
template <typename TString>
class StringSet<TString, Dependent<Generous> >
{
public:
    typedef String<TString *>                           TStrings;
    typedef typename Size<StringSet>::Type              TSize;
    typedef typename StringSetLimits<StringSet>::Type   TLimits;
    typedef typename Concatenator<StringSet>::Type      TConcatenator;

    TStrings        strings;
    TLimits         limits;
    bool            limitsValid;        // is true if limits contains the cumulative sum of the sequence lengths
    TConcatenator   concat;

    StringSet()
        : limitsValid(true)
    {
        SEQAN_CHECKPOINT;
        appendValue(limits, 0);
    }

    template <typename TDefault>
    StringSet(StringSet<TString, Owner<TDefault> > const& _other)
        : limitsValid(true)
    {
        SEQAN_CHECKPOINT;
        appendValue(limits, 0);
        for (unsigned int i = 0; i < length(_other); ++i)
            appendValue(*this, _other[i]);
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

template <typename TString, typename TExpand >
inline void appendValue(
    StringSet<TString, Dependent<Generous> > & me,
    TString const & obj,
    Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    appendValue(me.limits, lengthSum(me) + length(obj), tag);
    appendValue(me.strings, const_cast<TString*>(&obj), tag);
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename TString >
inline void clear(StringSet<TString, Dependent<Generous> > & me)
{
    SEQAN_CHECKPOINT;
    clear(me.strings);
    resize(me.limits, 1, Exact());
    me.limitsValid = true;
}

// --------------------------------------------------------------------------
// Function length()
// --------------------------------------------------------------------------

template <typename TString >
inline typename Size<StringSet<TString, Dependent<Generous> > >::Type
length(StringSet<TString, Dependent<Generous> > & me)
{
    return length(me.limits) - 1;
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TString, typename TPos >
inline typename Reference<StringSet<TString, Dependent<Generous> > >::Type
value(StringSet<TString, Dependent<Generous> >& me, TPos pos)
{
    SEQAN_CHECKPOINT;
    unsigned i = _findIthNonZeroValue(me.strings, pos);
    if (i <length(me.strings))
        return *me.strings[i];
    static TString tmp = "";
    return tmp;
}

template <typename TString, typename TPos >
inline typename Reference<StringSet<TString, Dependent<Generous> > const >::Type
value(StringSet<TString, Dependent<Generous> > const & me, TPos pos)
{
    SEQAN_CHECKPOINT;
    unsigned i = _findIthNonZeroValue(me.strings, pos);
    if (i < length(me.strings))
        return *me.strings[i];
    static TString tmp = "";
    return tmp;
}

// --------------------------------------------------------------------------
// Function getValueById()
// --------------------------------------------------------------------------

template <typename TString, typename TId>
inline typename Reference<StringSet<TString, Dependent<Generous> > >::Type
getValueById(StringSet<TString, Dependent<Generous> >& me,
            TId const id)
{
    SEQAN_CHECKPOINT;
    if (me.strings[id])
        return *me.strings[id];
    static TString tmp = "";
    return tmp;
}

// --------------------------------------------------------------------------
// Function assignValueById()
// --------------------------------------------------------------------------

template<typename TString, typename TId>
inline typename Id<StringSet<TString, Dependent<Generous> > >::Type
assignValueById(StringSet<TString, Dependent<Generous> >& me,
                TString& obj,
                TId id)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_EQ(length(stringSetLimits(me)), length(me) + 1);
    if (id >= (TId) length(me.strings)) resize(me.strings, id+1, (TString*) 0);
    if ((TString*) me.strings[id] == (TString*) 0)
        resize(me.limits, length(me.limits) + 1, Generous());
    me.strings[id] = &obj;
    me.limitsValid = false;
    SEQAN_ASSERT_EQ(length(stringSetLimits(me)), length(me) + 1);
    return id;
}

// --------------------------------------------------------------------------
// Function removeValueById()
// --------------------------------------------------------------------------

template<typename TString, typename TId>
inline void
removeValueById(StringSet<TString, Dependent<Generous> >& me, TId const id)
{
    SEQAN_CHECKPOINT;
    if (me.strings[id] != (TString*) 0)
    {
        resize(me.limits, length(me.limits) - 1, Generous());
        me.limitsValid = empty(me);
    }
    me.strings[id] = 0;
    while (!empty(me.strings) && !me.strings[length(me.strings) - 1])
        resize(me.strings, length(me.strings) - 1, Generous());
}

// --------------------------------------------------------------------------
// Function positionToId()
// --------------------------------------------------------------------------

template <typename TString, typename TPos>
inline typename Id<StringSet<TString, Dependent<Generous> > >::Type
positionToId(StringSet<TString, Dependent<Generous> >& me,
        TPos const pos)
{
    SEQAN_CHECKPOINT;
    return _findIthNonZeroValue(me.strings,pos);
}

template <typename TString, typename TPos>
inline typename Id<StringSet<TString, Dependent<Generous> > >::Type
positionToId(StringSet<TString, Dependent<Generous> > const& me,
            TPos const pos)
{
    SEQAN_CHECKPOINT;
    return _findIthNonZeroValue(me.strings,pos);
}

// --------------------------------------------------------------------------
// Function idToPosition()
// --------------------------------------------------------------------------

template <typename TString, typename TId>
inline typename Id<StringSet<TString, Dependent<Generous> > >::Type
idToPosition(StringSet<TString, Dependent<Generous> > const& me,
            TId const id)
{
    SEQAN_CHECKPOINT;
    return _countNonZeroValues(me.strings,id);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_SET_DEPENDENT_GENEROUSH_
