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
// Implementation of the StringSet specialization Dependent<Tight>, the
// the default specialization of Dependent<>.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_SET_DEPENDENT_TIGHT_H_
#define SEQAN_SEQUENCE_STRING_SET_DEPENDENT_TIGHT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Change name of specialization to Dependent StringSet.

/*!
 * @class DependentStringSet Dependent StringSet
 * @extends StringSet
 * @headerfile <seqan/sequence.h>
 * @brief StringSet implementation that only stores pointers to strings in other string sets.
 *
 * @signature template <typename TString, typename TSpec>
 *            class StringSet<TString, Depedent<TSpec> >;
 *
 * @tparam TString The type of the string to store in the string set.
 * @tparam TSpec   Tag for further specializing the string set.
 *
 * The class is not usable itself, only its subclasses @link TightDependentStringSet @endlink and
 * @link GenerousDependentStringSet @endlink are.
 */

/*!
 * @class TightDependentStringSet Tight Dependent StringSet
 * @extends DependentStringSet
 * @headerfile <seqan/sequence.h>
 * @brief Very space efficient Dependent StringSet implementation.
 *
 * @signature template <typename TString>
 *            class StringSet<TString, Depedent<Tight> >;
 *
 * @tparam TString The type of the string to store in the string set.
 *
 * See @link GenerousDependentStringSet @endlink for a Dependent StringSet implementation that allows for more
 * efficient access to strings in the container via ids at the cost of higher memory usage.
 */

/**
.Spec.Dependent:
..summary:A string set storing references of the strings.
..cat:Sequences
..general:Class.StringSet
..signature:StringSet<TString, Dependent<TSpec> >
..param.TString:The string type.
...type:Class.String
..param.TSpec:The specializing type for the dependent string set.
...default:$Tight$
...remarks:Possible values are $Tight$ or $Generous$
...remarks:$Tight$ is very space efficient whereas $Generous$ provides fast access to the strings in the container via ids.
..include:seqan/sequence.h
 */
// Default id holder string set
template <typename TSpec = Tight>
struct Dependent;

// StringSet with individual sequences in a tight string of string pointers and corr. IDs
template <typename TString>
class StringSet<TString, Dependent<Tight> >
{
public:
    typedef String<TString *>                           TStrings;
    typedef typename Id<StringSet>::Type                TIdType;
    typedef typename Position<StringSet>::Type          TPosition;
    typedef String<TIdType>                             TIds;
    typedef std::map<TIdType, TPosition>                TIdPosMap;
    typedef typename StringSetLimits<StringSet>::Type   TLimits;
    typedef typename Concatenator<StringSet>::Type      TConcatenator;

    TIdType         lastId;
    TStrings        strings;
    TIds            ids;
    TIdPosMap       id_pos_map;
    TLimits         limits;
    bool            limitsValid;        // is true if limits contains the cumulative sum of the sequence lengths
    TConcatenator   concat;

    StringSet()
        : lastId(0), limitsValid(true)
    {
        SEQAN_CHECKPOINT;
        appendValue(limits, 0);
    }

    template <typename TDefault>
    StringSet(StringSet<TString, Owner<TDefault> > const & _other)
        : lastId(0), limitsValid(true)
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
        SEQAN_CHECKPOINT;
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
    StringSet<TString, Dependent<Tight> > & me,
    TString const & obj,
    Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    typedef typename Position<StringSet<TString, Dependent<Tight> > >::Type TPos;
    appendValue(me.limits, lengthSum(me) + length(obj), tag);
    typedef typename StringSet<TString, Dependent<Tight> >::TIdType TIdType;
    appendValue(me.strings, const_cast<TString*>(&obj));
    TIdType last = me.lastId++;
    appendValue(me.ids, last, tag);
    me.id_pos_map.insert(std::make_pair(last, (TPos)(length(me.strings) - 1)));
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename TString >
inline void clear(StringSet<TString, Dependent<Tight> >& me)
{
    SEQAN_CHECKPOINT;
    clear(me.strings);
    me.id_pos_map.clear();
    resize(me.limits, 1, Exact());
    me.limitsValid = true;

    clear(me.ids);
    me.lastId = 0;
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TString, typename TPos >
inline typename Reference<StringSet<TString, Dependent<Tight> > >::Type
value(StringSet<TString, Dependent<Tight> >& me, TPos pos)
{
    SEQAN_CHECKPOINT;
    return *me.strings[pos];
}

template <typename TString, typename TPos >
inline typename Reference<StringSet<TString, Dependent<Tight> > const >::Type
value(StringSet<TString, Dependent<Tight> >const & me, TPos pos)
{
    return *me.strings[pos];
}

// --------------------------------------------------------------------------
// Function getValueById()
// --------------------------------------------------------------------------

template <typename TString, typename TId>
inline typename Reference<StringSet<TString, Dependent<Tight> > >::Type
getValueById(StringSet<TString, Dependent<Tight> > & me,
            TId const id)
{
    SEQAN_ASSERT_GT_MSG(me.id_pos_map.count(id), 0u, "String id must be known!");
    return (value(me, me.id_pos_map.find(id)->second));
}

// --------------------------------------------------------------------------
// Function assignValueById()
// --------------------------------------------------------------------------

template<typename TString, typename TString2>
inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
assignValueById(StringSet<TString, Dependent<Tight> >& me,
                TString2& obj)
{
    SEQAN_CHECKPOINT;
    appendValue(me, obj);
    SEQAN_ASSERT_EQ(length(me.limits), length(me) + 1);
    return positionToId(me, length(me.strings) - 1);
}


template<typename TString, typename TId1>
inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
assignValueById(StringSet<TString, Dependent<Tight> >& me,
                TString& obj,
                TId1 id)
{
    SEQAN_CHECKPOINT;
    typedef StringSet<TString, Dependent<Tight> > TStringSet;
    typedef typename TStringSet::TIdPosMap::const_iterator TIter;
    typedef typename Id<TStringSet>::Type TId;

    if (me.lastId < (TId) id) me.lastId = (TId) (id + 1);

    TIter pos = me.id_pos_map.find(id);
    if (pos != me.id_pos_map.end()) {
        me.strings[pos->second] = &obj;
        me.limitsValid = false;
        return id;
    }
    appendValue(me.strings, &obj);
    appendValue(me.ids, id);
    me.id_pos_map.insert(std::make_pair(id, length(me.strings) - 1));
    appendValue(me.limits, lengthSum(me) + length(obj));
    return id;
}

// --------------------------------------------------------------------------
// Function removeValueById()
// --------------------------------------------------------------------------

template<typename TString, typename TId>
inline void
removeValueById(StringSet<TString, Dependent<Tight> >& me, TId const id)
{
    SEQAN_CHECKPOINT;
    typedef StringSet<TString, Dependent<Tight> > TStringSet;
    typedef typename Size<TStringSet>::Type TSize;
    typedef typename TStringSet::TIdPosMap::iterator TIter;

    SEQAN_ASSERT_EQ(length(me.limits), length(me) + 1);
    TIter pos = me.id_pos_map.find(id);
    if (pos != me.id_pos_map.end()) {
        TSize remPos = pos->second;
        erase(me.strings, remPos);
        erase(me.ids, remPos);
        me.id_pos_map.erase(pos);
        resize(me.limits, length(me.limits) - 1, Generous());

        for(TIter itChange = me.id_pos_map.begin(); itChange != me.id_pos_map.end(); ++itChange) {
            if (itChange->second > remPos) --(itChange->second);
        }
    }
    SEQAN_ASSERT_EQ(length(me.limits), length(me) + 1);
}

// --------------------------------------------------------------------------
// Function positionToId()
// --------------------------------------------------------------------------

template <typename TString, typename TPos>
inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
positionToId(StringSet<TString, Dependent<Tight> > & me,
            TPos const pos)
{
    SEQAN_CHECKPOINT;
    return me.ids[pos];
}

template <typename TString, typename TPos>
inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
positionToId(StringSet<TString, Dependent<Tight> > const & me,
            TPos const pos)
{
    SEQAN_CHECKPOINT;
    return me.ids[pos];
}

// --------------------------------------------------------------------------
// Function idToPosition()
// --------------------------------------------------------------------------

template <typename TString, typename TId>
inline typename Id<StringSet<TString, Dependent<Tight> > >::Type
idToPosition(StringSet<TString, Dependent<Tight> > const & me,
            TId const id)
{
    SEQAN_CHECKPOINT;
    return me.id_pos_map.find(id)->second;
/*
    for(unsigned i = 0; i < length(me.ids); ++i)
        if ((TId) me.ids[i] == id)
            return i;
    return 0;
    */
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_SET_DEPENDENT_TIGHT_H_
