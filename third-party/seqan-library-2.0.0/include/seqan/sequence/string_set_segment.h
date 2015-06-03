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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Implementation of Segment StringSet, a string set storing segments of
// another container.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_SET_SEGMENT_H_
#define SEQAN_SEQUENCE_STRING_SET_SEGMENT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct HostMember_;

template <typename THost>
struct StringSetPositions;

// ============================================================================
// Classes
// ============================================================================

// --------------------------------------------------------------------------
// Class StringSet<Segment>
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
class StringSet<THost, Segment<TSpec> >
{
public:
    typedef typename Member<StringSet, HostMember_>::Type   THostMember;
    typedef typename StringSetPositions<THost>::Type        TPositions;
    typedef typename StringSetLimits<StringSet>::Type       TLimits;
    typedef typename Concatenator<StringSet>::Type          TConcatenator;


    THostMember     data_host;
    TPositions      positions;
    TLimits         limits;
    bool            limitsValid;
    TConcatenator   concat;

    StringSet() :
        data_host(),
        limitsValid(true)
    {
        _initStringSetLimits(*this);
    }

    StringSet(THost & _host) :
        data_host(_toPointer(_host)),
        limitsValid(true)
    {
        _initStringSetLimits(*this);
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
};

// ============================================================================
// Members
// ============================================================================

// --------------------------------------------------------------------------
// Member HostMember_
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Member<StringSet<THost, Segment<TSpec> >, HostMember_>
{
    typedef typename IfView<THost, THost, THost *>::Type Type;
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction StringSetPositions
// --------------------------------------------------------------------------

template <typename THost>
struct StringSetPositions
{
    typedef typename Position<THost>::Type                  TPos_;
    typedef typename StringSpec<THost>::Type                TSpec_;

    typedef String<TPos_, TSpec_>                           Type;
};

template <typename THost>
struct StringSetPositions<THost const>
{
    typedef typename StringSetPositions<THost>::Type const  Type;
};

template <typename THost, typename TSpec>
struct StringSetPositions<StringSet<THost, TSpec> >
{
    typedef StringSet<THost, TSpec>                         TStringSet_;
    typedef typename StringSetPosition<TStringSet_>::Type   TPos_;
    typedef typename StringSpec<TStringSet_>::Type          TSpec_;

    typedef String<TPos_, TSpec_>                           Type;
};

// --------------------------------------------------------------------------
// Metafunction Value
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Host<StringSet<THost, Segment<TSpec> > >
{
    typedef THost Type;
};

template <typename THost, typename TSpec>
struct Host<StringSet<THost, Segment<TSpec> > const> :
    Host<StringSet<THost, Segment<TSpec> > > {};

// --------------------------------------------------------------------------
// Metafunction Value
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Value<StringSet<THost, Segment<TSpec> > >
    : Infix<THost> {};

template <typename THost, typename TSpec>
struct Value<StringSet<THost, Segment<TSpec> > const>
    : Infix<THost const> {};

// --------------------------------------------------------------------------
// Metafunction GetValue
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct GetValue<StringSet<THost, Segment<TSpec> > >
{
    typedef typename Infix<THost>::Type const Type;
};

template <typename THost, typename TSpec>
struct GetValue<StringSet<THost, Segment<TSpec> > const>
{
    typedef typename Infix<THost const>::Type const Type;
};

// --------------------------------------------------------------------------
// Metafunction Reference
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Reference<StringSet<THost, Segment<TSpec> > >
    : Infix<THost> {};

template <typename THost, typename TSpec>
struct Reference<StringSet<THost, Segment<TSpec> > const>
    : Infix<THost const> {};

// --------------------------------------------------------------------------
// Metafunction Prefix
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Prefix<StringSet<THost, Segment<TSpec> > >
    : Infix<THost> {};

template <typename THost, typename TSpec>
struct Prefix<StringSet<THost, Segment<TSpec> > const>
    : Infix<THost const> {};

// --------------------------------------------------------------------------
// Metafunction Suffix
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Suffix<StringSet<THost, Segment<TSpec> > >
    : Infix<THost> {};

template <typename THost, typename TSpec>
struct Suffix<StringSet<THost, Segment<TSpec> > const>
    : Infix<THost const> {};

// --------------------------------------------------------------------------
// Metafunction Infix
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Infix<StringSet<THost, Segment<TSpec> > >
    : Infix<THost> {};

template <typename THost, typename TSpec>
struct Infix<StringSet<THost, Segment<TSpec> > const >
    : Infix<THost const> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _refreshStringSetLimits()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TThreading>
SEQAN_HOST_DEVICE void
_refreshStringSetLimits(StringSet<THost, Segment<TSpec> > & me, Tag<TThreading> const & tag)
{
    partialSum(me.limits, tag);
    me.limitsValid = true;
}

template <typename THost, typename TSpec>
SEQAN_HOST_DEVICE void
_refreshStringSetLimits(StringSet<THost, Segment<TSpec> > & me)
{
    _refreshStringSetLimits(me, Serial());
}

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec>
SEQAN_HOST_DEVICE inline typename Parameter_<THost>::Type
host(StringSet<THost, Segment<TSpec> > & me)
{
    return _toParameter<THost>(me.data_host);
}

template <typename THost, typename TSpec>
SEQAN_HOST_DEVICE inline typename Parameter_<THost>::Type
host(StringSet<THost, Segment<TSpec> > const & me)
{
    return _toParameter<THost>(me.data_host);
}

// ----------------------------------------------------------------------------
// Function setHost()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline void
setHost(StringSet<THost, Segment<TSpec> > & me, typename Parameter_<THost>::Type _host)
{
    me.data_host = _toPointer(_host);
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec>
typename View<StringSet<THost, Segment<TSpec> > >::Type
view(StringSet<THost, Segment<TSpec> > & me)
{
    typename View<StringSet<THost, Segment<TSpec> > >::Type stringSetView;

    stringSetView.data_host = view(host(me));
    stringSetView.positions = view(me.positions);
    stringSetLimits(stringSetView) = view(stringSetLimits(me));
    stringSetView.limitsValid = me.limitsValid;

    // TODO(esiragusa): define view of concat.
//    concat(stringSetView) = view(concat(me));

    return stringSetView;
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TString2, typename TSpec2>
void assign(StringSet<THost, Segment<TSpec> > & me,
            StringSet<TString2, Segment<TSpec2> > const & other)
{
    setHost(me, host(other));
    assign(me.positions, other.positions);
    assign(stringSetLimits(me), stringSetLimits(other));
    assign(me.limitsValid, other.limitsValid);
    assign(concat(me), concat(other));
}

template <typename THost, typename TSpec, typename TString2, typename TSpec2>
void assign(StringSet<THost, Segment<TSpec> > & me,
            StringSet<TString2, Segment<TSpec2> > & other)
{
    assign(me, reinterpret_cast<StringSet<TString2, Segment<TSpec2> > const &> (other));
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TPos, typename TSegment>
inline void assignValue(StringSet<THost, Segment<TSpec> > & /* me */,
                        TPos /* pos */,
                        TSegment const & /* seq */)
{
    // NOTE(esiragusa): The host global position would be lost in the Segment.
    SEQAN_FAIL("Not implemented");
}

// --------------------------------------------------------------------------
// Function appendValue()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TSegment, typename TExpand>
inline void appendValue(StringSet<THost, Segment<TSpec> > & /* me */,
                        TSegment const & /* obj */,
                        Tag<TExpand> /* tag */)
{
    // NOTE(esiragusa): The host global position would be lost in the Segment.
    SEQAN_FAIL("Not implemented");
}

// --------------------------------------------------------------------------
// Function assignInfixWithLength()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TPos, typename TInfixPos, typename TSize>
SEQAN_HOST_DEVICE inline void
assignInfixWithLength(StringSet<THost, Segment<TSpec> > & me,
                      TPos pos, TInfixPos infixPos, TSize len)
{
    SEQAN_ASSERT_LT(pos, length(me));
    assignValue(me.positions, pos, infixPos);
    me.limitsValid = false;
    // NOTE(esiragusa): this breaks the invariant on limits.
    assignValue(me.limits, pos + 1, len);
//    assignValue(me.limits, pos + 1, me.limits[pos] + len);
}

// --------------------------------------------------------------------------
// Function appendInfixWithLength()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TPos, typename TSize, typename TExpand>
SEQAN_HOST_DEVICE inline void
appendInfixWithLength(StringSet<THost, Segment<TSpec> > & me,
                      TPos pos, TSize length, Tag<TExpand> tag)
{
    appendValue(me.positions, pos, tag);
    appendValue(me.limits, lengthSum(me) + length, tag);
}

// --------------------------------------------------------------------------
// Function appendInfix()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TPos, typename TExpand>
SEQAN_HOST_DEVICE inline void
appendInfix(StringSet<THost, Segment<TSpec> > & me,
            TPos posBegin, TPos posEnd, Tag<TExpand> tag)
{
    SEQAN_ASSERT_EQ(getSeqNo(posBegin), getSeqNo(posEnd));
    SEQAN_ASSERT_LEQ(getSeqOffset(posBegin), getSeqOffset(posEnd));

    appendValue(me.positions, posBegin, tag);
    appendValue(me.limits, lengthSum(me) + getSeqOffset(posEnd) - getSeqOffset(posBegin), tag);
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline void clear(StringSet<THost, Segment<TSpec> > & me)
{
//    me.data_host = NULL;
    resize(me.positions, 0, Exact());
    resize(me.limits, 1, Exact());
    me.limitsValid = true;
}

// --------------------------------------------------------------------------
// Function length()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
SEQAN_HOST_DEVICE inline typename Size<StringSet<THost, Segment<TSpec> > >::Type
length(StringSet<THost, Segment<TSpec> > const & me)
{
    return length(me.limits) - 1;
}

// --------------------------------------------------------------------------
// Function resize()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TSize, typename TExpand>
inline typename Size<StringSet<THost, Segment<TSpec> > >::Type
resize(StringSet<THost, Segment<TSpec> > & me, TSize new_size, Tag<TExpand> tag)
{
    resize(me.positions, new_size, tag);

    if (new_size < length(me.limits))
        return resize(me.limits, new_size + 1, tag) - 1;
    else
        return resize(me.limits, new_size + 1, back(me.limits), tag) - 1;
}

// --------------------------------------------------------------------------
// Function reserve()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TSize, typename TExpand>
inline typename Size<StringSet<THost, Segment<TSpec> > >::Type
reserve(StringSet<THost, Segment<TSpec> > & me,
        TSize const & new_capacity,
        Tag<TExpand> tag)
{
    reserve(me.positions, new_capacity, tag);
    return reserve(me.limits, new_capacity + 1, tag) - 1;
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TPos >
SEQAN_HOST_DEVICE inline typename Infix<THost>::Type
value(StringSet<THost, Segment<TSpec> > & me, TPos pos)
{
    SEQAN_ASSERT_NOT(empty(me));
    return infixWithLength(host(me), me.positions[pos], me.limits[pos + 1] - me.limits[pos]);
}

template <typename THost, typename TSpec, typename TPos >
SEQAN_HOST_DEVICE inline typename Infix<THost const>::Type
value(StringSet<THost, Segment<TSpec> > const & me, TPos pos)
{
    SEQAN_ASSERT_NOT(empty(me));
    return infixWithLength(host(me), me.positions[pos], me.limits[pos + 1] - me.limits[pos]);
}

// --------------------------------------------------------------------------
// Function swap()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
void swap(StringSet<THost, Segment<TSpec> > & lhs,
          StringSet<THost, Segment<TSpec> > & rhs)
{
    using std::swap;

    swap(lhs.data_host, rhs.data_host);
    swap(lhs.positions, rhs.positions);
    swap(lhs.limits, rhs.limits);
    swap(lhs.limitsValid, rhs.limitsValid);
    swap(lhs.concat, rhs.concat);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_SET_SEGMENT_H_
