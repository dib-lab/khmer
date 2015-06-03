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
// Implementation of the Packed String class.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_PACKED_H_
#define SEQAN_SEQUENCE_STRING_PACKED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T>
struct HostIterator;

//struct Device_ {};
//typedef Tag<Device_> Device;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Specialization Packed String
// --------------------------------------------------------------------------

/*!
 * @class PackedString Packed String
 * @extends String
 * @headerfile <seqan/sequence.h>
 * @brief A string that stores as many values in one machine word as possible.
 *
 * @signature template <typename TValue, typename THostSpec>
 *            class String<TValue, Packed<THostSpec> >;
 *
 * @tparam TValue The value type, that is the type of the items/characters
 *                stored in the string.Use @link Value @endlink to get the value
 *                type for a given class.
 * @tparam THostSpec The specializing type.This is the specialization of the
 *                   host string that is used for storing the packed values.
 *                   Default: @link AllocString @endlink
 */

template <typename THostspec = Alloc<> >
struct Packed;

// --------------------------------------------------------------------------
// Metafunction PackedTraits_
// --------------------------------------------------------------------------

template <typename TValue, int BIT_WIDTH, int NUM>
struct FillMultiplierRecursion_
{
    static const TValue VALUE = (FillMultiplierRecursion_<TValue, BIT_WIDTH, NUM - 1>::VALUE << BIT_WIDTH) | (TValue)1;
};

template <typename TValue, int BIT_WIDTH>
struct FillMultiplierRecursion_<TValue, BIT_WIDTH, 0>
{
    static const TValue VALUE = 0;
};

template <typename TPackedString>
struct PackedTraits_
{
    typedef typename Size<TPackedString>::Type                          TSize;
    typedef typename Value<TPackedString>::Type                         TValue;
    typedef typename Value<typename Host<TPackedString>::Type>::Type    THostValue;

    enum
    {
        BITS_PER_VALUE = BitsPerValue<TValue>::VALUE,
        VALUES_PER_HOST_VALUE = LENGTH<THostValue>::VALUE,
        WASTED_BITS = 64 - BITS_PER_VALUE * VALUES_PER_HOST_VALUE
    };

    static inline
    typename Size<typename Host<TPackedString>::Type>::Type
    toHostLength(typename Size<TPackedString>::Type len)
    {
        return (len + VALUES_PER_HOST_VALUE - 1) / VALUES_PER_HOST_VALUE;
    }

    static inline
    typename Size<typename Host<TPackedString>::Type>::Type
    toHostLength1(typename Size<TPackedString>::Type len)
    {
        return (len == 0)? 0: 1 + (len + VALUES_PER_HOST_VALUE - 1) / VALUES_PER_HOST_VALUE;
    }
};

/*???TODO Optimierungsm�glichkeiten:
- _clearSpace kopiert Zeichenweise im Packed-String, und nicht im Host-String
- _clearSpace verwendet resize, um den Host zu vergr��ern, d.h. der Inhalt wird eventuell doppelt kopiert.
*/

template <typename TValue, typename THostspec>
class String<TValue, Packed<THostspec> >
{
public:
    typedef typename Host<String>::Type THost;
    typedef typename Size<String>::Type TSize;
    typedef PackedTraits_<String>       TTraits;

    THost data_host;

    String()
    {
    }

    template <typename TSource>
    String(TSource & source)
    {
        reserve(*this, capacity(source), Exact());
        assign(*this, source);
    }
    template <typename TSource>
    String(TSource const & source)
    {
        reserve(*this, capacity(source), Exact());
        assign(*this, source);
    }

    template <typename TSource>
    String & operator =(TSource const & source)
    {
        assign(*this, source);
        return *this;
    }
    String & operator =(String const & source)
    {
        assign(*this, source);
        return *this;
    }


    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    template <typename TPos>
    inline typename Reference<String>::Type
    operator[](TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<String const>::Type
    operator[](TPos pos) const
    {
        return data_host[1 + pos / TTraits::VALUES_PER_HOST_VALUE][pos % TTraits::VALUES_PER_HOST_VALUE];
    }
};

// ----------------------------------------------------------------------------
// Functor FunctorTestAllZeros
// ----------------------------------------------------------------------------

template <typename T>
struct FunctorTestAllZeros
{};

template <typename THostSpec>
struct FunctorTestAllZeros<String<bool, Packed<THostSpec> > >
{
    typedef String<bool, Packed<THostSpec> > TPackedString;
    typedef typename Host<TPackedString>::Type TPackedHost;
    typedef typename Value<TPackedHost>::Type TPackedHostValue;
    typedef typename Size<TPackedHostValue>::Type TSize;

    TSize _wastedBits;

    template <typename TShift>
    FunctorTestAllZeros(TShift const & shift) : _wastedBits(shift)
    {}

    template <typename TValue>
    inline bool operator()(TValue const & val, False const & /*tag*/) const
    {
        return testAllZeros(val);
    }

    template <typename TValue>
    inline bool operator()(TValue const & val, True const & /*tag*/) const
    {
        return testAllZeros(val >> _wastedBits);
    }
};

// ----------------------------------------------------------------------------
// Functor FunctorTestAllOnes
// ----------------------------------------------------------------------------

template <typename T>
struct FunctorTestAllOnes
{};

template <typename THostSpec>
struct FunctorTestAllOnes<String<bool, Packed<THostSpec> > >
{
    typedef String<bool, Packed<THostSpec> > TPackedString;
    typedef typename Host<TPackedString>::Type TPackedHost;
    typedef typename Value<TPackedHost>::Type TPackedHostValue;
    typedef typename TPackedHostValue::TBitVector TBitVector;
    typedef typename Size<TPackedHostValue>::Type TSize;

    TSize _wastedBits;

    template <typename TShift>
    FunctorTestAllOnes(TShift const & shift) : _wastedBits(shift)
    {}

    template <typename TValue>
    inline bool operator()(TValue const & val, False const & /*tag*/) const
    {
        return testAllOnes(val);
    }

    template <typename TValue>
    inline bool operator()(TValue const & val, True const & /*tag*/) const
    {
        TBitVector maskWastedBits = (1 << _wastedBits) - 1;  // We only need the mask if we are in the last value.
        return testAllOnes(val | maskWastedBits);
    }
};

// --------------------------------------------------------------------------
// Specialization Packed String Iter
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
class Iter<TPackedString, Packed<THostspec> >
{
public:
    typedef typename Host<Iter>::Type THostIterator;
    typedef typename Position<TPackedString>::Type TPosition;
    typedef PackedTraits_<TPackedString> TTraits;

    THostIterator data_iterator;
    unsigned char localPos;

    Iter()
    {
    }

    Iter(THostIterator data_iterator_):
          data_iterator(data_iterator_),
          localPos(0)
    {
    }

    Iter(THostIterator data_iterator_, TPosition localPos_):
          data_iterator(data_iterator_),
          localPos(localPos_)
    {
    }

    Iter(TPackedString &container):
          data_iterator(begin(host(container), Standard()) + 1),
          localPos(0)
    {
    }

    Iter(TPackedString &container, TPosition pos):
          data_iterator(begin(host(container), Standard()) + 1 + pos / TTraits::VALUES_PER_HOST_VALUE),
          localPos(pos % TTraits::VALUES_PER_HOST_VALUE)
    {
    }

//    inline
//    Iter const &
//    operator=(Iter const & other_)
//    {
//        data_iterator = other_.data_iterator;
//        localPos = other_.localPos;
//        return *this;
//    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction StringSpec
// --------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct StringSpec<String<TValue, Packed<TSpec> > > : StringSpec<String<TValue, TSpec> > {};

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct DefaultOverflowImplicit<String<TValue, Packed<THostspec> > >
        : DefaultOverflowImplicit<typename Host<String<TValue, Packed<THostspec> > >::Type>
{};

template <typename TValue, typename THostspec>
struct DefaultOverflowImplicit<String<TValue, Packed<THostspec> > const>
        : DefaultOverflowImplicit<typename Host<String<TValue, Packed<THostspec> > const>::Type>
{};

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowExplicit
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct DefaultOverflowExplicit<String<TValue, Packed<THostspec> > >
        : DefaultOverflowExplicit<typename Host<String<TValue, Packed<THostspec> > >::Type>
{};

template <typename TValue, typename THostspec>
struct DefaultOverflowExplicit<String<TValue, Packed<THostspec> > const>
        : DefaultOverflowExplicit<typename Host<String<TValue, Packed<THostspec> > const>::Type>
{};

// --------------------------------------------------------------------------
// Metafunction IsContiguous
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct IsContiguous<String<TValue, Packed<THostspec> > >
{
    typedef False Type;
    enum { VALUE = false };
};

// --------------------------------------------------------------------------
// Metafunction PackedHostValue_
// --------------------------------------------------------------------------

template <typename TString>
struct PackedHostValue_
{
    typedef typename Value<TString>::Type TValue;
    typedef Tuple<TValue, 64 / BitsPerValue<TValue>::VALUE, BitPacked<> > Type; // use 64bit words
};

// --------------------------------------------------------------------------
// Metafunction Host
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct Host<String<TValue, Packed<THostspec> > >
{
    typedef typename PackedHostValue_<String<TValue, Packed<THostspec> > >::Type TInternalValue;
    typedef String<TInternalValue, THostspec> Type;
};


#ifdef PLATFORM_CUDA
template <typename TValue, typename TSpec>
struct Host<String<TValue, Packed<Device<TSpec> > > >
{
    typedef typename PackedHostValue_<String<TValue, Packed<Device<TSpec> > > >::Type TInternalValue;
    typedef thrust::device_vector<TInternalValue> Type;
};

template <typename TValue, typename TSpec>
struct PackedHostValue_<String<TValue, Packed<Device<TSpec> > > >
{
    typedef Tuple<TValue, 32 / BitsPerValue<TValue>::VALUE, BitPacked<> > Type;
};
#endif

template <typename TValue, typename THostspec>
struct Host<String<TValue, Packed<THostspec> > const>
{
    typedef typename Host<String<TValue, Packed<THostspec> > >::Type const Type;
};

// --------------------------------------------------------------------------
// Metafunction GetValue
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct GetValue<String<TValue, Packed<THostspec> > > :
    public Value<String<TValue, Packed<THostspec> > > {};

template <typename TValue, typename THostspec>
struct GetValue<String<TValue, Packed<THostspec> > const> :
    public Value<String<TValue, Packed<THostspec> > const> {};

// --------------------------------------------------------------------------
// Metafunction Reference
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct Reference<String<TValue, Packed<THostspec> > >
{
    typedef typename Iterator<String<TValue, Packed<THostspec> >, Standard>::Type TIterator;
    typedef Proxy<IteratorProxy<TIterator> > Type;
};

template <typename TValue, typename THostspec>
struct Reference<String<TValue, Packed<THostspec> > const> :
    public GetValue<String<TValue, Packed<THostspec> > const> {};

// --------------------------------------------------------------------------
// Metafunction Size
// --------------------------------------------------------------------------

/*
template <typename TValue, typename THostspec>
struct Size<String<TValue, Packed<THostspec> > >
{
    typedef __int64 Type;
};
template <typename TValue, typename THostspec>
struct Size<String<TValue, Packed<THostspec> > const>
{
    typedef __int64 Type;
};
*/

// --------------------------------------------------------------------------
// Metafunction Iterator
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct Iterator<String<TValue, Packed<THostspec> >, Standard>
{
    typedef Iter<String<TValue, Packed<THostspec> >, Packed<THostspec> > Type;
};

template <typename TValue, typename THostspec>
struct Iterator<String<TValue, Packed<THostspec> > const, Standard>
{
    typedef Iter<String<TValue, Packed<THostspec> > const, Packed<THostspec> > Type;
};

// --------------------------------------------------------------------------
// Metafunction Host
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
struct Host<Iter<TPackedString, Packed<THostspec> > > :
    public Iterator<typename Host<TPackedString>::Type, Standard> {};

template <typename TPackedString, typename THostspec>
struct Host<Iter<TPackedString, Packed<THostspec> > const>
{
    typedef typename Host<TPackedString>::Type THost_;
    typedef typename Iterator<THost_, Standard>::Type const Type;
};

// --------------------------------------------------------------------------
// Internal Metafunction TempCopy_
// --------------------------------------------------------------------------

// Note: this works only, if the copy assignment is done without using TempCopy_.
template <typename TValue, typename THostspec>
struct TempCopy_<String<TValue, Packed<THostspec> > >
{
    typedef String<TValue, Packed<THostspec> > Type;
};

// ============================================================================
// Functions
// ============================================================================

// ****************************************************************************
// Functions for Packed String
// ****************************************************************************

// ----------------------------------------------------------------------------
// Function std::swap()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline void
swap(String<TValue, Packed<THostspec> > & a,
     String<TValue, Packed<THostspec> > & b)
{
    std::swap(a.data_host, b.data_host);
}

// --------------------------------------------------------------------------
// Function host
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline typename Host<String<TValue, Packed<THostspec> > >::Type &
host(String<TValue, Packed<THostspec> > & me)
{
    return me.data_host;
}

template <typename TValue, typename THostspec>
inline typename Host<String<TValue, Packed<THostspec> > const>::Type const &
host(String<TValue, Packed<THostspec> > const & me)
{
    return me.data_host;
}

// --------------------------------------------------------------------------
// Function length
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline typename Size<String<TValue, Packed<THostspec> > const>::Type
length(String<TValue, Packed<THostspec> > const & me)
{
    if (empty(host(me)))
        return 0;
    else
        return front(host(me)).i;
}

// --------------------------------------------------------------------------
// Function _setLength
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TSize>
inline void
_setLength(
    String<TValue, Packed<THostspec> > & me,
    TSize new_length)
{
    typedef String<TValue, Packed<THostspec> > TString;
    if (new_length == 0)
    {
        _setLength(host(me), 0);
        return;
    }
    _setLength(host(me), 1 + PackedTraits_<TString>::toHostLength(new_length));
    front(host(me)).i = new_length;
}

template <typename TValue, typename THostspec, typename TSize, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type
resize(
    String<TValue, Packed<THostspec> > & me,
    TSize new_length,
    Tag<TExpand> tag)
{
    typedef String<TValue, Packed<THostspec> > TString;
    typedef PackedTraits_<TString> TTraits;
    typedef typename Size<TString>::Type TStringSize;

    if (new_length == 0)
    {
        clear(host(me));
        return 0;
    }
    TStringSize max_length = (resize(host(me), TTraits::toHostLength(new_length) + 1, tag) - 1) * TTraits::VALUES_PER_HOST_VALUE;
    if ((TStringSize)new_length > max_length)
        new_length = max_length;
    return front(host(me)).i = new_length;
}

// --------------------------------------------------------------------------
// Function assign()
//
// Helpers: _assignCopyPackedString()
// --------------------------------------------------------------------------

// optimized variant for copy assignment. The host sequence is copied instead of
// copying the packed string value by value.
template <typename TTarget, typename TSource, typename TTag>
inline void
_assignCopyPackedString(TTarget & target,
                        TSource & source,
                        Tag<TTag> const & tag)
{
    typedef typename Size<TTarget>::Type TSize;

    assign(host(target), host(source), tag);
    if (empty(host(target)))
        return;
    TSize new_length_limit = (length(host(target)) - 1) * PackedTraits_<TTarget>::VALUES_PER_HOST_VALUE;
    _setLength(target, _min((TSize)length(source), new_length_limit));
}

template <typename TTarget, typename TSource, typename TSize, typename TTag>
inline void
_assignCopyPackedString(TTarget & target,
                        TSource & source,
                        TSize limit,
                        Tag<TTag> const & tag)
{
    typedef typename Size<TTarget>::Type TSize2;

    TSize2 host_limit = PackedTraits_<TTarget>::toHostLength1(limit);
    assign(host(target), host(source), host_limit, tag);
    if (empty(host(target)))
        return;
    TSize2 new_length_limit = (length(host(target)) - 1) * PackedTraits_<TTarget>::VALUES_PER_HOST_VALUE;
    _setLength(target, _min((TSize2)length(source), _min(new_length_limit, (TSize2)limit)));
}

template <typename TValue, typename THostspec, typename TTag>
inline void
assign(String<TValue, Packed<THostspec> > & target,
       String<TValue, Packed<THostspec> > & source,
       Tag<TTag> const & tag)
{
    _assignCopyPackedString(target, source, tag);
}

template <typename TValue, typename THostspec, typename TTag>
inline void
assign(String<TValue, Packed<THostspec> > & target,
       String<TValue, Packed<THostspec> > const & source,
       Tag<TTag> const & tag)
{
    _assignCopyPackedString(target, source, tag);
}

template <typename TValue, typename THostspec, typename TSize, typename TTag>
void assign(String<TValue, Packed<THostspec> > & target,
            String<TValue, Packed<THostspec> > & source,
            TSize limit,
            Tag<TTag> const & tag)
{
    _assignCopyPackedString(target, source, limit, tag);
}
template <typename TValue, typename THostspec, typename TSize, typename TTag>
void assign(String<TValue, Packed<THostspec> > & target,
            String<TValue, Packed<THostspec> > const & source,
            TSize limit,
            Tag<TTag> const & tag)
{
    _assignCopyPackedString(target, source, limit, tag);
}

// --------------------------------------------------------------------------
// Function getObjectId()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline void const *
getObjectId(String<TValue, Packed<THostspec> > const & me)
{
    return getObjectId(host(me));
}

// --------------------------------------------------------------------------
// Function iter()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TPos>
inline typename Iterator<String<TValue, Packed<THostspec> >, Standard>::Type
iter(String<TValue, Packed<THostspec> > & me,
     TPos pos,
     Standard)
{
    typedef typename Iterator<String<TValue, Packed<THostspec> >, Standard>::Type TIterator;
    return TIterator(me, pos);
}

template <typename TValue, typename THostspec, typename TPos>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Standard>::Type
iter(String<TValue, Packed<THostspec> > const & me,
     TPos pos,
     Standard)
{
    typedef typename Iterator<String<TValue, Packed<THostspec> > const, Standard>::Type TIterator;
    return TIterator(me, pos);
}

// --------------------------------------------------------------------------
// Function begin()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline typename Iterator<String<TValue, Packed<THostspec> >, Standard>::Type
begin(String<TValue, Packed<THostspec> > & me,
      Standard)
{
    typedef typename Iterator<String<TValue, Packed<THostspec> >, Standard>::Type TIterator;
    return TIterator(me);
}

template <typename TValue, typename THostspec>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Standard>::Type
begin(String<TValue, Packed<THostspec> > const & me,
      Standard)
{
    typedef typename Iterator<String<TValue, Packed<THostspec> > const, Standard>::Type TIterator;
    return TIterator(me);
}

// --------------------------------------------------------------------------
// Function end()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline typename Iterator<String<TValue, Packed<THostspec> >, Standard>::Type
end(String<TValue, Packed<THostspec> > & me,
    Standard)
{
    return iter(me, length(me), Standard());
}

template <typename TValue, typename THostspec>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Standard>::Type
end(String<TValue, Packed<THostspec> > const & me,
    Standard)
{
    return iter(me, length(me), Standard());
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TPos>
inline typename Reference<String<TValue, Packed<THostspec> > >::Type
value(String<TValue, Packed<THostspec> > & me,
      TPos const & pos)
{
    return *iter(me, pos, Standard());
}

template <typename TValue, typename THostspec, typename TPos>
inline typename Reference<String<TValue, Packed<THostspec> > const>::Type
value(String<TValue, Packed<THostspec> > const & me,
      TPos const & pos)
{
    typedef String<TValue, Packed<THostspec> > TPackedString;
    typedef PackedTraits_<TPackedString> TTraits;
    return me.data_host[1 + pos / TTraits::VALUES_PER_HOST_VALUE][pos % TTraits::VALUES_PER_HOST_VALUE];
}

// --------------------------------------------------------------------------
// Function capacity()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline typename Size<String<TValue, Packed<THostspec> > const>::Type
capacity(String<TValue, Packed<THostspec> > const & me)
{
    typedef String<TValue, Packed<THostspec> > TPackedString;
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename Size<TPackedString>::Type TSize;

    TSize cap = capacity(host(me));
    return (cap == 0)? 0: (cap - 1) * (TSize)TTraits::VALUES_PER_HOST_VALUE;
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline void
clear(String<TValue, Packed<THostspec> > & me)
{
    clear(host(me));
}

// --------------------------------------------------------------------------
// Function shrinkToFit()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline void
shrinkToFit(String<TValue, Packed<THostspec> > & me)
{
    shrinkToFit(host(me));
}

// --------------------------------------------------------------------------
// Function _clearUnusedBits()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline void
_clearUnusedBits(String<TValue, Packed<THostspec> > & me)
{
    typedef String<TValue, Packed<THostspec> > TPackedString;
    typedef PackedTraits_<TPackedString> TTraits;

    typedef typename Host<TPackedString>::Type THost;
    typedef typename Iterator<THost, Standard>::Type THostIterator;

    typedef typename TTraits::THostValue THostValue;
    typedef typename THostValue::TBitVector TBitVector;
    typedef typename Size<TPackedString>::Type TSize;

    static const TBitVector ALL_ONE = ~(TBitVector)0 >> TTraits::WASTED_BITS;

    THostIterator it = begin(host(me), Standard());
    THostIterator itEnd = end(host(me), Standard());

    if (it == itEnd)
        return;

    // zero all wasted bits
    if (TTraits::WASTED_BITS != 0)
        for (++it; it != itEnd; ++it)
            it->i &= ALL_ONE;

    // zero the last half-filled tuple
    TSize lastValues = length(me) % TTraits::VALUES_PER_HOST_VALUE;
    if (lastValues != 0)
        back(host(me)).i &= ALL_ONE & ~(ALL_ONE >> (lastValues * TTraits::BITS_PER_VALUE));
}

/*

template<typename TTarget, typename TSource1, typename TSource2>
inline void
_arrayConstructCopyDefault(TSource1 source_begin,
                           TSource2 source_end,
                           TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    while (source_begin != source_end)
    {
        // NOTE(holtgrew): getValue() is used here since value() could return
        // a proxy!
        valueConstruct(target_begin, getValue(source_begin));
        ++source_begin;
        ++target_begin;
    }
}

*/

// --------------------------------------------------------------------------
// Function arrayCopyForward()
// --------------------------------------------------------------------------

// Copy segment of a source packed strings into a target packed string
//
// ------------------------------------
// |      |      | aaaBB|CCCCdd|d     | source
// ------------------------------------
//
// ------------------------------------
// |      |   aaa|BBCCCC|ddd   |      | target
// ------------------------------------

template < typename TPackedString, typename TSpec >
inline void
arrayCopyForward(Iter<TPackedString, Packed<TSpec> > source_begin,
                 Iter<TPackedString, Packed<TSpec> > source_end,
                 Iter<TPackedString, Packed<TSpec> > target_begin)
{
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename TTraits::THostValue THostValue;
    typedef typename Size<TPackedString>::Type TSize;
    typedef typename Host<Iter<TPackedString, Packed<TSpec> > >::Type THostIter;

    TSize size = source_end - source_begin;

    // will we touch more than one word in the target string?
    if (size > (TSize)(TTraits::VALUES_PER_HOST_VALUE - target_begin.localPos))
    {
        // first, copy all values to a target word boundary (a's in figure above)
        if (target_begin.localPos != 0)
        {
            while (target_begin.localPos != TTraits::VALUES_PER_HOST_VALUE)
            {
                assignValue(target_begin, getValue(source_begin));
                ++source_begin;
                ++target_begin.localPos;
            }
            target_begin.localPos = 0;
            ++hostIterator(target_begin);
        }

        // now copy whole words
        int leftShift = source_begin.localPos;
        if (leftShift != 0)
        {
            THostIter source_lastWord = hostIterator(source_end);
            if (source_end.localPos < source_begin.localPos)
                --source_lastWord;

            typename THostValue::TBitVector prevWord = hostIterator(source_begin)->i;
            int rightShift = TTraits::VALUES_PER_HOST_VALUE - leftShift;
            leftShift *= TTraits::BITS_PER_VALUE;
            rightShift *= TTraits::BITS_PER_VALUE;
            for (; hostIterator(source_begin) != source_lastWord; ++hostIterator(target_begin))
            {
                // words must be shifted and or'ed (BB|CCCC in the figure)
                ++hostIterator(source_begin);
                typename THostValue::TBitVector curWord = hostIterator(source_begin)->i;
                hostIterator(target_begin)->i = (prevWord << leftShift) | (curWord >> rightShift);
                prevWord = curWord;
            }
        }
        else
        {
            // words need not to be shifted
            arrayCopyForward(hostIterator(source_begin), hostIterator(source_end), hostIterator(target_begin));
            hostIterator(target_begin) += hostIterator(source_end) - hostIterator(source_begin);
        }
    }

    // copy (less than VALUES_PER_HOST_VALUE many) remaining values (d's in figure above)
    while (source_begin.localPos != source_end.localPos)
    {
        assignValue(target_begin, getValue(source_begin));
        ++source_begin;
        ++target_begin.localPos;
    }
}

// --------------------------------------------------------------------------
// Function arrayCopyBackward()
// --------------------------------------------------------------------------

// Copy segment of a source packed strings into a target packed string
//
// ------------------------------------
// |      |      | aaaBB|CCCCdd|d     | source
// ------------------------------------
//
// ------------------------------------
// |      |      |   aaa|BBCCCC|ddd   | target
// ------------------------------------

template < typename TPackedString, typename TSpec >
inline void
arrayCopyBackward(Iter<TPackedString, Packed<TSpec> > source_begin,
                  Iter<TPackedString, Packed<TSpec> > source_end,
                  Iter<TPackedString, Packed<TSpec> > target_begin)
{
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename TTraits::THostValue THostValue;
    typedef typename Size<TPackedString>::Type TSize;
    typedef typename Host<Iter<TPackedString, Packed<TSpec> > >::Type THostIter;

    // iterator to the first whole word in the target
    THostIter target_firstWord = hostIterator(target_begin);
    if (target_begin.localPos != 0)
        ++target_firstWord;

    TSize size = source_end - source_begin;
    target_begin += size;

    // will we touch more than one word in the target string?
    if (size > (TSize)target_begin.localPos)
    {
        // first, copy all values to a target word boundary (a's in figure above)
        while (target_begin.localPos != 0)
        {
            --source_end;
            --target_begin.localPos;
            assignValue(target_begin, getValue(source_end));
        }

        // now copy whole words
        int leftShift = source_end.localPos;
        if (leftShift != 0)
        {
            typename THostValue::TBitVector prevWord = hostIterator(source_end)->i;
            int rightShift = TTraits::VALUES_PER_HOST_VALUE - leftShift;
            leftShift *= TTraits::BITS_PER_VALUE;
            rightShift *= TTraits::BITS_PER_VALUE;
            while (hostIterator(target_begin) != target_firstWord)
            {
                --hostIterator(source_end);
                --hostIterator(target_begin);

                // words must be shifted and or'ed (BB|CCCC in the figure)
                typename THostValue::TBitVector curWord = hostIterator(source_end)->i;
                hostIterator(target_begin)->i = (curWord << leftShift) | (prevWord >> rightShift);
                prevWord = curWord;
            }
        }
        else
        {
            // words need not to be shifted
            while (hostIterator(target_begin) != target_firstWord)
            {
                --hostIterator(source_end);
                --hostIterator(target_begin);
                assignValue(hostIterator(target_begin), getValue(hostIterator(source_end)));
            }
        }

        // now we are at the end of a word and the next --target_begin would
        // decrease hostIterator(target_begin) and set localPos to TTraits::VALUES_PER_HOST_VALUE - 1
        --hostIterator(target_begin);
        target_begin.localPos = TTraits::VALUES_PER_HOST_VALUE;
    }

    // copy (less than VALUES_PER_HOST_VALUE many) remaining values (d's in figure above)
    while (source_end.localPos != source_begin.localPos)
    {
        --source_end;
        --target_begin.localPos;
        assignValue(target_begin, getValue(source_end));
    }
}

template<typename TPackedString, typename TSpec, typename TValue2>
inline void
arrayFill(Iter<TPackedString, Packed<TSpec> > begin_,
          Iter<TPackedString, Packed<TSpec> > end_,
          TValue2 const & value)
{
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename TTraits::THostValue THostValue;
    typedef typename THostValue::TBitVector TBitVector;
    typedef typename Value<TPackedString>::Type TValue;

    THostValue fillValue;
    fillValue.i = FillMultiplierRecursion_<
        TBitVector,
        TTraits::BITS_PER_VALUE,
        TTraits::VALUES_PER_HOST_VALUE>::VALUE * ordValue((TValue)value);

    static const TBitVector ALL_ONE = ~(TBitVector)0 >> TTraits::WASTED_BITS;
    TBitVector maskB =   ALL_ONE >> (begin_.localPos * TTraits::BITS_PER_VALUE);        // begin mask (0..keep old value, 1..fill value)
    TBitVector maskE = ~(ALL_ONE >> (  end_.localPos * TTraits::BITS_PER_VALUE));       // end mask

    if (hostIterator(begin_) == hostIterator(end_))
    {
        maskE &= maskB;     // intersect masks if beginning and end are in the same word
    }
    else
    {
        if (begin_.localPos != 0)
        {
            hostIterator(begin_)->i = (hostIterator(begin_)->i & ~maskB) | (fillValue.i & maskB);   // blend with begin mask
            ++hostIterator(begin_);
        }
        arrayFill(hostIterator(begin_), hostIterator(end_), fillValue);                 // fill whole words
        if (end_.localPos == 0)                                                         // don't touch last word (if we don't need to)
            return;

    }
    hostIterator(end_)->i = (hostIterator(end_)->i & ~maskE) | (fillValue.i & maskE);   // blend with end mask
}


// TODO(weese): The following wrappers are boiler-plate code and should be not necessary in the new sequence interface

// --------------------------------------------------------------------------
// Function arrayMoveForward()
// --------------------------------------------------------------------------

// TODO(weese): There should be default wrappers using arrayCopyForward (so that these overloads are not necessary)
template < typename TPackedString, typename TSpec >
inline void
arrayMoveForward(Iter<TPackedString, Packed<TSpec> > source_begin,
                 Iter<TPackedString, Packed<TSpec> > source_end,
                 Iter<TPackedString, Packed<TSpec> > target_begin)
{
    arrayCopyForward(source_begin, source_end, target_begin);
}

// --------------------------------------------------------------------------
// Function arrayMoveBackward()
// --------------------------------------------------------------------------

// TODO(weese): There should be default wrappers using arrayCopyBackward (so that these overloads are not necessary)
template < typename TPackedString, typename TSpec >
inline void
arrayMoveBackward(Iter<TPackedString, Packed<TSpec> > source_begin,
                  Iter<TPackedString, Packed<TSpec> > source_end,
                  Iter<TPackedString, Packed<TSpec> > target_begin)
{
    arrayCopyBackward(source_begin, source_end, target_begin);
}

// --------------------------------------------------------------------------
// Function arrayConstructCopy()
// --------------------------------------------------------------------------

// TODO(weese): it should be not necessary to overload construct/destruct functions for POD/Simple types (IsSimple == true)
template < typename TPackedString, typename TSpec >
inline void
arrayConstructCopy(Iter<TPackedString, Packed<TSpec> > source_begin,
                   Iter<TPackedString, Packed<TSpec> > source_end,
                   Iter<TPackedString, Packed<TSpec> > target_begin)
{
    arrayCopyForward(source_begin, source_end, target_begin);
}

// TODO(weese): it should be not necessary to overload construct/destruct functions for POD/Simple types (IsSimple == true)
template < typename TPackedString, typename TSpec >
inline void
arrayConstruct(Iter<TPackedString, Packed<TSpec> > begin_,
               Iter<TPackedString, Packed<TSpec> > end_)
{
    // TODO(weese:) actually, I don't want this default zero-initialization for POD/Simple types
    // If someone still needs zero-fill resize, he/she should use the function below
    arrayFill(begin_, end_, typename Value<TPackedString>::Type());
}

template < typename TPackedString, typename TSpec, typename TParam >
inline void
arrayConstruct(Iter<TPackedString, Packed<TSpec> > begin_,
               Iter<TPackedString, Packed<TSpec> > end_,
               TParam const & param_)
{
    arrayFill(begin_, end_, param_);
}

// --------------------------------------------------------------------------
// Function arrayDestructCopy()
// --------------------------------------------------------------------------

// TODO(weese): it should be not necessary to overload construct/destruct functions for POD/Simple types (IsSimple == true)
template < typename TPackedString, typename TSpec >
inline void
arrayDestruct(Iter<TPackedString, Packed<TSpec> >,
              Iter<TPackedString, Packed<TSpec> >)
{
}

// --------------------------------------------------------------------------
// Function _clearSpace()
//
// Helper struct ClearSpaceStringPacked_.
// --------------------------------------------------------------------------

//implementation for all expand tags other than "limit"
template <typename TExpand>
struct ClearSpaceStringPacked_
{
    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size)
    {
        typedef typename Size<T>::Type TSize;
        TSize wanted_host_length = PackedTraits_<T>::toHostLength1(size);
        TSize new_host_length = resize(host(seq), wanted_host_length, TExpand());
        if (new_host_length == 0)
            return 0;
        if (new_host_length < wanted_host_length)
            size = (new_host_length - 1) * PackedTraits_<T>::VALUES_PER_HOST_VALUE;
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
        if (size > limit)
            size = limit;
        return _clearSpace_(seq, limit);
    }

    template <typename T>
    static inline typename Size<T>::Type
    _clearSpace_(
        T & seq,
        typename Size<T>::Type size,
        typename Size<T>::Type start,
        typename Size<T>::Type end)
    {
        return _clearSpace_(seq, size, start, end, maxValue<typename Size<T>::Type >());
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
//??? TODO: This function can be accelerated this way:
//              - move values in host
//              - avoid double moving of the rest-part if "resize" allocates a new block

        typedef typename Size<T>::Type TSize;

        TSize old_length = length(seq);
        TSize old_size = end - start;
        TSize wanted_new_length = _min(old_length + size - old_size, limit);
        TSize wanted_host_length = PackedTraits_<T>::toHostLength1(wanted_new_length);
        TSize new_host_length = resize(host(seq), wanted_host_length, TExpand());

        TSize new_length;
        if (new_host_length < wanted_host_length)
        {
            new_length = (new_host_length == 0)? 0: (new_host_length - 1) * PackedTraits_<T>::VALUES_PER_HOST_VALUE;
            if (new_length <= start + size)
                goto FINISH;
            old_length = new_length - size + old_size;
        }
        else
        {
            new_length = wanted_new_length;
        }
        arrayCopy(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq, end + size - old_size, Standard()));

FINISH:
        _setLength(seq, new_length);
        return size;
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

template<typename TValue, typename THostspec, typename TSize, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type
_clearSpace(String<TValue, Packed<THostspec> > & me,
            TSize size,
            Tag<TExpand>)
{
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size);
}

template<typename TValue, typename THostspec, typename TSize, typename TCapacity, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type
_clearSpace(String<TValue, Packed<THostspec> > & me,
            TSize size,
            TCapacity limit,
            Tag<TExpand>)
{
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size, limit);
}

template<typename TValue, typename THostspec, typename TSize, typename TPosition, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type
_clearSpace(String<TValue, Packed<THostspec> > & me,
            TSize size,
            TPosition pos_begin,
            TPosition pos_end,
            Tag<TExpand>)
{
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end);
}

template<typename TValue, typename THostspec, typename TSize, typename TPosition, typename TCapacity, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type
_clearSpace(String<TValue, Packed<THostspec> > & me,
            TSize size,
            TPosition pos_begin,
            TPosition pos_end,
            TCapacity limit,
            Tag<TExpand>)
{
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end, limit);
}

// --------------------------------------------------------------------------
// Function reserve()
// --------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline typename Size< String<TValue, Packed<TSpec> > >::Type
reserve(
    String<TValue, Packed<TSpec> > & seq,
    TSize_ new_capacity,
    Tag<TExpand> tag)
{
    reserve(host(seq), PackedTraits_<String<TValue, Packed<TSpec> > >::toHostLength1(new_capacity), tag);
    return capacity(seq);
}

// ****************************************************************************
// Functions for Packed String Iter
// ****************************************************************************

// --------------------------------------------------------------------------
// Function hostIterator()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline typename Host<Iter<TPackedString, Packed<THostspec> > >::Type &
hostIterator(Iter<TPackedString, Packed<THostspec> > & me)
{
    return me.data_iterator;
}

template <typename TPackedString, typename THostspec>
inline typename Host<Iter<TPackedString, Packed<THostspec> > const>::Type  &
hostIterator(Iter<TPackedString, Packed<THostspec> > const & me)
{
    return me.data_iterator;
}


// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline typename Reference<Iter<TPackedString, Packed<THostspec> > >::Type
value(Iter<TPackedString, Packed<THostspec> > & me)
{
    return typename Reference<Iter<TPackedString, Packed<THostspec> > >::Type(me);
}

template <typename TPackedString, typename THostspec>
inline typename Reference<Iter<TPackedString, Packed<THostspec> > const>::Type
value(Iter<TPackedString, Packed<THostspec> > const & me)
{
    return typename Reference<Iter<TPackedString, Packed<THostspec> > const>::Type(me);
}

template <typename TPackedString, typename THostspec>
inline typename Reference<Iter<TPackedString const, Packed<THostspec> > >::Type
value(Iter<TPackedString const, Packed<THostspec> > & me)
{
    return getValue(hostIterator(me))[me.localPos];
}

template <typename TPackedString, typename THostspec>
inline typename Reference<Iter<TPackedString const, Packed<THostspec> > const>::Type
value(Iter<TPackedString const, Packed<THostspec> > const & me)
{
    return getValue(hostIterator(me))[me.localPos];
}

// --------------------------------------------------------------------------
// Function getValue()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline typename GetValue<Iter<TPackedString, Packed<THostspec> > >::Type
getValue(Iter<TPackedString, Packed<THostspec> > & me)
{
    return getValue(hostIterator(me))[me.localPos];
}

template <typename TPackedString, typename THostspec>
inline typename GetValue<Iter<TPackedString, Packed<THostspec> > const>::Type
getValue(Iter<TPackedString, Packed<THostspec> > const & me)
{
    return getValue(hostIterator(me))[me.localPos];
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TValue>
inline void
assignValue(Iter<TPackedString, Packed<THostspec> > const & me,
            TValue const & _value)
{
    typedef Iter<TPackedString, Packed<THostspec> > const TIterator;
    assignValue(value(hostIterator(me)), me.localPos, (typename Value<TIterator>::Type)_value);
}

template <typename TPackedString, typename THostspec, typename TValue>
inline void
assignValue(Iter<TPackedString, Packed<THostspec> > & me,
            TValue const & _value)
{
    typedef Iter<TPackedString, Packed<THostspec> > TIterator;
    assignValue(value(hostIterator(me)), me.localPos, (typename Value<TIterator>::Type)_value);
}

// --------------------------------------------------------------------------
// Function moveValue()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TValue>
inline void
moveValue(Iter<TPackedString, Packed<THostspec> > & me,
          TValue const & _value)
{
    assignValue(me, _value);
}

template <typename TPackedString, typename THostspec, typename TValue>
inline void
moveValue(Iter<TPackedString, Packed<THostspec> > const & me,
          TValue const & _value)
{
    assignValue(me, _value);
}

// --------------------------------------------------------------------------
// Function valueConstruct()
// --------------------------------------------------------------------------

//emulate construction and destruction

template <typename TPackedString, typename THostspec>
inline void
valueConstruct(Iter<TPackedString, Packed<THostspec> > const & /*it*/)
{
    // TODO(holtgrew): Why not assign default-constructed?
}

template <typename TPackedString, typename THostspec, typename TParam>
inline void
valueConstruct(Iter<TPackedString, Packed<THostspec> > const & it,
               TParam SEQAN_FORWARD_CARG param_)
{
    assignValue(it, SEQAN_FORWARD(TParam, param_));
}

#ifndef SEQAN_CXX11_STANDARD
template <typename TPackedString, typename THostspec, typename TParam>
inline void
valueConstruct(Iter<TPackedString, Packed<THostspec> > const & it,
               TParam const & param_,
               Move const & /*tag*/)
{
    moveValue(it, param_);
}
#endif

// --------------------------------------------------------------------------
// Function valueDestruct()
// --------------------------------------------------------------------------

// Packed strings cannot contain non-POD data types.

template <typename TPackedString, typename THostspec>
inline void
valueDestruct(Iter<TPackedString, Packed<THostspec> > const & /*it*/)
{
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool
operator==(Iter<TPackedString, Packed<THostspec> > const & left,
           Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) == hostIterator(right)) && (left.localPos == right.localPos);
}

// --------------------------------------------------------------------------
// Function operator!=()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool
operator!=(Iter<TPackedString, Packed<THostspec> > const & left,
           Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) != hostIterator(right)) || (left.localPos != right.localPos);
}

// --------------------------------------------------------------------------
// Function operator>()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool
operator>(Iter<TPackedString, Packed<THostspec> > const & left,
          Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) > hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (left.localPos > right.localPos));
}

// --------------------------------------------------------------------------
// Function operator>=()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool
operator>=(Iter<TPackedString, Packed<THostspec> > const & left,
           Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) > hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (left.localPos >= right.localPos));
}

// --------------------------------------------------------------------------
// Function operator<()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool
operator<(Iter<TPackedString, Packed<THostspec> > const & left,
          Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) < hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (left.localPos < right.localPos));
}

// --------------------------------------------------------------------------
// Function operator<=()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline bool
operator <= (Iter<TPackedString, Packed<THostspec> > const & left,
             Iter<TPackedString, Packed<THostspec> > const & right)
{
    return (hostIterator(left) < hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (left.localPos <= right.localPos));
}

// --------------------------------------------------------------------------
// Function operator++()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline Iter<TPackedString, Packed<THostspec> > &
operator++(Iter<TPackedString, Packed<THostspec> > & me)
{
    if (++me.localPos == PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE)
    {
        me.localPos = 0;
        ++hostIterator(me);
    }
    return me;
}

// --------------------------------------------------------------------------
// Function operator--()
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline Iter<TPackedString, Packed<THostspec> > &
operator--(Iter<TPackedString, Packed<THostspec> > & me)
{
    if (me.localPos-- == 0)
    {
        me.localPos = PackedTraits_<TPackedString>::VALUES_PER_HOST_VALUE - 1;
        --hostIterator(me);
    }
    return me;
}

// --------------------------------------------------------------------------
// Function operator+() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> >
operator+(Iter<TPackedString, Packed<THostspec> > const & iter,
          TIntegral const & delta)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (isNegative(delta))
        return iter - (-(typename MakeSigned<TIntegral>::Type)delta);

    TIntegral ofs = (TIntegral)iter.localPos + delta;
    return Iter<TPackedString, Packed<THostspec> >(
        hostIterator(iter) + ofs / TTraits::VALUES_PER_HOST_VALUE,
        ofs % TTraits::VALUES_PER_HOST_VALUE);
}

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> >
operator+(TIntegral const & delta,
          Iter<TPackedString, Packed<THostspec> > const & iter)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (isNegative(delta))
        return iter - (-(typename MakeSigned<TIntegral>::Type)delta);

    TIntegral ofs = (TIntegral)iter.localPos + delta;
    return Iter<TPackedString, Packed<THostspec> >(
        hostIterator(iter) + ofs / TTraits::VALUES_PER_HOST_VALUE,
        ofs % TTraits::VALUES_PER_HOST_VALUE);
}

// --------------------------------------------------------------------------
// Function operator+=() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> > &
operator+=(Iter<TPackedString, Packed<THostspec> > & iter,
           TIntegral const & delta)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (isNegative(delta))
        return iter -= -(typename MakeSigned<TIntegral>::Type)delta;

    TIntegral ofs = (TIntegral)iter.localPos + delta;
    hostIterator(iter) += ofs / TTraits::VALUES_PER_HOST_VALUE;
    iter.localPos = ofs % TTraits::VALUES_PER_HOST_VALUE;
    return iter;
}

// --------------------------------------------------------------------------
// Function operator-() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> >
operator-(Iter<TPackedString, Packed<THostspec> > const & iter,
          TIntegral const & delta)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (isNegative(delta))
        return iter + (-(typename MakeSigned<TIntegral>::Type)delta);

    TIntegral ofs = delta + (TIntegral)(TTraits::VALUES_PER_HOST_VALUE - 1) - (TIntegral)iter.localPos;
    return Iter<TPackedString, Packed<THostspec> >(
        hostIterator(iter) - ofs / TTraits::VALUES_PER_HOST_VALUE,
        (TTraits::VALUES_PER_HOST_VALUE - 1) - (ofs % TTraits::VALUES_PER_HOST_VALUE));
}

// --------------------------------------------------------------------------
// Function operator-=() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec, typename TIntegral>
inline Iter<TPackedString, Packed<THostspec> > &
operator-=(Iter<TPackedString, Packed<THostspec> > & iter,
           TIntegral const & delta)
{
    typedef PackedTraits_<TPackedString> TTraits;

    if (isNegative(delta))
        return iter += -(typename MakeSigned<TIntegral>::Type)delta;

    TIntegral ofs = delta + (TIntegral)(TTraits::VALUES_PER_HOST_VALUE - 1) - (TIntegral)iter.localPos;
    hostIterator(iter) -= ofs / TTraits::VALUES_PER_HOST_VALUE;
    iter.localPos = (TTraits::VALUES_PER_HOST_VALUE - 1) - (ofs % TTraits::VALUES_PER_HOST_VALUE);
    return iter;
}

// --------------------------------------------------------------------------
// Function operator-() for (iter, iter)
// --------------------------------------------------------------------------

template <typename TPackedString, typename THostspec>
inline typename Difference<Iter<TPackedString, Packed<THostspec> > >::Type
operator-(Iter<TPackedString, Packed<THostspec> > const & iterLeft,
          Iter<TPackedString, Packed<THostspec> > const & iterRight)
{
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename Difference<Iter<TPackedString, Packed<THostspec> > >::Type TDiff;
    return (TDiff)(hostIterator(iterLeft) - hostIterator(iterRight)) * (TDiff)TTraits::VALUES_PER_HOST_VALUE +
           (TDiff)iterLeft.localPos - (TDiff)iterRight.localPos;
}

// ----------------------------------------------------------------------------
// Function bitScanReverse()
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline typename Position<String<bool, Packed<THostSpec> > >::Type
bitScanReverse(String<bool, Packed<THostSpec> > const & obj)
{
    typedef String<bool, Packed<THostSpec> > TPackedString;
    typedef typename Position<TPackedString>::Type TPosition;
    typedef typename Host<TPackedString>::Type TPackedHost;
    typedef typename Iterator<TPackedHost const, Standard>::Type TConstPackedHostIterator;
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename TTraits::THostValue THostValue;  // Is a tuple.
    typedef typename THostValue::TBitVector TBitVector;

    if (empty(host(obj)))
        return MaxValue<TPosition>::VALUE;

    TConstPackedHostIterator it = end(host(obj), Standard()) - 1;
    TConstPackedHostIterator itBegin = begin(host(obj), Standard());

    // We need to treat the last value differently, because it's possible not all bits are in use.
    TBitVector lastVal = it->i & (~static_cast<TBitVector>(0) <<
                                  (BitsPerValue<TBitVector>::VALUE - (length(obj) % BitsPerValue<TBitVector>::VALUE)));

    if (!testAllZeros(lastVal))
        return (((it - itBegin) - 1) * BitsPerValue<TBitVector>::VALUE) +
               (BitsPerValue<TBitVector>::VALUE - 1) - bitScanForward(lastVal);
    --it;
    for(;it != itBegin && testAllZeros(*it); --it)
    {}

    if (it != itBegin)
        return (((it - itBegin) - 1) * BitsPerValue<TBitVector>::VALUE) +
               (BitsPerValue<TBitVector>::VALUE - 1) - bitScanForward(it->i);
    return length(obj);
}

// ----------------------------------------------------------------------------
// Function bitScanForward()
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline typename Position<String<bool, Packed<THostSpec> > >::Type
bitScanForward(String<bool, Packed<THostSpec> > const & obj)
{
    typedef String<bool, Packed<THostSpec> > TPackedString;
    typedef typename Position<TPackedString>::Type TPosition;
    typedef typename Host<TPackedString>::Type TPackedHost;
    typedef typename Iterator<TPackedHost const, Standard>::Type TConstPackedHostIterator;
    typedef PackedTraits_<TPackedString> TTraits;
    typedef typename TTraits::THostValue THostValue;
    typedef typename THostValue::TBitVector TBitVector;

    if (empty(host(obj)))
        return MaxValue<TPosition>::VALUE;

    TConstPackedHostIterator itBegin = begin(host(obj), Standard()) + 1;
    TConstPackedHostIterator it = itBegin;
    TConstPackedHostIterator itEnd = end(host(obj), Standard()) - 1;

    for (; it != itEnd && testAllZeros(*it); ++it)
    {}

    // If last element is not 0, we return the last position. Note, that we do not check for the returned index to be
    // bigger than the length of the string. The caller has to do this.
    TBitVector lastVal = (it != itEnd) ? it->i :
             it->i & (~static_cast<TBitVector>(0) << (BitsPerValue<TBitVector>::VALUE -
                                                      (length(obj) % BitsPerValue<TBitVector>::VALUE)));
    if (testAllZeros(lastVal))
        return length(obj);
    return ((it - itBegin) * BitsPerValue<TBitVector>::VALUE) + (BitsPerValue<TBitVector>::VALUE - 1) - bitScanReverse(lastVal);
}

// ----------------------------------------------------------------------------
// Function transform()
// ----------------------------------------------------------------------------

template <typename THostSpec, typename TBinaryFunctor>
inline void transform(String<bool, Packed<THostSpec> > & target,
                      String<bool, Packed<THostSpec> > const & lhs,
                      String<bool, Packed<THostSpec> > const & rhs,
                      TBinaryFunctor const & binaryFunctor)
{
    typedef String<bool, Packed<THostSpec> > TBitVector;
    typedef typename Host<TBitVector>::Type THost;
    typedef typename Iterator<THost const, Standard>::Type TConstIterator;
    typedef typename Iterator<THost, Standard>::Type TIterator;

    if (empty(host(lhs)) || empty(host(rhs)))
        return;

    SEQAN_ASSERT_EQ(length(lhs), length(rhs));

    resize(target, length(lhs), Exact());

    TConstIterator itLhsBegin = begin(host(lhs), Standard()) + 1;
    TConstIterator itLhsEnd = end(host(lhs), Standard());
    TConstIterator itRhsBegin = begin(host(rhs), Standard()) + 1;
    TIterator itTarget = begin(host(target), Standard()) + 1;

    std::transform(itLhsBegin, itLhsEnd, itRhsBegin, itTarget, binaryFunctor);
}

// ----------------------------------------------------------------------------
// Function transform()
// ----------------------------------------------------------------------------

template <typename THostSpec, typename TUnaryFunctor>
inline void transform(String<bool, Packed<THostSpec> > & target,
                      String<bool, Packed<THostSpec> > const & src,
                      TUnaryFunctor const & unaryFunctor)
{
    typedef String<bool, Packed<THostSpec> > TBitVector;
    typedef typename Host<TBitVector>::Type THost;
    typedef typename Iterator<THost const, Standard>::Type TConstIterator;
    typedef typename Iterator<THost, Standard>::Type TIterator;

    if (empty(host(src)))
        return;

    resize(target, length(src), Exact());

    TConstIterator itSrcBegin = begin(host(src), Standard()) + 1;
    TConstIterator itSrcEnd = end(host(src), Standard());
    TIterator itTarget = begin(host(target), Standard()) + 1;

    std::transform(itSrcBegin, itSrcEnd, itTarget, unaryFunctor);
}

// ----------------------------------------------------------------------------
// Function operator|
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline String<bool, Packed<THostSpec> >
operator|(String<bool, Packed<THostSpec> > const & lhs, String<bool, Packed<THostSpec> > const & rhs)
{
    String<bool, Packed<THostSpec> > tmp;
    transform(tmp, lhs, rhs, FunctorBitwiseOr());
    return tmp;
}

// ----------------------------------------------------------------------------
// Function operator|=
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline String<bool, Packed<THostSpec> > &
operator|=(String<bool, Packed<THostSpec> > & lhs, String<bool, Packed<THostSpec> > const & rhs)
{
    transform(lhs, lhs, rhs, FunctorBitwiseOr());
    return lhs;
}

// ----------------------------------------------------------------------------
// Function operator&
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline String<bool, Packed<THostSpec> >
operator&(String<bool, Packed<THostSpec> > const & lhs, String<bool, Packed<THostSpec> > const & rhs)
{
    String<bool, Packed<THostSpec> > tmp;
    transform(tmp, lhs, rhs, FunctorBitwiseAnd());
    return tmp;
}

// ----------------------------------------------------------------------------
// Function operator&=
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline String<bool, Packed<THostSpec> > &
operator&=(String<bool, Packed<THostSpec> > & lhs, String<bool, Packed<THostSpec> > const & rhs)
{
    transform(lhs, lhs, rhs, FunctorBitwiseAnd());
    return lhs;
}

// ----------------------------------------------------------------------------
// Function operator~
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline String<bool, Packed<THostSpec> >
operator~(String<bool, Packed<THostSpec> > const & vec)
{
    String<bool, Packed<THostSpec> > tmp;
    transform(tmp, vec, FunctorBitwiseNot());
    return tmp;
}

// ----------------------------------------------------------------------------
// Function operator^
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline String<bool, Packed<THostSpec> >
operator^(String<bool, Packed<THostSpec> > const & lhs, String<bool, Packed<THostSpec> > const & rhs)
{
    String<bool, Packed<THostSpec> > tmp;
    transform(tmp, lhs, rhs, FunctorBitwiseXor());
    return tmp;
}

// ----------------------------------------------------------------------------
// Function operator^=
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline String<bool, Packed<THostSpec> > &
operator^=(String<bool, Packed<THostSpec> > & lhs, String<bool, Packed<THostSpec> > const & rhs)
{
    transform(lhs, lhs, rhs, FunctorBitwiseXor());
    return lhs;
}

// ----------------------------------------------------------------------------
// Function _packedStringTestAll()
// ----------------------------------------------------------------------------

template <typename THostSpec, typename TTester>
inline bool
_packedStringTestAll(String<bool, Packed<THostSpec> > const & obj, TTester const & tester)
{
    typedef String<bool, Packed<THostSpec> > TPackedString;
    typedef typename Host<TPackedString>::Type TPackedHost;
    typedef typename Iterator<TPackedHost const>::Type TConstIterator;

    if (empty(host(obj)))
        return false;

    TConstIterator itFirst = begin(host(obj), Standard()) + 1;
    TConstIterator itLast = end(host(obj), Standard()) - 1;

    while (itFirst != itLast && tester(*itFirst, False()))
        ++itFirst;

    if (itFirst != itLast)
        return false;

    return tester(*(itFirst), True());  // Test last case.
}

// ----------------------------------------------------------------------------
// Function testAllZeros()
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline bool
testAllZeros(String<bool, Packed<THostSpec> > const & obj)
{
    typedef String<bool, Packed<THostSpec> > TPackedString;
    typedef PackedTraits_<TPackedString> TTraits;

    return _packedStringTestAll(obj, FunctorTestAllZeros<TPackedString>((TTraits::VALUES_PER_HOST_VALUE -
                                     (length(obj) % TTraits::VALUES_PER_HOST_VALUE))));
}

// ----------------------------------------------------------------------------
// Function testAllOnes()
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline bool
testAllOnes(String<bool, Packed<THostSpec> > const & obj)
{
    typedef String<bool, Packed<THostSpec> > TPackedString;
    typedef PackedTraits_<TPackedString> TTraits;

    return _packedStringTestAll(obj, FunctorTestAllOnes<TPackedString>((TTraits::VALUES_PER_HOST_VALUE -
                                     (length(obj) % TTraits::VALUES_PER_HOST_VALUE))));
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline bool open(String<TValue, Packed<THostspec> > & me, const char *fileName, int openMode)
{
    return open(host(me), fileName, openMode);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline bool save(String<TValue, Packed<THostspec> > const & me, const char *fileName, int openMode)
{
    // the visible part of the string is kept untouched and the function is thread-safe
    _clearUnusedBits(const_cast<String<TValue, Packed<THostspec> > &>(me));
    return save(host(me), fileName, openMode);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_PACKED_H_
