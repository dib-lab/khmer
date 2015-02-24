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
// Implementation of the specialization ConcatVirtual of class Iter that
// allows the iteration of arbitrary StringSet objects as if they were
// the concatenation of all strings.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_ITER_CONCAT_VIRTUAL_H_
#define SEQAN_SEQUENCE_ITER_CONCAT_VIRTUAL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// The Metafunction Concatenator is actually defined in string_set_base.h.
template <typename T>
struct Concatenator;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TDelimiter = void >
struct ConcatVirtual;

/*!
 * @class ConcatVirtualIterator
 * @extends Iter
 * @headerfile <seqan/sequence.h>
 * @brief Iterator that sequentially iterates through the elements of a @link StringSet @endlink as if they were
 *        directly concatenated, also see @link ConcatDirectStringSet @endlink.
 *
 * @signature template <typename TStringSet[, typename TDelimiter]>
 *            class Iter<TStringSet, ConcatVirtual<TDelimiter> >;
 *
 * @tparam TStringSet The @link StringSet @endlink to iterate.
 * @tparam TDelimiter The delimiter character.  Default: <tt>void</tt>.
 *
 *
 * @fn ConcatVirtualIterator::Iter
 * @brief Constructor.
 *
 * @signature Iter::Iter(host[, objNo, offset]);
 *
 * @param[in] host   Container to iterate.
 * @param[in] objNo  Sequence number to set the iterator to.  Default: 0.
 * @param[in] offset Offset in the object (specified by objNo) to point to.  Defaults 0.
 *
 * If <tt>objNo</tt> and <tt>offset</tt> are not given then the iterator will point to the first element with offset 0.
 */

template <typename TStringSet, typename TSpec >
class Iter<TStringSet, ConcatVirtual<TSpec> >
{
public:
    typedef typename Value<TStringSet>::Type    TString;
    typedef typename Value<TString>::Type       TValue;
    typedef typename Size<TString>::Type        TSize;

    // TODO(holtgrew): obj_iterator and const_obj_iterator do not appear to be in any C++/STL standard and do not conform to SeqAn's naming scheme.
    typedef typename Iterator<TString, Standard>::Type            obj_iterator;
    typedef typename Iterator<TString const, Standard>::Type      const_obj_iterator;

    // ----------------------------------------------------------------------
    // STL compatible public iterator interface
    typedef Iter                                iterator;
    typedef std::bidirectional_iterator_tag   iterator_category;
    typedef TValue                              value_type;
    typedef TValue &                            reference;
    typedef TValue const &                      const_reference;
    typedef TValue *                            pointer;
    typedef TSize                               size_type;
    typedef typename Difference<TString>::Type  difference_type;
    // ----------------------------------------------------------------------

    TStringSet *    host;
    unsigned        objNo;
    obj_iterator    _begin, _cur, _end;

    Iter() {}

    Iter(TStringSet &_host):
        host(&_host),
        objNo(0)
    {
        if (!empty(_host))
        {
            _begin = _cur = begin(_host[0]);
            _end = end(_host[0]);
            _testEnd();
        }
        else
        {
            _begin = _cur = _end = obj_iterator();
        }
    }

    Iter(TStringSet &_host, unsigned _objNo, difference_type _offset):
        host(&_host),
        objNo(0)
    {
        _setPosition(_objNo, _offset);
    }

    // ----------------------------------------------------------------------
    // Conversion operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    inline operator obj_iterator()
    {
        return _cur;
    }

    inline operator void * ()
    {
        return _cur;
    }

    // ----------------------------------------------------------------------
    // Shared methods for methods that *have* to be defined in the class def.
    // ----------------------------------------------------------------------

    // _testEnd() is used in the constructors only, so having it here as a
    // non-public function is allowed.
    inline void _testEnd()
    {
        if (objNo >= length(*host))
        {
            SEQAN_ASSERT(_begin == obj_iterator());
            SEQAN_ASSERT(_cur == obj_iterator());
            SEQAN_ASSERT(_end == obj_iterator());
            return;
        }

        while (true)
        {
            if (_cur != _end)
                return;
            if (++objNo == length(*host))
                break;
            _begin = _cur = begin((*host)[objNo]);
            _end = end((*host)[objNo]);
        }

        _begin = _cur = _end = obj_iterator();
    }

    inline void _setPosition(unsigned _objNo, difference_type _offset)
    {
        objNo = _objNo;
        if (_objNo < length(*host))
        {
            _begin = _cur = begin((*host)[objNo]);
            _end = end((*host)[objNo]);
            goFurther(_cur, _offset);
            _testEnd();
        }
        else
        {
            _begin = _cur = _end = obj_iterator();
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction Value
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec>
struct Value<Iter<TStringSet, ConcatVirtual<TSpec> > >
    : Value<typename Value<TStringSet>::Type> {};

template <typename TStringSet, typename TSpec>
struct Value<Iter<TStringSet, ConcatVirtual<TSpec> > const>
    : Value<typename Value<TStringSet>::Type> {};

// --------------------------------------------------------------------------
// Metafunction GetValue
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec>
struct GetValue<Iter<TStringSet, ConcatVirtual<TSpec> > >
    : GetValue<typename Value<TStringSet>::Type> {};

template <typename TStringSet, typename TSpec>
struct GetValue<Iter<TStringSet, ConcatVirtual<TSpec> > const>
    : GetValue<typename Value<TStringSet>::Type> {};

// --------------------------------------------------------------------------
// Metafunction Size
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec>
struct Size<Iter<TStringSet, ConcatVirtual<TSpec> > >
    : Size<typename Value<TStringSet>::Type> {};
// Default implementation Size<T const> redirects to non-const variant.

// --------------------------------------------------------------------------
// Metafunction Reference
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec>
struct Reference<Iter<TStringSet, ConcatVirtual<TSpec> > >
    : Reference<typename Value<TStringSet>::Type> {};

template <typename TStringSet, typename TSpec>
struct Reference<Iter<TStringSet, ConcatVirtual<TSpec> > const >
    : Reference<typename Value<TStringSet>::Type> {};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Functions value(), operator*()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec>
inline typename Reference<Iter<TStringSet, ConcatVirtual<TSpec> > const>::Type
value(Iter<TStringSet, ConcatVirtual<TSpec> > const & me)
{
    return *me._cur;
}

template <typename TStringSet, typename TSpec>
inline typename Reference<Iter<TStringSet, ConcatVirtual<TSpec> > >::Type
value(Iter<TStringSet, ConcatVirtual<TSpec> > & me)
{
    return *me._cur;
}

template <typename TStringSet, typename TSpec>
inline typename Reference<Iter<TStringSet, ConcatVirtual<TSpec> > const>::Type
operator*(Iter<TStringSet, ConcatVirtual<TSpec> > const & me)
{
    return *me._cur;
}

template <typename TStringSet, typename TSpec>
inline typename Reference<Iter<TStringSet, ConcatVirtual<TSpec> > >::Type
operator*(Iter<TStringSet, ConcatVirtual<TSpec> > & me)
{
    return *me._cur;
}

// --------------------------------------------------------------------------
// Functions goNext(), operator++()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec>
inline void
goNext(Iter<TStringSet, ConcatVirtual<TSpec> > & me)
{
    ++me._cur;
    me._testEnd();
}

template <typename TStringSet, typename TSpec>
inline Iter<TStringSet, ConcatVirtual<TSpec> > const &
operator++(Iter<TStringSet, ConcatVirtual<TSpec> > & me)
{
    goNext(me);
    return me;
}

template <typename TStringSet, typename TSpec>
inline Iter<TStringSet, ConcatVirtual<TSpec> >
operator++(Iter<TStringSet, ConcatVirtual<TSpec> > & me, int)
{
    Iter<TStringSet, ConcatVirtual<TSpec> > before = me;
    goNext(me);
    return before;
}

// --------------------------------------------------------------------------
// Functions goPrevious(), operator++()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec>
inline void
goPrevious(Iter<TStringSet, ConcatVirtual<TSpec> > & me)
{
    while (me._cur == me._begin && me.objNo > 0)
    {
        --me.objNo;
        me._begin = begin((*me.host)[me.objNo]);
        me._end = me._cur = end((*me.host)[me.objNo]);
    }
    --me._cur;
}

template <typename TStringSet, typename TSpec>
inline Iter<TStringSet, ConcatVirtual<TSpec> > const &
operator--(Iter<TStringSet, ConcatVirtual<TSpec> > & me)
{
    goPrevious(me);
    return me;
}

template <typename TStringSet, typename TSpec>
inline Iter<TStringSet, ConcatVirtual<TSpec> >
operator--(Iter<TStringSet, ConcatVirtual<TSpec> > & me, int)
{
    Iter<TStringSet, ConcatVirtual<TSpec> > before = me;
    goPrevious(me);
    return before;
}

// --------------------------------------------------------------------------
// Helper function _tell()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec>
inline typename Size<typename Value<TStringSet>::Type >::Type
_tell(Iter<TStringSet, ConcatVirtual<TSpec> > const & me)
{
    typedef typename Size<typename Value<TStringSet>::Type >::Type TSize;
    typedef Pair<unsigned, TSize> TPair;
    return posGlobalize(TPair(me.objNo, difference(me._begin, me._cur)), stringSetLimits(*me.host));
}

// --------------------------------------------------------------------------
// Function operator+()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec, typename TDelta>
inline Iter<TStringSet, ConcatVirtual<TSpec> >
operator+(Iter<TStringSet, ConcatVirtual<TSpec> > const & me, TDelta delta)
{
    Pair<unsigned, typename Size<typename Value<TStringSet>::Type>::Type> pos;
    posLocalize(pos, _tell(me) + delta, stringSetLimits(*me.host));
    return Iter<TStringSet, ConcatVirtual<TSpec> > (*me.host, getValueI1(pos), getValueI2(pos));
}

template <typename TStringSet, typename TSpec, typename T1, typename T2, typename TPack>
inline Iter<TStringSet, ConcatVirtual<TSpec> >
operator+(Iter<TStringSet, ConcatVirtual<TSpec> > const & me, Pair<T1, T2, TPack> delta)
{
    Pair<unsigned, typename Size<typename Value<TStringSet>::Type>::Type> pos;
    posLocalize(pos, _tell(me) + delta, stringSetLimits(*me.host));
    return Iter<TStringSet, ConcatVirtual<TSpec> > (*me.host, getValueI1(pos), getValueI2(pos));
}

// --------------------------------------------------------------------------
// Function operator+=()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec, typename TDelta>
inline Iter<TStringSet, ConcatVirtual<TSpec> >
operator+=(Iter<TStringSet, ConcatVirtual<TSpec> > & me, TDelta delta)
{
    Pair<unsigned, typename Size<typename Value<TStringSet>::Type>::Type> pos;
    posLocalize(pos, _tell(me) + delta, stringSetLimits(*me.host));
    me._setPosition(getValueI1(pos), getValueI2(pos));
    return me;
}

// --------------------------------------------------------------------------
// Function operator-()
// --------------------------------------------------------------------------

template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
typename Difference<Iter<TSSetL, ConcatVirtual<TSpecL> > >::Type
operator-(
    Iter<TSSetL, ConcatVirtual<TSpecL> > const & L,
    Iter<TSSetR, ConcatVirtual<TSpecR> > const & R)
{
    return _tell(L) - _tell(R);
}

template <typename TStringSet, typename TSpec, typename TDelta>
inline Iter<TStringSet, ConcatVirtual<TSpec> >
operator-(Iter<TStringSet, ConcatVirtual<TSpec> > const & me, TDelta delta)
{
    Pair<unsigned, typename Size<typename Value<TStringSet>::Type>::Type> pos;
    posLocalize(pos, _tell(me) - delta, stringSetLimits(*me.host));
    return Iter<TStringSet, ConcatVirtual<TSpec> > (*me.host, getValueI1(pos), getValueI2(pos));
}

// --------------------------------------------------------------------------
// Function operator-=()
// --------------------------------------------------------------------------

template <typename TStringSet, typename TSpec, typename TDelta>
inline Iter<TStringSet, ConcatVirtual<TSpec> >
operator-=(Iter<TStringSet, ConcatVirtual<TSpec> > & me, TDelta delta)
{
    Pair<unsigned, typename Size<typename Value<TStringSet>::Type>::Type> pos;
    posLocalize(pos, _tell(me) - delta, stringSetLimits(*me.host));
    me._setPosition(getValueI1(pos), getValueI2(pos));
    return me;
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
inline bool
operator==(
    Iter<TSSetL, ConcatVirtual<TSpecL> > const & L,
    Iter<TSSetR, ConcatVirtual<TSpecR> > const & R)
{
    SEQAN_ASSERT_EQ(L.host, R.host);
    return L.objNo == R.objNo && L._cur == R._cur;
}

template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
inline bool
operator!=(
    Iter<TSSetL, ConcatVirtual<TSpecL> > const & L,
    Iter<TSSetR, ConcatVirtual<TSpecR> > const & R)
{
    SEQAN_ASSERT(L.host == R.host);
    return L.objNo != R.objNo || L._cur != R._cur;
}

// --------------------------------------------------------------------------
// Function operator<()
// --------------------------------------------------------------------------

template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
inline bool
operator<(
    Iter<TSSetL, ConcatVirtual<TSpecL> > const & L,
    Iter<TSSetR, ConcatVirtual<TSpecR> > const & R)
{
    SEQAN_ASSERT(L.host == R.host);
    return L.objNo < R.objNo || (L.objNo == R.objNo && L._cur < R._cur);
}

template <typename TSSetL, typename TSpecL, typename TSSetR, typename TSpecR>
inline bool
operator > (
    Iter<TSSetL, ConcatVirtual<TSpecL> > const & L,
    Iter<TSSetR, ConcatVirtual<TSpecR> > const & R)
{
    SEQAN_ASSERT(L.host == R.host);
    return L.objNo > R.objNo || (L.objNo == R.objNo && L._cur > R._cur);
}

// --------------------------------------------------------------------------
// Function container()
// --------------------------------------------------------------------------

template <typename TSSet, typename TSpec>
inline typename Concatenator<TSSet>::Type &
container(Iter<TSSet, ConcatVirtual<TSpec> > & me)
{
    return concat(*me.host);
}

template <typename TSSet, typename TSpec>
inline typename Concatenator<TSSet>::Type &
container(Iter<TSSet, ConcatVirtual<TSpec> > const & me)
{
    return concat(*me.host);
}

// --------------------------------------------------------------------------
// Function atBegin()
// --------------------------------------------------------------------------

template <typename TSSet, typename TSpec>
inline bool
atBegin(Iter<TSSet, ConcatVirtual<TSpec> > & me)
{
    return me._cur == me._begin && me.objNo == 0;
}

template <typename TSSet, typename TSpec>
inline bool
atBegin(Iter<TSSet, ConcatVirtual<TSpec> > const & me)
{
    return me._cur == me._begin && me.objNo == 0;
}

// --------------------------------------------------------------------------
// Function atEnd()
// --------------------------------------------------------------------------

template <typename TSSet, typename TSpec>
inline bool
atEnd(Iter<TSSet, ConcatVirtual<TSpec> > & me)
{
    return me._cur == me._end;
}

template <typename TSSet, typename TSpec>
inline bool
atEnd(Iter<TSSet, ConcatVirtual<TSpec> > const & me)
{
    return me._cur == me._end;
}

// --------------------------------------------------------------------------
// Function atEndOfSequence()
// --------------------------------------------------------------------------

// TODO(holtgrew): Specifying a catch-all implementation appears a bit too generous, what about concept checking?

template <typename TIterator>
inline bool
atEndOfSequence(TIterator const & me)
{
    return atEnd(me);
}

template <typename TSSet, typename TSpec>
inline bool
atEndOfSequence(Iter<TSSet, ConcatVirtual<TSpec> > const & me)
{
    return me._cur == me._begin && me.objNo > 0;
}

template <typename TIterator>
inline bool
atEndOfSequence(TIterator & me)
{
    return atEndOfSequence(reinterpret_cast<TIterator const &>(me));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_ITER_CONCAT_VIRTUAL_H_
