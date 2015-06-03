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
// Implementation of the Suffix Segment specialization.
// ==========================================================================

#ifndef SEQAN_HEADER_SEGMENT_SUFFIX_H
#define SEQAN_HEADER_SEGMENT_SUFFIX_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// SuffixSegment
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class SuffixSegment Suffix Segment
 * @extends Segment
 * @headerfile <seqan/sequence.h>
 * @brief A suffix of a sequence.
 *
 * @signature template <typename THost>
 *            class Segment<THost, SuffixSegment>;
 *
 * @tparam THost The underlying @link ContainerConcept sequence @endlink type.
 *
 * @see InfixSegment
 * @see PrefixSegment
 * @aka substring
 */

struct SuffixSegment {};

template <typename THost_>
class Segment<THost_, SuffixSegment>
{
public:
    typedef typename Host<Segment>::Type THost;

    typename Pointer_<THost>::Type data_host;
    typename Position<THost>::Type data_begin_position;

//____________________________________________________________________________

public:

    SEQAN_HOST_DEVICE
    Segment():
        data_host(),
        data_begin_position(0)
    {
SEQAN_CHECKPOINT
    }

    Segment(THost & _host):
        data_host(& _host),
        data_begin_position(0)
    {
SEQAN_CHECKPOINT
    }

    Segment(typename Parameter_<THost>::Type _host, typename Position<THost>::Type _begin_index):
        data_host(_toPointer(_host)),
        data_begin_position(_begin_index)
    {
SEQAN_CHECKPOINT
    }
/*
    Segment(typename Parameter_<THost>::Type _host, typename Iterator<THost, Rooted>::Type _begin):
        data_host(_toPointer(_host)),
        data_begin_position(position(_begin))
    {
SEQAN_CHECKPOINT
    }
*/
    Segment(typename Parameter_<THost>::Type _host, typename Iterator<THost, Standard>::Type _begin):
        data_host(_toPointer(_host)),
        data_begin_position(position(_begin, _host))
    {
SEQAN_CHECKPOINT
    }
/*
    Segment(Segment const & _other):
        data_host(_other.data_host),
        data_begin_position(_other.data_begin_position)
    {
SEQAN_CHECKPOINT
    }
*/
    template <typename THost2, typename TSpec2>
    Segment(Segment<THost2, TSpec2> const & _other):
        data_host(_toPointer(host(_other))),
        data_begin_position(beginPosition(_other))
    {
SEQAN_CHECKPOINT
    }

    SEQAN_HOST_DEVICE
    ~ Segment()
    {
SEQAN_CHECKPOINT
    }

    inline Segment &
    operator = (Segment const & source)
    {
        assign(*this, source);
        return *this;
    }
//____________________________________________________________________________

public:


//____________________________________________________________________________

    template <typename TPos>
    inline typename Reference<Segment>::Type
    operator [] (TPos pos)
    {
SEQAN_CHECKPOINT
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<Segment const>::Type
    operator [] (TPos pos) const
    {
SEQAN_CHECKPOINT
        return value(*this, pos);
    }

//____________________________________________________________________________
};

// template <typename THost>
// inline void
// clear(Segment<THost, SuffixSegment> & target)
// {
//     replace(host(target), beginPosition(target), endPosition(target), "");
// }

//////////////////////////////////////////////////////////////////////////////

template <typename THost_>
SEQAN_HOST_DEVICE inline typename Parameter_<THost_>::Type
host(Segment<THost_, SuffixSegment> & me)
{
SEQAN_CHECKPOINT
    return _toParameter<THost_>(me.data_host);
}

template <typename THost_>
SEQAN_HOST_DEVICE inline typename Parameter_<THost_>::Type
host(Segment<THost_, SuffixSegment> const & me)
{
SEQAN_CHECKPOINT
    return _toParameter<THost_>(me.data_host);
}

//____________________________________________________________________________

template <typename THost_>
inline void
setHost(Segment<THost_, SuffixSegment> & me, typename Parameter_<THost_>::Type _host)
{
SEQAN_CHECKPOINT
    me.data_host = _toPointer(_host);
}

template <typename THost_>
inline void
setHost(Segment<THost_ const, SuffixSegment> & me, typename Parameter_<THost_>::Type _host)
{
SEQAN_CHECKPOINT
    me.data_host = _toPointer(_host);
}

//____________________________________________________________________________

template <typename THost_>
SEQAN_HOST_DEVICE inline typename Iterator<Segment<THost_, SuffixSegment>, Standard>::Type
begin(Segment<THost_, SuffixSegment> & me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(host(me), Standard()) + me.data_begin_position;
}
template <typename THost_>
SEQAN_HOST_DEVICE inline typename Iterator<Segment<THost_, SuffixSegment> const, Standard>::Type
begin(Segment<THost_, SuffixSegment> const & me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(host(me), Standard()) + me.data_begin_position;
}

//____________________________________________________________________________

template <typename THost_>
inline typename Position<Segment<THost_, SuffixSegment> const>::Type
beginPosition(Segment<THost_, SuffixSegment> const & me)
{
SEQAN_CHECKPOINT
    return me.data_begin_position;
}
template <typename THost_>
inline typename Position<Segment<THost_, SuffixSegment> >::Type
beginPosition(Segment<THost_, SuffixSegment> & me)
{
SEQAN_CHECKPOINT
    return me.data_begin_position;
}
//____________________________________________________________________________

template <typename THost_, typename TIterator>
inline void
setBegin(Segment<THost_, SuffixSegment> & me, TIterator new_begin)
{
SEQAN_CHECKPOINT
    me.data_begin_position = new_begin - begin(host(me));//, Standard());
}

template <typename THost_>
inline void
setBegin(typename Iterator<Segment<THost_, SuffixSegment>, Rooted>::Type new_begin)
{
SEQAN_CHECKPOINT
    container(new_begin).data_begin_position = hostIterator(new_begin) - begin(host(container(new_begin)));//, Standard());
}

//____________________________________________________________________________

template <typename THost_, typename TPosition>
inline void
setBeginPosition(Segment<THost_, SuffixSegment> & me, TPosition new_begin)
{
SEQAN_CHECKPOINT
    me.data_begin_position = new_begin;
}

//____________________________________________________________________________

template <typename THost_>
inline typename Iterator<Segment<THost_, SuffixSegment>, Standard>::Type
end(Segment<THost_, SuffixSegment> & me,
    Standard)
{
SEQAN_CHECKPOINT
    return end(host(me), Standard());
}
template <typename THost_>
inline typename Iterator<Segment<THost_, SuffixSegment> const, Standard>::Type
end(Segment<THost_, SuffixSegment> const & me,
    Standard)
{
SEQAN_CHECKPOINT
    return end(host(me), Standard());
}

//____________________________________________________________________________


template <typename THost_>
inline typename Position<Segment<THost_, SuffixSegment> >::Type
endPosition(Segment<THost_, SuffixSegment> & me)
{
SEQAN_CHECKPOINT
    return length(host(me));
}

template <typename THost_>
inline typename Position<Segment<THost_, SuffixSegment> const>::Type
endPosition(Segment<THost_, SuffixSegment> const & me)
{
SEQAN_CHECKPOINT
    return length(host(me));
}

//____________________________________________________________________________

template <typename TIterator, typename THost_>
inline void
setEnd(Segment<THost_, SuffixSegment> &, TIterator)
{
}

template <typename THost_, typename TSize>
inline void
_setLength(Segment<THost_, SuffixSegment> &, TSize)
{
}

//////////////////////////////////////////////////////////////////////////////

struct PrefixSegment;
struct InfixSegment;

template <typename THost>
struct Suffix
{
    typedef Segment<THost, SuffixSegment> Type;
};

template <typename THost>
struct Suffix< Segment<THost, InfixSegment> >
{
    typedef Segment<THost, InfixSegment> Type;
};
template <typename THost>
struct Suffix< Segment<THost, SuffixSegment> >
{
    typedef Segment<THost, SuffixSegment> Type;
};
template <typename THost>
struct Suffix< Segment<THost, PrefixSegment> >
{
    typedef Segment<THost, InfixSegment> Type;
};

template <typename THost, typename TSpec>
struct Suffix< Segment<THost, TSpec> const >:
    Suffix< Segment<THost, TSpec> > {};

template <typename THost>
struct Suffix<THost &>:
    Suffix<THost> {};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TPosition>
inline void
set(Segment<THost, SuffixSegment> & me,
    THost & host_,
    TPosition begin_)
{
SEQAN_CHECKPOINT
    setHost(me, host_);
    setBeginPosition(me, begin_);
}
//____________________________________________________________________________

template <typename THost>
inline void
set(Segment<THost, SuffixSegment> & me,
    THost & host_)
{
SEQAN_CHECKPOINT
    setHost(me, host_);
    setBegin(me, begin(host_, Standard()));
}

//____________________________________________________________________________

template <typename THost, typename TSpec>
inline void
set(Segment<THost, SuffixSegment> & me,
    Segment<THost, TSpec> & source)
{
SEQAN_CHECKPOINT
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
}

template <typename THost, typename TSpec>
inline void
set(Segment<THost, SuffixSegment> & me,
    Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atBegin(Segment<THost, SuffixSegment> const & segment)
{
SEQAN_CHECKPOINT
    return (beginPosition(segment) == 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atEnd(Segment<THost, SuffixSegment> const & segment)
{
SEQAN_CHECKPOINT
    return (beginPosition(segment) == length(host(segment)));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goBegin(Segment<THost, SuffixSegment> & segment,
        THost &)
{
SEQAN_CHECKPOINT
    goBegin(segment);
}

template <typename THost>
inline void
goBegin(Segment<THost, SuffixSegment> & segment)
{
    setBegin(segment);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goEnd(Segment<THost, SuffixSegment> & segment,
      THost &)
{
SEQAN_CHECKPOINT
    goEnd(segment);
}

template <typename THost>
inline void
goEnd(Segment<THost, SuffixSegment> & segment)
{
    setBegin(segment, length(host(segment))-1);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, SuffixSegment> &
operator ++(Segment<THost, SuffixSegment> & segment)
{
    setBegin(segment, beginPosition(segment) + 1);
    return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, SuffixSegment> &
operator --(Segment<THost, SuffixSegment> & segment)
{
    setBegin(segment, beginPosition(segment) - 1);
    return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename TPosBegin>
inline typename Suffix<T>::Type
suffix(T & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
    return typename Suffix<T>::Type(t, pos_begin);
}
template <typename T, typename TPosBegin>
inline typename Suffix<T const>::Type
suffix(T const & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
    return typename Suffix<T const>::Type(t, pos_begin);
}

template <typename T, typename TPosBegin>
inline typename Suffix<T *>::Type
suffix(T * t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
    return typename Suffix<T *>::Type (t, pos_begin);
}

//////////////////////////////////////////////////////////////////////////////
// A suffix of a prefix -> is an infix
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, PrefixSegment> >::Type
suffix(Segment<T, PrefixSegment> & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, PrefixSegment> >::Type (
        host(t),
        beginPosition(t) + pos_begin,
        endPosition(t));
}
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, PrefixSegment> const>::Type
suffix(Segment<T, PrefixSegment> const & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, PrefixSegment> const>::Type (
        host(t),
        beginPosition(t) + pos_begin,
        endPosition(t));
}

//////////////////////////////////////////////////////////////////////////////
// A suffix of a infix -> is an infix
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, InfixSegment> >::Type
suffix(Segment<T, InfixSegment> & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, InfixSegment> >::Type (
        host(t),
        beginPosition(t) + pos_begin,
        endPosition(t));
}
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, InfixSegment> const>::Type
suffix(Segment<T, InfixSegment> const & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, InfixSegment> const>::Type (
        host(t),
        beginPosition(t) + pos_begin,
        endPosition(t));
}


//////////////////////////////////////////////////////////////////////////////
// A suffix of a suffix -> is a suffix
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, SuffixSegment> >::Type
suffix(Segment<T, SuffixSegment> & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, SuffixSegment> >::Type (
        host(t),
        beginPosition(t) + pos_begin);
}
template <typename T, typename TPosBegin>
inline typename Suffix<Segment<T, SuffixSegment> const>::Type
suffix(Segment<T, SuffixSegment> const & t, TPosBegin pos_begin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, SuffixSegment> const>::Type (
        host(t),
        beginPosition(t) + pos_begin);
}

//////////////////////////////////////////////////////////////////////////////

// suffix() with iterators

template <typename T, typename TIterSpec>
inline typename Suffix<T>::Type
suffix(T & t,
       Iter<Segment<T, PrefixSegment>, TIterSpec> const & iterBegin)
{
SEQAN_CHECKPOINT
    return typename Suffix<T>::Type(t, iterBegin);
}
template <typename T, typename TIterSpec>
inline typename Suffix<T const>::Type
suffix(T const & t,
       Iter<Segment<T, PrefixSegment> const, TIterSpec> const & iterBegin)
{
SEQAN_CHECKPOINT
    return typename Suffix<T const>::Type(t, iterBegin);
}

template <typename T>
inline typename Suffix<T *>::Type
suffix(T * t,
       T * & iterBegin)
{
SEQAN_CHECKPOINT
    return typename Suffix<T *>::Type (t, iterBegin);
}

//////////////////////////////////////////////////////////////////////////////
// A suffix of a prefix -> is an infix
template <typename T, typename TIterSpec>
inline typename Suffix<Segment<T, PrefixSegment> >::Type
suffix(Segment<T, PrefixSegment> & t,
       Iter<Segment<T, PrefixSegment>, TIterSpec> const & iterBegin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, PrefixSegment> >::Type (
        host(t),
        iterBegin,
        end(t));
}
template <typename T, typename TIterSpec>
inline typename Suffix<Segment<T, PrefixSegment> const>::Type
suffix(Segment<T, PrefixSegment> const & t,
       Iter<Segment<T, PrefixSegment> const, TIterSpec> const & iterBegin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, PrefixSegment> const>::Type (
        host(t),
        iterBegin,
        end(t));
}

//////////////////////////////////////////////////////////////////////////////
// A suffix of a infix -> is an infix
template <typename T, typename TIterSpec>
inline typename Suffix<Segment<T, InfixSegment> >::Type
suffix(Segment<T, InfixSegment> & t,
       Iter<Segment<T, InfixSegment>, TIterSpec> const & iterBegin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, InfixSegment> >::Type (
        host(t),
        iterBegin,
        end(t));
}
template <typename T, typename TIterSpec>
inline typename Suffix<Segment<T, InfixSegment> const>::Type
suffix(Segment<T, InfixSegment> const & t,
       Iter<Segment<T, InfixSegment> const, TIterSpec> const & iterBegin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, InfixSegment> const>::Type (
        host(t),
        iterBegin,
        end(t));
}


//////////////////////////////////////////////////////////////////////////////
// A suffix of a suffix -> is a suffix
template <typename T, typename TIterSpec>
inline typename Suffix<Segment<T, SuffixSegment> >::Type
suffix(Segment<T, SuffixSegment> & t,
       Iter<Segment<T, SuffixSegment>, TIterSpec> const & iterBegin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, SuffixSegment> >::Type (
        host(t),
        iterBegin);
}
template <typename T, typename TIterSpec>
inline typename Suffix<Segment<T, SuffixSegment> const>::Type
suffix(Segment<T, SuffixSegment> const & t,
       Iter<Segment<T, SuffixSegment> const, TIterSpec> const & iterBegin)
{
SEQAN_CHECKPOINT
    return typename Suffix<Segment<T, SuffixSegment> const>::Type (
        host(t),
        iterBegin);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
