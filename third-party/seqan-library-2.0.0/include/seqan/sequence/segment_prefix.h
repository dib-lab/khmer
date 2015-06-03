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
// Implementation of the Prefix Segment specialization.
// ==========================================================================

#ifndef SEQAN_HEADER_SEGMENT_PREFIX_H
#define SEQAN_HEADER_SEGMENT_PREFIX_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// PrefixSegment
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class PrefixSegment Prefix Segment
 * @extends Segment
 * @headerfile <seqan/sequence.h>
 * @brief A prefix of a sequence.
 *
 * @signature template <typename THost>
 *            class Segment<THost, PrefixSegment>;
 *
 * @tparam THost The underlying @link ContainerConcept sequence @endlink type.
 *
 * @see InfixSegment
 * @see SuffixSegment
 * @aka substring
 */

struct PrefixSegment {};

template <typename THost_>
class Segment<THost_, PrefixSegment>
{
public:
    typedef typename Host<Segment>::Type THost;

    typename Pointer_<THost>::Type data_host;
    typename Position<THost>::Type data_end_position;

//____________________________________________________________________________

public:

    Segment():
        data_host(),
        data_end_position(0)
    {}

    Segment(THost & _host):
        data_host(& _host),
        data_end_position(length(_host))
    {}

    Segment(typename Parameter_<THost>::Type _host, typename Position<THost>::Type _end_index):
        data_host(_toPointer(_host)),
        data_end_position(_end_index)
    {}
/*
    Segment(typename Parameter_<THost>::Type _host, typename Iterator<THost, Rooted>::Type _end):
        data_host(_toPointer(_host)),
        data_end_position(position(_end))
    {}
*/
    Segment(typename Parameter_<THost>::Type _host, typename Iterator<THost const, Standard>::Type _end):
        data_host(_toPointer(_host)),
        data_end_position(position(_end, _host))
    {}

/*
    Segment(Segment const & _other):
        data_host(_other.data_host),
        data_end_position(_other.data_end_position)
    {}
*/
    template <typename THost2, typename TSpec2>
    Segment(Segment<THost2, TSpec2> const & _other) :
        data_host(_toPointer(host(_other))),
        data_end_position(endPosition(_other))
    {}

    inline Segment &
    operator = (Segment const & source)
    {
        assign(*this, source);
        return *this;
    }
//____________________________________________________________________________

public:


    template <typename TPos>
    inline typename Reference<Segment>::Type
    operator [] (TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<Segment const>::Type
    operator [] (TPos pos) const
    {
        return value(*this, pos);
    }
//____________________________________________________________________________
};
//////////////////////////////////////////////////////////////////////////////

// template <typename THost>
// inline void
// clear(Segment<THost, PrefixSegment> & target)
// {
//     replace(host(target), beginPosition(target), endPosition(target), "");
//     setEndPosition(target, 0);
// }

template <typename THost_>
SEQAN_HOST_DEVICE inline typename Parameter_<THost_>::Type
host(Segment<THost_, PrefixSegment> & me)
{
SEQAN_CHECKPOINT
    return _toParameter<THost_>(me.data_host);
}

template <typename THost_>
SEQAN_HOST_DEVICE inline typename Parameter_<THost_>::Type
host(Segment<THost_, PrefixSegment> const & me)
{
SEQAN_CHECKPOINT
    return _toParameter<THost_>(me.data_host);
}

//____________________________________________________________________________

template <typename THost_>
inline void
setHost(Segment<THost_, PrefixSegment> & me, typename Parameter_<THost_>::Type _host)
{
SEQAN_CHECKPOINT
    me.data_host = _toPointer(_host);
}

template <typename THost_>
inline void
setHost(Segment<THost_ const, PrefixSegment> & me, typename Parameter_<THost_>::Type _host)
{
SEQAN_CHECKPOINT
    me.data_host = _toPointer(_host);
}

//____________________________________________________________________________

template <typename THost_>
SEQAN_HOST_DEVICE inline typename Iterator<Segment<THost_, PrefixSegment>, Standard>::Type
begin(Segment<THost_, PrefixSegment> & me,
    Standard)
{
    // TODO(holtgrew): An update at another point is required to make the call chain "return begin(host(me), Standard())" possible again. Same below.
    typename Parameter_<THost_>::Type tmpHost = host(me);
    return begin(tmpHost, Standard());
}
template <typename THost_>
SEQAN_HOST_DEVICE inline typename Iterator<Segment<THost_, PrefixSegment> const, Standard>::Type
begin(Segment<THost_, PrefixSegment> const & me,
    Standard)
{
    typename Parameter_<THost_>::Type tmpHost = host(me);
    return begin(tmpHost, Standard());
}

//____________________________________________________________________________

template <typename THost_>
SEQAN_HOST_DEVICE inline typename Position<Segment<THost_, PrefixSegment> const>::Type
beginPosition(Segment<THost_, PrefixSegment> const & /*me*/)
{
SEQAN_CHECKPOINT
    return 0;
}
template <typename THost_>
SEQAN_HOST_DEVICE inline typename Position<Segment<THost_, PrefixSegment> >::Type
beginPosition(Segment<THost_, PrefixSegment> & /*me*/)
{
SEQAN_CHECKPOINT
    return 0;
}

//____________________________________________________________________________

template <typename THost_, typename TIterator>
inline void
setBegin(Segment<THost_, PrefixSegment> &, TIterator)
{
}

//____________________________________________________________________________

template <typename THost_>
SEQAN_HOST_DEVICE inline typename Iterator<Segment<THost_, PrefixSegment>, Standard>::Type
end(Segment<THost_, PrefixSegment> & me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(host(me), Standard()) + me.data_end_position;
}
template <typename THost_>
SEQAN_HOST_DEVICE inline typename Iterator<Segment<THost_, PrefixSegment> const, Standard>::Type
end(Segment<THost_, PrefixSegment> const & me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(host(me), Standard()) + me.data_end_position;
}

//____________________________________________________________________________


template <typename THost_, typename TPosition>
inline void
setEndPosition(Segment<THost_, PrefixSegment> & me, TPosition new_end)
{
SEQAN_CHECKPOINT
    me.data_end_position = new_end;
}

template <typename THost_, typename TIterator>
inline void
setEnd(Segment<THost_, PrefixSegment> & me, TIterator new_end)
{
    SEQAN_CHECKPOINT;
    // me.data_end_position = new_end - begin(host(me));//, Standard());
    me.data_end_position = new_end - TIterator(begin(host(me)));
}

template <typename THost_>
inline void
setEnd(typename Iterator<Segment<THost_, PrefixSegment>, Rooted>::Type new_end)
{
SEQAN_CHECKPOINT
    container(new_end).data_end_position = hostIterator(new_end) - begin(host(container(new_end)));//, Standard());
}

//____________________________________________________________________________

template <typename THost_, typename TSize>
inline void
_setLength(
    Segment<THost_, PrefixSegment> & me,
    TSize new_length)
{
SEQAN_CHECKPOINT
    me.data_end_position = new_length;
}

//____________________________________________________________________________

template <typename THost_>
SEQAN_HOST_DEVICE inline typename Position<Segment<THost_, PrefixSegment> >::Type
endPosition(Segment<THost_, PrefixSegment> & me)
{
SEQAN_CHECKPOINT
    return me.data_end_position;
}
template <typename THost_>
SEQAN_HOST_DEVICE inline typename Position<Segment<THost_, PrefixSegment> const>::Type
endPosition(Segment<THost_, PrefixSegment> const & me)
{
SEQAN_CHECKPOINT
    return me.data_end_position;
}

//////////////////////////////////////////////////////////////////////////////

struct InfixSegment;
struct SuffixSegment;

template <typename THost>
struct Prefix
{
    typedef Segment<THost, PrefixSegment> Type;
};

template <typename THost>
struct Prefix< Segment<THost, InfixSegment> >
{
    typedef Segment<THost, InfixSegment> Type;
};
template <typename THost>
struct Prefix< Segment<THost, SuffixSegment> >
{
    typedef Segment<THost, InfixSegment> Type;
};
template <typename THost>
struct Prefix< Segment<THost, PrefixSegment> >
{
    typedef Segment<THost, PrefixSegment> Type;
};

template <typename THost, typename TSpec>
struct Prefix< Segment<THost, TSpec> const >:
    Prefix< Segment<THost, TSpec> > {};

template <typename THost>
struct Prefix<THost &>:
    Prefix<THost> {};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TPosition>
inline void
set(Segment<THost, PrefixSegment> & me,
    THost & host_,
    TPosition end_)
{
SEQAN_CHECKPOINT
    setHost(me, host_);
    setEndPosition(me, end_);
}
//____________________________________________________________________________

template <typename THost>
inline void
set(Segment<THost, PrefixSegment> & me,
    THost & host_)
{
SEQAN_CHECKPOINT
    setHost(me, host_);
    setEnd(me, end(host_));
}

//____________________________________________________________________________

template <typename THost, typename TSpec>
inline void
set(Segment<THost, PrefixSegment> & me,
    Segment<THost, TSpec> & source)
{
SEQAN_CHECKPOINT
    setHost(me, host(source));
    setEndPosition(me, endPosition(source));
}

template <typename THost, typename TSpec>
inline void
set(Segment<THost, PrefixSegment> & me,
    Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
    setHost(me, host(source));
    setEndPosition(me, endPosition(source));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atBegin(Segment<THost, PrefixSegment> const & segment)
{
SEQAN_CHECKPOINT
    return (endPosition(segment) == length(host(segment)));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atEnd(Segment<THost, PrefixSegment> const & segment)
{
SEQAN_CHECKPOINT
    return (endPosition(segment) == 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goBegin(Segment<THost, PrefixSegment> & segment,
        THost &)
{
SEQAN_CHECKPOINT
    goBegin(segment);
}

template <typename THost>
inline void
goBegin(Segment<THost, PrefixSegment> & segment)
{
    setEnd(segment);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goEnd(Segment<THost, PrefixSegment> & segment,
      THost &)
{
SEQAN_CHECKPOINT
    goEnd(segment);
}

template <typename THost>
inline void
goEnd(Segment<THost, PrefixSegment> & segment)
{
    setEnd(segment, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, PrefixSegment> &
operator ++(Segment<THost, PrefixSegment> & segment)
{
    setEnd(segment, endPosition(segment) - 1);
    return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, PrefixSegment> &
operator --(Segment<THost, PrefixSegment> & segment)
{
    setEnd(segment, endPosition(segment) + 1);
    return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename T, typename TPosEnd>
inline typename Prefix<T>::Type
prefix(T & t, TPosEnd pos_end)
{
    return typename Prefix<T>::Type(t, pos_end);
}

template <typename T, typename TPosEnd>
inline typename Prefix<T const>::Type
prefix(T const & t, TPosEnd pos_end)
{
    return typename Prefix<T const>::Type(t, pos_end);
}

template <typename T, typename TPosEnd>
inline typename Prefix<T *>::Type
prefix(T * t, TPosEnd pos_end)
{
    return typename Prefix<T *>::Type (t, pos_end);
}

//////////////////////////////////////////////////////////////////////////////
// A prefix of a prefix -> is a prefix
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, PrefixSegment> >::Type
prefix(Segment<T, PrefixSegment> & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, PrefixSegment> >::Type (
        host(t),
        beginPosition(t) + pos_end);
}
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, PrefixSegment> const>::Type
prefix(Segment<T, PrefixSegment> const & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, PrefixSegment> const>::Type (
        host(t),
        beginPosition(t) + pos_end);
}

//////////////////////////////////////////////////////////////////////////////
// A prefix of an infix -> is an infix
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, InfixSegment> >::Type
prefix(Segment<T, InfixSegment> & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, InfixSegment> >::Type (
        host(t),
        beginPosition(t),
        beginPosition(t) + pos_end);
}
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, InfixSegment> const>::Type
prefix(Segment<T, InfixSegment> const & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, InfixSegment> const>::Type (
        host(t),
        beginPosition(t),
        beginPosition(t) + pos_end);
}


//////////////////////////////////////////////////////////////////////////////
// A prefix of an suffix -> is an infix
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, SuffixSegment> >::Type
prefix(Segment<T, SuffixSegment> & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, SuffixSegment> >::Type (
        host(t),
        beginPosition(t),
        beginPosition(t) + pos_end);
}
template <typename T, typename TPosEnd>
inline typename Prefix<Segment<T, SuffixSegment> const>::Type
prefix(Segment<T, SuffixSegment> const & t, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, SuffixSegment> const>::Type (
        host(t),
        beginPosition(t),
        beginPosition(t) + pos_end);
}

//////////////////////////////////////////////////////////////////////////////

// prefix() with iterators

template <typename T, typename TIterSpec>
inline typename Prefix<T>::Type
prefix(T & t,
       Iter<Segment<T, PrefixSegment>, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Prefix<T>::Type(t, iterEnd);
}
template <typename T, typename TIterSpec>
inline typename Prefix<T const>::Type
prefix(T const & t,
       Iter<Segment<T, PrefixSegment>, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Prefix<T const>::Type(t, iterEnd);
}

template <typename T, typename TIterSpec>
inline typename Prefix<T *>::Type
prefix(T * t,
       Iter<Segment<T, PrefixSegment>, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Prefix<T *>::Type (t, iterEnd);
}

//////////////////////////////////////////////////////////////////////////////
// A prefix of a prefix -> is a prefix
template <typename T, typename TIterSpec>
inline typename Prefix<Segment<T, PrefixSegment> >::Type
prefix(Segment<T, PrefixSegment> & t,
       Iter<Segment<T, PrefixSegment>, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, PrefixSegment> >::Type (
        host(t),
        iterEnd);
}
template <typename T, typename TIterSpec>
inline typename Prefix<Segment<T, PrefixSegment> const>::Type
prefix(Segment<T, PrefixSegment> const & t,
       Iter<Segment<T, PrefixSegment> const, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, PrefixSegment> const>::Type (
        host(t),
        iterEnd);
}

//////////////////////////////////////////////////////////////////////////////
// A prefix of an infix -> is an infix
template <typename T, typename TIterSpec>
inline typename Prefix<Segment<T, InfixSegment> >::Type
prefix(Segment<T, InfixSegment> & t,
       Iter<Segment<T, InfixSegment>, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, InfixSegment> >::Type (
        host(t),
        begin(t),
        iterEnd);
}
template <typename T, typename TIterSpec>
inline typename Prefix<Segment<T, InfixSegment> const>::Type
prefix(Segment<T, InfixSegment> const & t,
       Iter<Segment<T, InfixSegment> const, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, InfixSegment> const>::Type (
        host(t),
        begin(t),
        iterEnd);
}


//////////////////////////////////////////////////////////////////////////////
// A prefix of an suffix -> is an infix
template <typename T, typename TIterSpec>
inline typename Prefix<Segment<T, SuffixSegment> >::Type
prefix(Segment<T, SuffixSegment> & t,
       Iter<Segment<T, SuffixSegment> const, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, SuffixSegment> >::Type (
        host(t),
        begin(t),
        iterEnd);
}
template <typename T, typename TIterSpec>
inline typename Prefix<Segment<T, SuffixSegment> const>::Type
prefix(Segment<T, SuffixSegment> const & t,
       Iter<Segment<T, SuffixSegment> const, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Prefix<Segment<T, SuffixSegment> const>::Type (
        host(t),
        begin(t),
        iterEnd);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
