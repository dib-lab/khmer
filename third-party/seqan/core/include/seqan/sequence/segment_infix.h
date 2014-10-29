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
// Implementation of the Infix Segment specialization.
// ==========================================================================

#ifndef SEQAN_HEADER_SEGMENT_INFIX_H
#define SEQAN_HEADER_SEGMENT_INFIX_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// InfixSegment
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class InfixSegment Infix Segment
 * @extends Segment
 * @headerfile <seqan/sequence.h>
 * @brief An infix of a sequence.
 *
 * @signature template <typename THost>
 *            class Segment<THost, InfixSegment>;
 *
 * @tparam THost The underlying @link SequenceConcept sequence@ type.
 */

/**
.Spec.InfixSegment:
..cat:Segments
..summary:An arbitrary segment.
..general:Class.Segment
..signature:Segment<THost, InfixSegment>
..param.THost:Type of the whole sequence.
...text:Instances of $Segment<THost, InfixSegment>$ are infixes of $THost$ objects.
...remarks:Use @Metafunction.Host@ to get the host type for a given class.
..remarks.note:Since the appropriate segment type depends on the host sequence type,
    it is recommended to use the metafunction @Metafunction.Infix@ instead of explicitely
    choose a specialization of @Class.Segment@.
..see:Metafunction.Infix
..example.file:demos/sequence/infix.cpp
..example.text:The output is as follows:
..example.output:
Infix: CGCG
..include:seqan/sequence.h
*/

template <typename THost_>
class Segment<THost_, InfixSegment>
{
public:
    typedef typename Host<Segment>::Type THost;

    typename Pointer_<THost>::Type data_host;
    typename Position<THost>::Type data_begin_position;
    typename Position<THost>::Type data_end_position;


//____________________________________________________________________________

public:
    // Check member variables with assertions.  This is called in the
    // constructors.
    void _checkMemberVariables() const {
        SEQAN_ASSERT_LEQ(data_begin_position, data_end_position);
    }

/**
.Memfunc.InfixSegment#Segment:
..class:Spec.InfixSegment
..summary:Constructor
..signature:Segment<THost, InfixSegment> ()
..signature:Segment<THost, InfixSegment> (infix)
..signature:Segment<THost, InfixSegment> (host [, begin, end])
..param.infix:Other infix object. (copy constructor)
..param.host:The whole sequence.
..param.begin:Position/iterator in $host$ of the first item in segment.
...type:Metafunction.Position.$Position<THost>::Type$
...type:Metafunction.Iterator.$Iterator<THost>::Type$
..param.end:Position/iterator behind the end of the segment.
...type:Metafunction.Position.$Position<THost>::Type$
...type:Metafunction.Iterator.$Iterator<THost>::Type$
..remarks:
...text:A Segment object cannot work without a host. If the object is default constructed,
the host must be set by @Function.setHost@ before the segment can be used.
...text:If a segment object is constructed by the copy constructor, the
members of the new constructed object are set to the same values as the members in the
source object; the host object is not modified.
Note that this is a special case, since all other copy operations result in changes
of the host object.
...text:$begin$ and $end$ must be valid positions/iterators in $host$.
...text:The predicate $begin <= end$ must be true.
*/
    Segment():
        data_host(),
        data_begin_position(0),
        data_end_position(0)
    {
SEQAN_CHECKPOINT
        _checkMemberVariables();
    }

    Segment(typename Parameter_<THost>::Type _host):
        data_host(_toPointer(_host)),
        data_begin_position(0),
        data_end_position(length(value(data_host)))
    {
SEQAN_CHECKPOINT
        _checkMemberVariables();
    }

    Segment(typename Parameter_<THost>::Type _host, typename Position<THost>::Type _begin_index, typename Position<THost>::Type _end_index):
        data_host(_toPointer(_host)),
        data_begin_position(_begin_index),
        data_end_position(_end_index)
    {
SEQAN_CHECKPOINT
        _checkMemberVariables();
    }
/*
    Segment(typename Parameter_<THost>::Type _host, typename Iterator<THost, Rooted>::Type _begin, typename Iterator<THost, Rooted>::Type _end):
        data_host(_toPointer(_host)),
        data_begin_position(position(_begin)),
        data_end_position(position(_end))
    {
SEQAN_CHECKPOINT
    }
*/
    Segment(typename Parameter_<THost>::Type _host, typename Iterator<THost, Standard>::Type _begin, typename Iterator<THost, Standard>::Type _end):
        data_host(_toPointer(_host)),
        data_begin_position(position(_begin, _host)),
        data_end_position(position(_end, _host))
    {
SEQAN_CHECKPOINT
        _checkMemberVariables();
    }
    template <typename THost2, typename TSpec2>
    Segment(Segment<THost2, TSpec2> const & _other):
        data_host(_toPointer(host(_other))),
        data_begin_position(beginPosition(_other)),
        data_end_position(endPosition(_other))
    {
SEQAN_CHECKPOINT
        _checkMemberVariables();
    }

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

};
//////////////////////////////////////////////////////////////////////////////


// template <typename THost>
// inline void
// clear(Segment<THost, InfixSegment> & target)
// {
//     replace(host(target), beginPosition(target), endPosition(target), "");
//     setEndPosition(target, beginPosition(target));
// }

///Function.host.param.object.type:Class.Segment

template <typename THost_>
inline typename Parameter_<THost_>::Type
host(Segment<THost_, InfixSegment> & me)
{
SEQAN_CHECKPOINT
    return _toParameter<THost_>(me.data_host);
}

template <typename THost_>
inline typename Parameter_<THost_>::Type
host(Segment<THost_, InfixSegment> const & me)
{
SEQAN_CHECKPOINT
    return _toParameter<THost_>(me.data_host);
}


//____________________________________________________________________________

///.Function.begin.param.object.type:Class.Segment
///.Function.begin.class:Class.Segment

template <typename THost_>
inline typename Iterator<Segment<THost_, InfixSegment>, Standard>::Type
begin(Segment<THost_, InfixSegment> & me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(host(me), Standard()) + me.data_begin_position;
}
template <typename THost_>
inline typename Iterator<Segment<THost_, InfixSegment> const, Standard>::Type
begin(Segment<THost_, InfixSegment> const & me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(host(me), Standard()) + me.data_begin_position;
}

//____________________________________________________________________________

///.Function.beginPosition.param.object.type:Class.Segment
///.Function.beginPosition.class:Class.Segment

template <typename THost_>
inline typename Position<Segment<THost_, InfixSegment> >::Type
beginPosition(Segment<THost_, InfixSegment> & me)
{
SEQAN_CHECKPOINT
    return me.data_begin_position;
}
template <typename THost_>
inline typename Position<Segment<THost_, InfixSegment> const>::Type
beginPosition(Segment<THost_, InfixSegment> const & me)
{
SEQAN_CHECKPOINT
    return me.data_begin_position;
}

//____________________________________________________________________________

/**
.Function.setBegin:
..class:Class.Segment
..summary:Sets begin of object in host.
..cat:Dependent Objects
..signature:setBegin(object, new_begin)
..param.object:An object.
...type:Spec.InfixSegment
...type:Spec.SuffixSegment
..param.new_begin:iterator to the new first item in $host(object)$ that belongs of $object$.
...type:Metafunction.Iterator
..see:Function.begin
..see:Function.beginPosition
..include:seqan/sequence.h
*/
template <typename THost_, typename TIterator>
inline void
setBegin(Segment<THost_, InfixSegment> & me, TIterator new_begin)
{
SEQAN_CHECKPOINT
    me.data_begin_position = new_begin - begin(host(me));//, Standard());
}


//____________________________________________________________________________

/**
.Function.setBeginPosition:
..class:Class.Segment
..summary:Sets begin position of object in host.
..cat:Dependent Objects
..signature:setBeginPosition(object, new_begin)
..param.object:An object.
...type:Spec.InfixSegment
...type:Spec.SuffixSegment
..param.new_begin:position of the new first item in $host(object)$ that belongs of $object$.
...type:Metafunction.Position
..see:Function.begin
..see:Function.beginPosition
..see:Function.setBegin
..include:seqan/sequence.h
*/

template <typename THost_, typename TPosition>
inline void
setBeginPosition(Segment<THost_, InfixSegment> & me, TPosition new_begin)
{
SEQAN_CHECKPOINT
    me.data_begin_position = new_begin;
}

//____________________________________________________________________________

///.Function.begin.param.object.type:Class.Segment
///.Function.begin.class:Class.Segment

template <typename THost_>
inline typename Iterator<Segment<THost_, InfixSegment>, Standard>::Type
end(Segment<THost_, InfixSegment> & me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(host(me), Standard()) + me.data_end_position;
}
template <typename THost_>
inline typename Iterator<Segment<THost_, InfixSegment> const, Standard>::Type
end(Segment<THost_, InfixSegment> const & me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(host(me), Standard()) + me.data_end_position;
}

//____________________________________________________________________________

///.Function.endPosition.param.object.type:Class.Segment
///.Function.endPosition.class:Class.Segment

template <typename THost_>
inline typename Position<Segment<THost_, InfixSegment> >::Type
endPosition(Segment<THost_, InfixSegment> & me)
{
SEQAN_CHECKPOINT
    return me.data_end_position;
}
template <typename THost_>
inline typename Position<Segment<THost_, InfixSegment> >::Type
endPosition(Segment<THost_, InfixSegment> const & me)
{
SEQAN_CHECKPOINT
    return me.data_end_position;
}

//____________________________________________________________________________

/**
.Function.setEnd:
..class:Class.Segment
..summary:Sets end of object in host.
..cat:Dependent Objects
..signature:setEnd(object, new_end)
..param.object:An object.
...type:Spec.InfixSegment
...type:Spec.PrefixSegment
..param.new_end:Iterator behind the last item in $host(object)$ belongs of $object$.
...type:Metafunction.Iterator
..see:Function.end
..see:Function.endPosition
..see:Function.setBegin
..include:seqan/sequence.h
*/

template <typename THost_, typename TIterator>
inline void
setEnd(Segment<THost_, InfixSegment> & me, TIterator new_end)
{
    SEQAN_CHECKPOINT;
    // me.data_end_position = new_end - begin(host(me)); //, Standard());
    me.data_end_position = new_end - TIterator(begin(host(me)));
}

//____________________________________________________________________________


/**
.Function.setEndPosition:
..class:Class.Segment
..summary:Sets begin position of object in host.
..cat:Dependent Objects
..signature:setEndPosition(object, new_end)
..param.object:An object.
...type:Spec.InfixSegment
...type:Spec.PrefixSegment
..param.new_end:position behind the last item in $host(object)$ that belongs of $object$.
...type:Metafunction.Position
..see:Function.end
..see:Function.endPosition
..see:Function.setBeginPosition
..see:Function.setEnd
..include:seqan/sequence.h
*/

template <typename THost_, typename TPosition>
inline void
setEndPosition(Segment<THost_, InfixSegment> & me, TPosition new_end)
{
SEQAN_CHECKPOINT
    me.data_end_position = new_end;
}

//____________________________________________________________________________

template <typename THost_>
inline void
_setLength(
    Segment<THost_, InfixSegment> & me,
    typename Size<THost_>::Type new_length)
{
SEQAN_CHECKPOINT
    me.data_end_position = me.data_begin_position + new_length;
}


//____________________________________________________________________________

/**
.Function.setHost:
..class:Class.Segment
..summary:Sets the host of an object.
..cat:Dependent Objects
..signature:setHost(object, host)
..param.object:The object that will get a new host.
...type:Class.Segment
..param.host:The new host.
..remarks:After this operation, $object$ depends on $host$.
...text:Note that setting the host can invalidate $object$.
For example, if one changes the host of a @Class.Segment@ object, it is possible
that begin- and end-position of the segment does not fit into the new host sequence.
..see:Function.host
..include:seqan/sequence.h
*/
template <typename THost_>
inline void
setHost(Segment<THost_, InfixSegment> & me, typename Parameter_<THost_>::Type _host)
{
SEQAN_CHECKPOINT
    me.data_host = _toPointer(_host);
}

template <typename THost_>
inline void
setHost(Segment<THost_ const, InfixSegment> & me, typename Parameter_<THost_>::Type _host)
{
SEQAN_CHECKPOINT
    me.data_host = _toPointer(_host);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.Infix:
..cat:Segments
..class:Class.String
..summary:Infix sequence type.
..signature:Infix<T>::Type
..remarks:Note that an infix of a @Class.Segment@ object is an @Spec.InfixSegment@ object having the same host type.
..param.T:A sequence type.
...type:Class.String
..returns.param.Type:The infix type.
..see:Spec.InfixSegment
..include:seqan/sequence.h
*/

template <typename THost>
struct Infix
{
    typedef Segment<THost, InfixSegment> Type;
};

template <typename THost, typename TSpec>
struct Infix< Segment<THost, TSpec> >
{
    typedef Segment<THost, InfixSegment> Type;
};

template <typename THost, typename TSpec>
struct Infix< Segment<THost, TSpec> const >:
    Infix< Segment<THost, TSpec> > {};

template <typename THost>
struct Infix<THost &>:
    Infix<THost> {};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TPosition1, typename TPosition2>
inline void
set(Segment<THost, InfixSegment> & me,
    THost & host_,
    TPosition1 begin_,
    TPosition2 end_)
{
SEQAN_CHECKPOINT
    setHost(me, host_);
    setBeginPosition(me, begin_);
    setEndPosition(me, end_);
}
//____________________________________________________________________________

template <typename THost>
inline void
set(Segment<THost, InfixSegment> & me,
    THost & host_)
{
SEQAN_CHECKPOINT
    setHost(me, host_);
    setBegin(me, begin(host_, Standard()));
    setEnd(me, end(host_, Standard()));
}
template <typename THost>
inline void
set(Segment<THost, InfixSegment> & me,
    THost const & host_)
{
SEQAN_CHECKPOINT
    setHost(me, host_);
    setBegin(me, begin(host_, Standard()));
    setEnd(me, end(host_, Standard()));
}

//____________________________________________________________________________

template <typename THost, typename TSpec>
inline void
set(Segment<THost, InfixSegment> & me,
    Segment<THost, TSpec> & source)
{
SEQAN_CHECKPOINT
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
    setEndPosition(me, endPosition(source));
}
template <typename THost, typename TSpec>
inline void
set(Segment<THost const, InfixSegment> & me,
    Segment<THost, TSpec> & source)
{
SEQAN_CHECKPOINT
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
    setEndPosition(me, endPosition(source));
}
template <typename THost, typename TSpec>
inline void
set(Segment<THost, InfixSegment> & me,
    Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
    setEndPosition(me, endPosition(source));
}
template <typename THost, typename TSpec>
inline void
set(Segment<THost const, InfixSegment> & me,
    Segment<THost, TSpec> const & source)
{
SEQAN_CHECKPOINT
    setHost(me, host(source));
    setBeginPosition(me, beginPosition(source));
    setEndPosition(me, endPosition(source));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atBegin(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
    return (beginPosition(segment) == endPosition(segment));
}
template <typename THost>
inline bool
atBegin(Segment<THost, InfixSegment> const & segment)
{
SEQAN_CHECKPOINT
    return (beginPosition(segment) == endPosition(segment));
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline bool
atEnd(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
    return (endPosition(segment) - beginPosition(segment)) > length(host(segment));
}
template <typename THost>
inline bool
atEnd(Segment<THost, InfixSegment> const & segment)
{
SEQAN_CHECKPOINT
    return (endPosition(segment) - beginPosition(segment)) > length(host(segment));
}


//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline void
goBegin(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
    setBeginPosition(segment, 0);
    setEndPosition(segment, 1);
}
template <typename THost, typename THost2>
inline void
goBegin(Segment<THost, InfixSegment> & segment,
        THost2 &)
{
    goBegin(segment);
}
template <typename THost, typename THost2>
inline void
goBegin(Segment<THost, InfixSegment> & segment,
        THost2 const &)
{
    goBegin(segment);
}

//////////////////////////////////////////////////////////////////////////////


template <typename THost>
inline void
goEnd(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
    setBeginPosition(segment, 0);
    setEndPosition(segment, length(host(segment)));
}
template <typename THost, typename THost2>
inline void
goEnd(Segment<THost, InfixSegment> & segment,
      THost2 &)
{
    goEnd(segment);
}
template <typename THost, typename THost2>
inline void
goEnd(Segment<THost, InfixSegment> & segment,
      THost2 const &)
{
    goEnd(segment);
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, InfixSegment> &
operator ++(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
    if (endPosition(segment) == length(host(segment)))
    {
        setEndPosition(segment, endPosition(segment) - beginPosition(segment) + 1);
        setBeginPosition(segment, 0);
    }
    else
    {
        setBeginPosition(segment, beginPosition(segment) + 1);
        setEndPosition(segment, endPosition(segment) + 1);
    }
    return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost>
inline Segment<THost, InfixSegment> &
operator --(Segment<THost, InfixSegment> & segment)
{
SEQAN_CHECKPOINT
    if (!beginPosition(segment))
    {
        typename Size<THost>::Type host_length = length(host(segment));

        setBeginPosition(segment, host_length - endPosition(segment) + beginPosition(segment) + 1);
        setEndPosition(segment, host_length);
    }
    else
    {
        setBeginPosition(segment, beginPosition(segment) - 1);
        setEndPosition(segment, endPosition(segment) - 1);
    }
    return segment;
}

//////////////////////////////////////////////////////////////////////////////

template <typename THost, typename TSpec, typename TPos>
inline typename Reference< Segment<THost, TSpec> >::Type
value(Segment<THost, TSpec> & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LT_MSG(pos, static_cast<TPos>(length(me)), "Trying to acces an element behind the last one!");
    return *(begin(me, Standard()) + pos);
}

template <typename THost, typename TSpec, typename TPos>
inline typename Reference< Segment<THost, TSpec> const >::Type
value(Segment<THost, TSpec> const & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LT_MSG(pos, static_cast<TPos>(length(me)), "Trying to acces an element behind the last one!");
    return *(begin(me, Standard()) + pos);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.infix:
..cat:Containers
..class:Class.String
..class:Adaption.char array
..summary:Creates infix object.
..signature:infix(host, begin, end)
..param.host:The complete sequence.
...type:Class.String
...type:Adaption.char array
..param.begin:Position or iterator of the first element of the segment.
...type:Metafunction.Position
...type:Metafunction.Iterator
..param.end:Position or iterator behind the last element of the segment.
...remarks:$end$ must have the same type as $begin$.
..returns:The infix of $host$ between $begin$ and $end-1$.
...remarks:The type of the infix is given by @Metafunction.Infix@.
..remarks:Notational sugar.
..see:Spec.InfixSegment
..see.Metafunction.Infix
..include:seqan/sequence.h
*/

template <typename T, typename TPosBegin, typename TPosEnd>
inline typename Infix<T>::Type
infix(T & t, TPosBegin pos_begin, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Infix<T>::Type(t, pos_begin, pos_end);
}

template <typename T, typename TPosBegin, typename TPosEnd>
inline typename Infix<T *>::Type
infix(T * t, TPosBegin pos_begin, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Infix<T *>::Type (t, pos_begin, pos_end);
}

template <typename T, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Segment<T, TSpec> >::Type
infix(Segment<T, TSpec> & t, TPosBegin pos_begin, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Infix<Segment<T, TSpec> >::Type (
        host(t),
        beginPosition(t) + pos_begin,
        beginPosition(t) + pos_end);
}

template <typename T, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<Segment<T, TSpec> const>::Type
infix(Segment<T, TSpec> const & t, TPosBegin pos_begin, TPosEnd pos_end)
{
SEQAN_CHECKPOINT
    return typename Infix<Segment<T, TSpec> const>::Type (
        host(t),
        beginPosition(t) + pos_begin,
        beginPosition(t) + pos_end);
}

// infix() with iterators

template <typename T, typename TSpec, typename TIterSpec>
inline typename Infix<Segment<T, TSpec> >::Type
infix(Segment<T, TSpec> & t,
      Iter<Segment<T, TSpec>, TIterSpec> const & iterBegin,
      Iter<Segment<T, TSpec>, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Infix<Segment<T, TSpec> >::Type (
        host(t),
        iterBegin,
        iterEnd);
}

template <typename T, typename TSpec, typename TIterSpec>
inline typename Infix<Segment<T, TSpec> const>::Type
infix(Segment<T, TSpec> const & t,
      Iter<Segment<T, TSpec> const, TIterSpec> const & iterBegin,
      Iter<Segment<T, TSpec> const, TIterSpec> const & iterEnd)
{
SEQAN_CHECKPOINT
    return typename Infix<Segment<T, TSpec> >::Type (
        host(t),
        iterBegin,
        iterEnd);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.infixWithLength:
..cat:Containers
..class:Class.String
..class:Adaption.char array
..summary:Creates infix object.
..signature:infixWithLength(host, begin, length)
..param.host:The complete sequence.
...type:Class.String
...type:Adaption.char array
..param.begin:Position or iterator of the first element of the segment.
...type:Metafunction.Position
...type:Metafunction.Iterator
..param.length:Length of the returned infix.
..returns:The infix of $host$ between $begin$ and $begin+length-1$.
...remarks:The type of the infix is given by @Metafunction.Infix@.
..remarks:Notational sugar.
..see:Spec.InfixSegment
..include:seqan/sequence.h
*/

template <typename T, typename TPosBegin, typename TSize>
inline typename Infix<T>::Type
infixWithLength(T & t, TPosBegin pos_begin, TSize length)
{
SEQAN_CHECKPOINT
    return typename Infix<T>::Type(t, pos_begin, pos_begin + length);
}

template <typename T, typename TPosBegin, typename TSize>
inline typename Infix<T *>::Type
infixWithLength(T * t, TPosBegin pos_begin, TSize length)
{
SEQAN_CHECKPOINT
    return typename Infix<T *>::Type (t, pos_begin, pos_begin + length);
}

template <typename T, typename TSpec, typename TPosBegin, typename TSize>
inline typename Infix<Segment<T, TSpec> >::Type
infixWithLength(Segment<T, TSpec> & t, TPosBegin pos_begin, TSize length)
{
SEQAN_CHECKPOINT
    return typename Infix<Segment<T, TSpec> >::Type (
        host(t),
        beginPosition(t) + pos_begin,
        beginPosition(t) + pos_begin + length);
}

template <typename T, typename TSpec, typename TPosBegin, typename TSize>
inline typename Infix<Segment<T, TSpec> const>::Type
infixWithLength(Segment<T, TSpec> const & t, TPosBegin pos_begin, TSize length)
{
SEQAN_CHECKPOINT
    return typename Infix<Segment<T, TSpec> const>::Type (
        host(t),
        beginPosition(t) + pos_begin,
        beginPosition(t) + pos_begin + length);
}

//////////////////////////////////////////////////////////////////////////////
//setBegin


template <typename TIterator>
inline void
setBegin(TIterator new_begin)
{
SEQAN_CHECKPOINT
    setBegin(container(new_begin), hostIterator(new_begin));
}


//////////////////////////////////////////////////////////////////////////////
//setEnd

template <typename TIterator>
inline void
setEnd(TIterator new_end)
{
SEQAN_CHECKPOINT
    setEnd(container(new_end), new_end);
}

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
