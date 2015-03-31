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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_MODIFIER_REVERSE_H
#define SEQAN_HEADER_MODIFIER_REVERSE_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace seqan
{

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Classes, Enums, Typedefs
// ==========================================================================

// --------------------------------------------------------------------------
// Class ModReverse Iterator
// --------------------------------------------------------------------------

/**
.Spec.ModReverse:
..summary:Mirrors the characters from begin to end.
..cat:Modifier
..general:Class.ModifiedIterator
..general:Class.ModifiedString
..signature:ModifiedIterator<THost, ModReverse>
..signature:ModifiedString<THost, ModReverse>
..param.THost:Original string/iterator.
...type:Concept.RandomAccessIteratorConcept
..include:seqan/modifier.h
*/

struct ModReverse_;
typedef Tag<ModReverse_> ModReverse;

template <typename THost>
class ModifiedIterator<THost, ModReverse>
{
public:
    typedef typename Cargo<ModifiedIterator>::Type TCargo_;

    Holder<THost, Simple> _host;
    TCargo_ _cargo;

    ModifiedIterator() : _host(), _cargo()
    {}

    ModifiedIterator(ModifiedIterator &_origin) :
			_host(_origin._host), _cargo(_origin._cargo)
    {}

    ModifiedIterator(ModifiedIterator const & _origin) :
			_host(_origin._host), _cargo(_origin._cargo)
    {}

    template <typename T>
    explicit
    ModifiedIterator(T & host) : _host(host)
    {}

    template <typename T>
    explicit
    ModifiedIterator(T const & host) : _host(host)
    {}
};


template <typename THost>
class ModifiedString<THost, ModReverse>
{
public:
    typedef typename Pointer_<THost>::Type       THostPointer_;
	typedef typename Cargo<ModifiedString>::Type TCargo_;

    typedef typename InnermostHost_<ModifiedString>::Type TInnermostHost_;

    mutable THostPointer_ _host;
    TCargo_ _cargo;

    // Default constructor.
    ModifiedString() : _host(), _cargo()
    {}

    // Construct with the actual host.
    explicit
    ModifiedString(THost & host) : _host(_toPointer(host)), _cargo()
    {}

    // Constructor for creating a ModifiedString with const host with a non-const host.
    template <typename THost_>
    explicit ModifiedString(THost_ const & host,
                            SEQAN_CTOR_ENABLE_IF(IsSameType<THost, THost_>)) :
            _host(_toPointer(host)), _cargo()
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Non-const variant.
    template <typename THost_>
    explicit
    ModifiedString(THost_ & host,
                   SEQAN_CTOR_ENABLE_IF(And<Not<IsSameType<TInnermostHost_, THost> >,
                                            IsSameType<TInnermostHost_, THost_> >)) :
            _host(host), _cargo()
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Const variant.
    template <typename THost_>
    explicit
    ModifiedString(THost_ const & host,
                   SEQAN_CTOR_ENABLE_IF(And<Not<IsSameType<TInnermostHost_, THost> >,
                                            IsSameType<TInnermostHost_, THost_> >)) :
            _host(host), _cargo()
    {
        ignoreUnusedVariableWarning(dummy);
    }

    template <typename TPos>
    inline typename Reference<ModifiedString>::Type 
    operator[](TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<ModifiedString const>::Type 
    operator[](TPos pos) const
    {
        return value(*this, pos);
    }
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction Cargo                           [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
struct Cargo<ModifiedIterator<THost, ModReverse> >
{
    typedef Cargo Type;		// to reduce namespace pollution
    bool _atEnd;

    Cargo() : _atEnd(false)
    {}
};

// --------------------------------------------------------------------------
// Metafunction Iterator                          [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost>
struct Iterator<ModifiedString<THost, ModReverse>, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost const, Rooted>::Type, ModReverse> Type;
};

template <typename THost>
struct Iterator<ModifiedString<THost, ModReverse> const, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost const, Rooted>::Type, ModReverse> Type;
};

// --------------------------------------------------------------------------
// Metafunction DefaultIteratorSpec               [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost>
struct DefaultIteratorSpec< ModifiedString<THost, ModReverse> >
{
    typedef Rooted Type;
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function goNext()                            [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline void
goNext(ModifiedIterator<THost, ModReverse> & me)
{
    if (atBegin(host(me)))
        cargo(me)._atEnd = true;
    else
        goPrevious(host(me));
}

// --------------------------------------------------------------------------
// Function goPrevious()                        [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline void
goPrevious(ModifiedIterator<THost, ModReverse> & me)
{
    if (cargo(me)._atEnd)
        cargo(me)._atEnd = false;
    else
        goNext(host(me));
}

// --------------------------------------------------------------------------
// Function goEnd()                             [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline void
goEnd(ModifiedIterator<THost, ModReverse> & me)
{
    goBegin(host(me));
    cargo(me)._atEnd = true;
}

// --------------------------------------------------------------------------
// Function goBegin()                           [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline void
goBegin(ModifiedIterator<THost, ModReverse> & me)
{
    goEnd(host(me));
    if (atBegin(host(me)))
    {
        cargo(me)._atEnd = true;
    }
    else
    {
        cargo(me)._atEnd = false;
        goPrevious(host(me));
    }
}

// --------------------------------------------------------------------------
// Function operator+=()                        [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TDelta>
inline ModifiedIterator<THost, ModReverse> &
operator+=(ModifiedIterator<THost, ModReverse> & me, TDelta delta_) 
{
    typedef ModifiedIterator<THost, ModReverse> TIterator;
    typedef typename Position<TIterator>::Type TPosition;
    TPosition delta = delta_;
    
    if (delta == 0)
    {
        return me;
    }
    if (delta > 0)
    {
        if (position(host(me)) < delta)
        {
            cargo(me)._atEnd = true;
            --delta;
        }
        host(me) -= delta;
    }
    else
    {
        if (cargo(me)._atEnd)
        {
            cargo(me)._atEnd = false;
            ++delta;
        }
        host(me) -= delta;
    } 
    return me;
}

// --------------------------------------------------------------------------
// Function operator-=()                        [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TDelta>
inline ModifiedIterator<THost, ModReverse> &
operator-=(ModifiedIterator<THost, ModReverse> & me, TDelta delta)
{
    if (delta > 0)
    {
        if (cargo(me)._atEnd)
        {
            cargo(me)._atEnd = false;
            --delta;
        }
        host(me) += delta;
    }
    else
    {
        if (position(host(me)) < -delta)
        {
            cargo(me)._atEnd = true;
            ++delta;
        }
        host(me) -= -delta;
    }
    return me;
}

// --------------------------------------------------------------------------
// Function operator-()                         [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline typename Difference< ModifiedIterator<THost, ModReverse> >::Type
operator-(ModifiedIterator<THost, ModReverse> const & a,
          ModifiedIterator<THost, ModReverse> const & b)
{
    typename Difference< ModifiedIterator<THost, ModReverse> >::Type diff = host(b) - host(a);
    if (cargo(a)._atEnd)
        ++diff;
    if (cargo(b)._atEnd)
        --diff;
    return diff;
}

// --------------------------------------------------------------------------
// Function position()                          [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline typename Position<ModifiedIterator<THost, ModReverse> const>::Type 
position(ModifiedIterator<THost, ModReverse> const & me)
{
    if (cargo(me)._atEnd)
        return length(container(host(me)));
    else
        return length(container(host(me))) - 1 - position(host(me));
}

template <typename THost, typename TContainer>
inline typename Position<ModifiedIterator<THost, ModReverse> const>::Type 
position(ModifiedIterator<THost, ModReverse> const & me, TContainer const &cont)
{
    if (cargo(me)._atEnd)
        return length(cont);
    else
        return length(cont) - 1 - position(host(me), cont);
}

// --------------------------------------------------------------------------
// Function setPosition()                       [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TPosition>
inline void
setPosition(ModifiedIterator<THost, ModReverse> const & me, TPosition pos)
{
    setPosition(host(me), length(container(host(me))) - 1 - pos);
}

// --------------------------------------------------------------------------
// Function operator==()                        [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline bool
operator==(ModifiedIterator<THost, ModReverse> const & a,
           ModifiedIterator<THost, ModReverse> const & b)
{
    return cargo(a)._atEnd == cargo(b)._atEnd && host(a) == host(b);
}

// --------------------------------------------------------------------------
// Function operator<()                         [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline bool
operator<(ModifiedIterator<THost, ModReverse> const & a,
          ModifiedIterator<THost, ModReverse> const & b)
 {
    return (!cargo(a)._atEnd && cargo(b)._atEnd) ||
            (!cargo(a)._atEnd && !cargo(b)._atEnd && host(a) > host(b));
}

// --------------------------------------------------------------------------
// Function atEnd()                             [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TContainer>
inline bool
atBegin(ModifiedIterator<THost, ModReverse> const & me,
        TContainer const & container)
{
    return position(me, container) == 0;
}

template <typename THost>
inline bool
atBegin(ModifiedIterator<THost, ModReverse> const & me)
{
    return position(me) == 0;
}

// --------------------------------------------------------------------------
// Function atEnd()                             [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TContainer>
inline bool
atEnd(ModifiedIterator<THost, ModReverse> const & me,
      TContainer const & /*container*/)
{
            return cargo(me)._atEnd;
}

template <typename THost>
inline bool
atEnd(ModifiedIterator<THost, ModReverse> const & me)
{
            return cargo(me)._atEnd;
}

// --------------------------------------------------------------------------
// Function value()                               [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost, typename TPos>
inline typename Reference<ModifiedString<THost, ModReverse> >::Type 
value(ModifiedString<THost, ModReverse> & me, TPos pos)
{
    return value(host(me), (length(host(me)) - 1) - pos);
}

template <typename THost, typename TPos>
inline typename Reference<ModifiedString<THost, ModReverse> const>::Type 
value(ModifiedString<THost, ModReverse> const & me, TPos pos)
{
    return value(host(me), (length(host(me)) - 1) - pos);
}

// --------------------------------------------------------------------------
// Function begin()                               [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template < typename THost, typename TTag >
inline typename Iterator< ModifiedString<THost, ModReverse> const >::Type 
begin(ModifiedString<THost, ModReverse> const & me)
{
    typename Iterator< ModifiedString<THost, ModReverse> const >::Type temp_(end(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template < typename THost >
inline typename Iterator< ModifiedString<THost, ModReverse> >::Type 
begin(ModifiedString<THost, ModReverse> & me)
{
    typename Iterator< ModifiedString<THost, ModReverse> >::Type temp_(end(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template < typename THost, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type 
begin(ModifiedString<THost, ModReverse> const & me, Tag<TTagSpec> const)
{
    typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type temp_(end(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template < typename THost, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type 
begin(ModifiedString<THost, ModReverse> & me, Tag<TTagSpec> const)
{
    typedef typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type TIterator;
    TIterator temp_(end(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

// --------------------------------------------------------------------------
// Function end()                                 [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost >
inline typename Iterator<ModifiedString<THost, ModReverse> const >::Type 
end(ModifiedString<THost, ModReverse> const & me)
{
    typename Iterator<ModifiedString<THost, ModReverse> const >::Type temp_(begin(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template <typename THost >
inline typename Iterator<ModifiedString<THost, ModReverse> >::Type 
end(ModifiedString<THost, ModReverse> & me)
{
    typename Iterator<ModifiedString<THost, ModReverse> >::Type temp_(begin(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template <typename THost, typename TTagSpec >
inline typename Iterator<ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const>::Type 
end(ModifiedString<THost, ModReverse> const & me, Tag<TTagSpec> const)
{
    typename Iterator<ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type temp_(begin(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template <typename THost, typename TTagSpec >
inline typename Iterator<ModifiedString<THost, ModReverse>, Tag<TTagSpec> const>::Type 
end(ModifiedString<THost, ModReverse> & me, Tag<TTagSpec> const)
{
    typename Iterator<ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type temp_(begin(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

// --------------------------------------------------------------------------
// Function reverse()
// --------------------------------------------------------------------------

/**
.Function.reverse
..summary:Reverse an object/container in-place.
..cat:Modifier
..signature:reverse(object)
..param.object:The object/container whose elements to reverse.
...type:Concept.ContainerConcept
...type:Adaption.std::list
..include:seqan/modifier.h
*/

template < typename TSequence >
inline void
reverse(TSequence & sequence) 
{
    typedef typename Value<TSequence>::Type					TValue;

#if defined (_OPENMP) && defined (SEQAN_PARALLEL)
    // OpenMP does not support for loop with iterators. Therefore use index variables.
    typedef typename Position<TSequence>::Type				TPos;
    typedef typename MakeSigned_<TPos>::Type				TSignedPos;

    TSignedPos pMid = length(sequence) / 2;

#pragma omp parallel for if(length(sequence) > 1000000)
    for(TSignedPos p1 = 0; p1 < pMid; ++p1) {
        TPos p2 = length(sequence) - 1 - p1;
        TValue tmp = sequence[p1];
        sequence[p1] = sequence[p2];
        sequence[p2] = tmp;
    }
#else
    typedef typename Iterator<TSequence, Standard>::Type	TIter;
    TIter it1 = begin(sequence, Standard());
    TIter it2 = it1 + (length(sequence) - 1);
    TIter itMid = it1 + length(sequence) / 2;

    for(; it1 != itMid; ++it1, --it2) {
        TValue tmp = *it1;
        *it1 = *it2;
        *it2 = tmp;
    }
#endif
}

template < typename TSequence >
inline void
reverse(TSequence const & sequence) 
{
    reverse(const_cast<TSequence &>(sequence));
}

template < typename TSequence, typename TSpec >
inline void
reverse(StringSet<TSequence, TSpec> & stringSet) 
{
    unsigned seqCount = length(stringSet);
    for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
        reverse(stringSet[seqNo]);
}

template < typename TSequence, typename TSpec >
inline void
reverse(StringSet<TSequence, TSpec> const & stringSet) 
{
    unsigned seqCount = length(stringSet);
    for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
        reverse(stringSet[seqNo]);
}

template <typename TValue>
inline void
reverse(std::list<TValue> & list)
{
    list.reverse();
}

// --------------------------------------------------------------------------
// Function reverseString()
// --------------------------------------------------------------------------

template <typename THost>
inline ModifiedString<THost, ModReverse>
reverseString(THost & host)
{
	return ModifiedString<THost, ModReverse>(host);
}

}

#endif
