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
// Author: Knut Reinert <knut.reinert@fu-berlin.de>
// ==========================================================================
// Adaptions for STL vectors to SeqAn strings.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS: No forwards are generated for this file.

#ifndef SEQAN_SEQUENCE_ADAPT_STD_VECTOR_H_
#define SEQAN_SEQUENCE_ADAPT_STD_VECTOR_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Concepts
// ===========================================================================

// ----------------------------------------------------------------------------
// Concept StringConcept
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::vector<TChar, TAlloc>), (StringConcept));          // resizable container

template <typename TChar, typename TAlloc>
SEQAN_CONCEPT_IMPL((std::vector<TChar, TAlloc> const), (ContainerConcept)); // read-only container

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TChar, typename TAlloc>
struct IsContiguous< std::vector<TChar, TAlloc> >
{
    enum { VALUE = true };
};

template <typename  TChar, typename TAlloc>
struct IsContiguous< std::vector<TChar, TAlloc> const>
        : IsContiguous< std::vector<TChar, TAlloc> > {};

template <typename TChar, typename TAlloc>
struct Value< std::vector<TChar, TAlloc> >
{
    typedef typename std::vector<TChar, TAlloc>::value_type Type;
};

template <typename TChar, typename TAlloc>
struct Value< std::vector<TChar, TAlloc> const>
        : Value< std::vector<TChar, TAlloc> > {};

// TODO(holtgrew): GetValue is a reference?! I thought the reverse was true in respect to Value<>.
template <typename TChar, typename TAlloc>
struct GetValue< std::vector<TChar, TAlloc> >
{
    typedef typename std::vector<TChar, TAlloc>::reference Type;
};

template <typename TChar, typename TAlloc>
struct GetValue< std::vector<TChar,  TAlloc> const>
{
    typedef typename std::vector<TChar, TAlloc>::const_reference Type;
};

template <typename TChar, typename TAlloc>
struct Reference< std::vector<TChar, TAlloc> >
{
    typedef typename std::vector<TChar, TAlloc>::reference Type;
};

template <typename TChar,  typename TAlloc>
struct Reference< std::vector<TChar, TAlloc> const>
{
    typedef typename std::vector<TChar,  TAlloc>::const_reference Type;
};

template <typename TChar, typename TAlloc>
struct Iterator< std::vector<TChar, TAlloc>, Rooted>
{
    typedef std::vector<TChar, TAlloc> TVector_;
    typedef Iter<TVector_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TVector_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar, typename TAlloc>
struct Iterator< std::vector<TChar, TAlloc> const, Rooted>
{
    typedef std::vector<TChar, TAlloc> const TVector_;
    typedef Iter<TVector_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TVector_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar,  typename TAlloc>
struct Iterator< std::vector<TChar, TAlloc>, Standard >
{
    typedef Iter< std::vector<TChar,  TAlloc>, StdIteratorAdaptor > Type;
};

template <typename TChar,  typename TAlloc>
struct Iterator< std::vector<TChar,  TAlloc> const, Standard>
{
    typedef Iter< std::vector<TChar, TAlloc> const, StdIteratorAdaptor > Type;
};

template <typename TChar,  typename TAlloc>
struct Position< std::vector<TChar, TAlloc> >
{
    typedef typename std::vector<TChar,  TAlloc>::size_type Type;
};

template <typename TChar,  typename TAlloc>
struct Position< std::vector<TChar,  TAlloc> const>
        : Position< std::vector<TChar,  TAlloc> > {};

template <typename TChar,  typename TAlloc>
struct Size< std::vector<TChar, TAlloc> >
{
    typedef typename std::vector<TChar, TAlloc>::size_type Type;
};

template <typename TChar, typename TAlloc>
struct Size< std::vector<TChar, TAlloc> const>
        : Size< std::vector<TChar, TAlloc> > {};

template <typename TChar, typename TAlloc>
struct DefaultOverflowImplicit< std::vector<TChar, TAlloc> >
{
    typedef Generous Type;
};

template <typename TChar, typename TAlloc>
struct StdContainerIterator< std::vector<TChar, TAlloc> >
{
    typedef std::vector<TChar, TAlloc> TContainer_;
    typedef typename TContainer_::iterator Type;
};

template <typename TChar, typename TAlloc>
struct StdContainerIterator< std::vector<TChar, TAlloc> const>
{
    typedef std::vector<TChar, TAlloc> TContainer_;
    typedef typename TContainer_::const_iterator Type;
};

template <typename TChar, typename TAlloc>
struct IsSequence<std::vector<TChar, TAlloc> > : True {};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TChar,  typename TAlloc>
inline void const *
getObjectId(std::vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    if (me.empty())
        return NULL;
    else
        return (& *(me.end() - 1)) + 1;
}

template <typename TChar,  typename TAlloc>
inline typename Iterator< std::vector<TChar,  TAlloc>, Standard>::Type
begin(std::vector<TChar,  TAlloc> & me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::vector<TChar,  TAlloc>, Standard>::Type(me.begin());
}
template <typename TChar,  typename TAlloc>
inline typename Iterator< std::vector<TChar,  TAlloc> const, Standard>::Type
begin(std::vector<TChar, TAlloc> const & me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::vector<TChar,  TAlloc> const, Standard>::Type(me.begin());
}

template <typename TChar, typename TAlloc>
inline typename Iterator< std::vector<TChar, TAlloc>, Standard>::Type
end(std::vector<TChar,  TAlloc> & me,
    Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::vector<TChar, TAlloc>, Standard>::Type(me.end());
}
template <typename TChar,  typename TAlloc>
inline typename Iterator< std::vector<TChar,  TAlloc> const, Standard>::Type
end(std::vector<TChar,  TAlloc> const & me,
    Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::vector<TChar,  TAlloc> const, Standard>::Type(me.end());
}

template <typename TChar,  typename TAlloc, typename TPos>
inline typename GetValue< std::vector<TChar, TAlloc> >::Type
value(std::vector<TChar,  TAlloc> & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    return me[pos];
}
template <typename TChar,  typename TAlloc, typename TPos>
inline typename GetValue< std::vector<TChar,  TAlloc> const>::Type
value(std::vector<TChar, TAlloc> const & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    return me[pos];
}

template <typename TChar, typename TAlloc>
inline typename Size< std::vector<TChar, TAlloc> >::Type
length(std::vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    return me.size();
}

template <typename TChar, typename TAlloc>
inline typename Size< std::vector<TChar, TAlloc> >::Type
capacity(std::vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    return me.capacity();
}

template <typename TChar, typename TAlloc>
inline bool
empty(std::vector<TChar, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    return me.empty();
}

template <typename TChar,  typename TAlloc>
inline void
clear(std::vector<TChar, TAlloc> & me)
{
    SEQAN_CHECKPOINT;
    me.clear();
}

template <typename TChar>
inline typename Reference<std::vector<TChar> >::Type
front(std::vector<TChar> & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

template <typename TChar>
inline typename Reference<std::vector<TChar> const>::Type
front(std::vector<TChar> const & list)
{
    SEQAN_CHECKPOINT;
    return list.front();
}

template <typename TChar>
inline typename Reference<std::vector<TChar> >::Type
back(std::vector<TChar> & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

template <typename TChar>
inline typename Reference<std::vector<TChar> const>::Type
back(std::vector<TChar> const & list)
{
    SEQAN_CHECKPOINT;
    return list.back();
}

//////////////////////////////////////////////////////////////////////////////
//assign to std::vector

template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar,  TAlloc> & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source, Generous());
}
template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar,  TAlloc> & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source, Generous());
}

template <typename TChar,  typename TAlloc, typename TSource, typename TSize>
inline void
assign(std::vector<TChar,  TAlloc> & target,
       TSource & source,
       TSize limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, Generous());
}
template <typename TChar,  typename TAlloc, typename TSource, typename TSize>
inline void
assign(std::vector<TChar,  TAlloc> & target,
       TSource const & source,
       TSize limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, Generous());
}

//____________________________________________________________________________

template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar, TAlloc> & target,
       TSource & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.assign(begin(source, Standard()), end(source, Standard()));
}
template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar, TAlloc> & target,
       TSource const & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.assign(begin(source, Standard()), end(source, Standard()));
}


template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign_std_vector_Generous_impl(std::vector<TChar,  TAlloc> & target,
                                TSource & source,
                                typename Size< std::vector<TChar,  TAlloc> >::Type limit)
{
    SEQAN_CHECKPOINT;
    typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
    typename Size<TSource const>::Type source_length = length(source);
    if (source_length > limit)
    {
        source_length = limit;
    }
    target.assign(source_begin, source_begin + source_length);
}
template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar,  TAlloc> & target,
       TSource & source,
       typename Size< std::vector<TChar, TAlloc> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    assign_std_vector_Generous_impl(target, source, limit);
}
template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size< std::vector<TChar, TAlloc> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    assign_std_vector_Generous_impl(target, source, limit);
}

//____________________________________________________________________________

template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar,  TAlloc> & target,
       TSource & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, target.capacity(), Generous());
}
template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar,  TAlloc> & target,
       TSource const & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar, TAlloc> & target,
       TSource & source,
       typename Size< std::vector<TChar, TAlloc> >::Type limit,
       Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    assign(target, source, limit, Generous());
}
template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(std::vector<TChar,  TAlloc> & target,
       TSource const & source,
       typename Size< std::vector<TChar,  TAlloc> >::Type limit,
       Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    assign(target, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
//append to std::vector

template <typename TChar, typename TAlloc, typename TSource>
inline void
append(std::vector<TChar,  TAlloc> & target,
       TSource const & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.insert(target.end(), begin(source, Standard()), end(source, Standard()));
}

template <typename TChar,  typename TAlloc, typename TSource>
inline void
append(std::vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size< std::vector<TChar, TAlloc> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    typename Size< std::vector<TChar, TAlloc> >::Type target_length = target.length();
    if (target_length > limit)
    {
        target.resize(limit);
    }
    else
    {
        limit -= target_length;
        typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
        typename Size<TSource const>::Type source_length = length(source);
        if (source_length > limit)
        {
            source_length = limit;
        }

        target.insert(target.end(), source_begin, source_begin + source_length);
    }
}

//____________________________________________________________________________

template <typename TChar, typename TAlloc, typename TSource>
inline void
append(std::vector<TChar,  TAlloc> & target,
       TSource const & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
append(std::vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size< std::vector<TChar, TAlloc> >::Type limit,
       Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    append(target, source, limit, Generous());
}

//////////////////////////////////////////////////////////////////////////////
template <typename TChar, typename TAlloc, typename TValue, typename TTag>
inline void
appendValue(std::vector<TChar, TAlloc> & me,
            TValue const & _value,
            TTag)
{
    SEQAN_CHECKPOINT;
    me.push_back(_value);
}

template <typename TChar, typename TAlloc, typename TValue>
inline void
appendValue(std::vector<TChar,  TAlloc> & me,
            TValue const & _value,
            Limit)
{
    SEQAN_CHECKPOINT;
    if (capacity(me) > length(me)) me.push_back(_value);
}

//////////////////////////////////////////////////////////////////////////////
//replace to std::vector

template <typename TChar,  typename TAlloc, typename TSource>
inline void
replace(std::vector<TChar, TAlloc> & target,
        typename Position< std::vector<TChar, TAlloc> >::Type pos_begin,
        typename Position< std::vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        Generous)
{
    SEQAN_CHECKPOINT;
    typename Size< std::vector<TChar, TAlloc> >::Type target_size = pos_end-pos_begin;
    typename Size< std::vector<TChar, TAlloc> >::Type source_size =length(source);

    if(target_size >= source_size)
        {
            copy(source.begin(),source.end(), target.begin()+pos_begin);
            if( target_size > source_size )
                target.erase(target.begin()+pos_begin+source_size,target.begin()+pos_end);
        }
    else
        {
            copy(source.begin(),source.begin()+target_size,target.begin()+pos_begin);
            target.insert(target.begin()+pos_end,source.begin()+target_size,source.end());
        }
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
replace(std::vector<TChar, TAlloc> & target,
        typename Position< std::vector<TChar, TAlloc> >::Type pos_begin,
        typename Position< std::vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        typename Size< std::vector<TChar, TAlloc> >::Type limit,
        Generous)
{
    SEQAN_CHECKPOINT;
    if (pos_begin >= limit)
    {
        target.resize(limit);
    }
    else
    {
        typename Size<TSource const>::Type source_length = length(source);
        typename Size< std::vector<TChar, TAlloc> >::Type pos_mid = pos_begin + source_length;
        typename Size< std::vector<TChar, TAlloc> >::Type pos_limit(limit);
        if (pos_mid > limit)
        {
            target.resize(limit);
            replace(target,pos_begin,pos_limit,source);
            target.resize(limit);
        }
        else
        {
            replace(target,pos_begin,pos_end,source);
            if (target.size() > limit)
            {
                target.resize(limit);
            }
        }
    }

}

template <typename TChar,  typename TAlloc, typename TSource>
inline void
replace(std::vector<TChar,  TAlloc> & target,
        typename Position< std::vector<TChar, TAlloc> >::Type pos_begin,
        typename Position< std::vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        Limit)
{
    SEQAN_CHECKPOINT;
    replace(target, pos_begin, pos_end, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
replace(std::vector<TChar, TAlloc> & target,
        typename Position< std::vector<TChar,  TAlloc> >::Type pos_begin,
        typename Position< std::vector<TChar,  TAlloc> >::Type pos_end,
        TSource const & source,
        typename Size< std::vector<TChar, TAlloc> >::Type limit,
        Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }
    replace(target, pos_begin, pos_end, source, limit, Generous());
}


//////////////////////////////////////////////////////////////////////////////
// handling of iterators as begin and end

template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand>
inline void
replace(std::vector<TChar, TAlloc> & target,
        typename Iterator< std::vector<TChar, TAlloc>, Rooted>::Type pos_begin,
        typename Iterator< std::vector<TChar, TAlloc>, Rooted>::Type pos_end,
        TSource & source,
        Tag<TExpand> tag)
{
    replace(target, position(pos_begin), position(pos_end), source, tag);
}

/*
template<typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
replace(std::vector<TChar, TAlloc> & target,
        typename Iterator< std::vector<TChar, TAlloc>, Rooted>::Type pos_begin,
        typename Iterator< std::vector<TChar, TAlloc>, Rooted>::Type pos_end,
        TSource & source,
        typename Size< std::vector<TChar, TAlloc> >::Type limit,
        Tag<TExpand> tag)
{
    replace(target,  position(pos_begin),  position(pos_end), source, tag);
}
*/


template <typename TChar,  typename TAlloc, typename TSize, typename TExpand>
inline typename Size< std::vector<TChar, TAlloc> >::Type
reserve(
    std::vector<TChar, TAlloc> & seq,
    TSize new_capacity,
    Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    seq.reserve(new_capacity);
    return _capacityReturned(seq, new_capacity, tag);
}

template <typename TChar, typename TAlloc, typename TSize>
inline typename Size< std::vector<TChar, TAlloc> >::Type
reserve(
    std::vector<TChar, TAlloc> & seq,
    TSize new_capacity,
    Insist const &)
{
    SEQAN_CHECKPOINT;
    // do nothing
    return _capacityReturned(seq, new_capacity, Insist());
}

template <typename TChar,  typename TAlloc, typename TSize>
inline typename Size< std::vector<TChar, TAlloc> >::Type
reserve(
    std::vector<TChar,  TAlloc> & seq,
    TSize new_capacity,
    Limit const &)
{
    SEQAN_CHECKPOINT;
    // do nothing
    return _capacityReturned(seq, new_capacity, Limit());
}

template <typename TChar,  typename TAlloc, typename TSize, typename TExpand>
inline typename Size< std::vector<TChar,  TAlloc> >::Type
resize(
    std::vector<TChar, TAlloc> & me,
    TSize new_length,
    Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    me.resize(new_length);
    return me.size();
}

template <typename TChar, typename TAlloc, typename TSize, typename TExpand>
inline typename Size< std::vector<TChar,  TAlloc> >::Type
fill(
    std::vector<TChar, TAlloc> & me,
    TSize new_length,
    TChar const & val,
    Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    me.resize(new_length, val);
    return me.length();
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_ADAPT_STD_VECTOR_H_
