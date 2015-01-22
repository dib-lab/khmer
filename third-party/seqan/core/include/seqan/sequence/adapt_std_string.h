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
// Adaptions for STL strings to SeqAn strings.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_ADAPT_STD_STRING_H_
#define SEQAN_SEQUENCE_ADAPT_STD_STRING_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Adaption."std::basic_string"
..summary:Adaption for STL strings.
 */

// ===========================================================================
// Metafunctions
// ===========================================================================

///.Metafunction.StdContainerIterator.param.T.type:Adaption."std::basic_string"
///.Metafunction.StdContainerIterator.class:Adaption."std::basic_string"

template <typename TChar, typename TCharTraits, typename TAlloc>
struct StdContainerIterator< std::basic_string<TChar, TCharTraits, TAlloc> >
{
    typedef std::basic_string<TChar, TCharTraits, TAlloc> TContainer;
    typedef typename TContainer::iterator Type;
};
    
template <typename TChar, typename TCharTraits, typename TAlloc>
struct StdContainerIterator< std::basic_string<TChar, TCharTraits, TAlloc> const>
{
    typedef std::basic_string<TChar, TCharTraits, TAlloc> TContainer;
    typedef typename TContainer::const_iterator Type;
};

///.Metafunction.IsContiguous.param.T.type:Adaption.std::basic_string
///.Metafunction.IsContiguous.class:Adaption."std::basic_string"

template <typename TChar, typename TCharTraits, typename TAlloc>
struct IsContiguous< std::basic_string<TChar, TCharTraits, TAlloc> >
{
    enum { VALUE = false };
};

template <typename  TChar, typename TCharTraits, typename TAlloc>
struct IsContiguous< std::basic_string<TChar, TCharTraits, TAlloc> const>
        : IsContiguous< std::basic_string<TChar, TCharTraits, TAlloc> > {};

///.Metafunction.Value.param.T.type:Adaption.std::basic_string
///.Metafunction.Value.class:Adaption."std::basic_string"

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Value< std::basic_string<TChar, TCharTraits, TAlloc> >
{
    typedef typename std::basic_string<TChar, TCharTraits, TAlloc>::value_type Type;
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Value< std::basic_string<TChar, TCharTraits, TAlloc> const>
        : Value< std::basic_string<TChar, TCharTraits, TAlloc> > {};

///.Metafunction.GetValue.param.T.type:Adaption.std::basic_string
///.Metafunction.GetValue.class:Adaption."std::basic_string"
// TODO(holtgrew): GetValue is a reference?! I thought the reverse was true in respect to Value<>.

template <typename TChar, typename TCharTraits, typename TAlloc>
struct GetValue< std::basic_string<TChar, TCharTraits, TAlloc> >
{
    typedef typename std::basic_string<TChar, TCharTraits, TAlloc>::reference Type;
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct GetValue< std::basic_string<TChar, TCharTraits, TAlloc> const>
{
    typedef typename std::basic_string<TChar, TCharTraits, TAlloc>::const_reference Type;
};

///.Metafunction.Reference.param.T.type:Adaption.std::basic_string
///.Metafunction.Reference.class:Adaption."std::basic_string"

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Reference< std::basic_string<TChar, TCharTraits, TAlloc> >
{
    typedef typename std::basic_string<TChar, TCharTraits, TAlloc>::reference Type;
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Reference< std::basic_string<TChar, TCharTraits, TAlloc> const>
{
    typedef typename std::basic_string<TChar, TCharTraits, TAlloc>::const_reference Type;
};

///.Metafunction.Iterator.param.T.type:Adaption.std::basic_string
///.Metafunction.Iterator.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>
{
    typedef std::basic_string<TChar, TCharTraits, TAlloc> TString_;
    typedef Iter<TString_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TString_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator< std::basic_string<TChar, TCharTraits, TAlloc> const, Rooted>
{
    typedef std::basic_string<TChar, TCharTraits, TAlloc> const TString_;
    typedef Iter<TString_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TString_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Standard >
{
    typedef Iter< std::basic_string<TChar, TCharTraits, TAlloc>, StdIteratorAdaptor > Type;
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Iterator< std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>
{
    typedef Iter< std::basic_string<TChar, TCharTraits, TAlloc> const, StdIteratorAdaptor > Type;
};

///.Metafunction.Position.param.T.type:Adaption.std::basic_string
///.Metafunction.Position.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Position< std::basic_string<TChar, TCharTraits, TAlloc> >
{
    typedef typename std::basic_string<TChar, TCharTraits, TAlloc>::size_type Type;
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Position< std::basic_string<TChar, TCharTraits, TAlloc> const>
        : Position< std::basic_string<TChar, TCharTraits, TAlloc> > {};

///.Metafunction.Size.param.T.type:Adaption.std::basic_string
///.Metafunction.Size.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Size< std::basic_string<TChar, TCharTraits, TAlloc> >
{
    typedef typename std::basic_string<TChar, TCharTraits, TAlloc>::size_type Type;
};

template <typename TChar, typename TCharTraits, typename TAlloc>
struct Size< std::basic_string<TChar, TCharTraits, TAlloc> const>
        : Size< std::basic_string<TChar, TCharTraits, TAlloc> > {};

///.Metafunction.Size.param.T.type:Adaption.std::basic_string
///.Metafunction.Size.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
struct DefaultOverflowImplicit< std::basic_string<TChar, TCharTraits, TAlloc> >
{
    typedef Generous Type;
};

// ===========================================================================
// Functions
// ===========================================================================

///.Function.getObjectId.param.object.type:Adaption.std::basic_string
///.Function.getObjectId.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline void const *
getObjectId(std::basic_string<TChar, TCharTraits, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    if (me.empty())
        return NULL;
    else
        return (& *(me.end() - 1)) + 1;
}

///.Function.begin.param.object.type:Adaption.std::basic_string
///.Function.begin.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type
begin(std::basic_string<TChar, TCharTraits, TAlloc> & me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type(me.begin());
}
template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type
begin(std::basic_string<TChar, TCharTraits, TAlloc> const & me,
      Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type(me.begin());
}

///.Function.end.param.object.type:Adaption.std::basic_string
///.Function.end.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type
end(std::basic_string<TChar, TCharTraits, TAlloc> & me,
    Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Standard>::Type(me.end());
}
template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type
end(std::basic_string<TChar, TCharTraits, TAlloc> const & me,
    Standard)
{
    SEQAN_CHECKPOINT;
    return typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc> const, Standard>::Type(me.end());
}

///.Function.value.param.container.type:Adaption.std::basic_string
///.Function.value.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TPos>
inline typename GetValue< std::basic_string<TChar, TCharTraits, TAlloc> >::Type
value(std::basic_string<TChar, TCharTraits, TAlloc> & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    return me[pos];
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TPos>
inline typename GetValue< std::basic_string<TChar, TCharTraits, TAlloc> const>::Type
value(std::basic_string<TChar, TCharTraits, TAlloc> const & me,
      TPos pos)
{
    SEQAN_CHECKPOINT;
    return me[pos];
}

///.Function.length.param.object.type:Adaption.std::basic_string
///.Function.length.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type
length(std::basic_string<TChar, TCharTraits, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    return me.length();
}

///.Function.capacity.param.object.type:Adaption.std::basic_string
///.Function.capacity.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type
capacity(std::basic_string<TChar, TCharTraits, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    return me.capacity();
}

///.Function.empty.param.object.type:Adaption.std::basic_string
///.Function.empty.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline bool
empty(std::basic_string<TChar, TCharTraits, TAlloc> const & me)
{
    SEQAN_CHECKPOINT;
    return me.empty();
}

///.Function.clear.param.object.type:Adaption.std::basic_string
///.Function.clear.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc>
inline void
clear(std::basic_string<TChar, TCharTraits, TAlloc> & me)
{
    SEQAN_CHECKPOINT;
    me.clear();
}

//////////////////////////////////////////////////////////////////////////////
//assign to std::basic_string

///.Function.assign.param.target.type:Adaption.std::basic_string
///.Function.assign.param.source.type:Adaption.std::basic_string
///.Function.assign.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source, Generous());
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source, Generous());
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TSize>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource & source,
       TSize limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, Generous());
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TSize>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source,
       TSize limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, limit, Generous());
}

//____________________________________________________________________________

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.assign(begin(source, Standard()), end(source, Standard()));
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.assign(begin(source, Standard()), end(source, Standard()));
}


template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign_std_string_Generous_impl(std::basic_string<TChar, TCharTraits, TAlloc> & target,
                                TSource & source,
                                typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit)
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
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource & source,
       typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    assign_std_string_Generous_impl(target, source, limit);
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source,
       typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    assign_std_string_Generous_impl(target, source, limit);
}

//____________________________________________________________________________

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, target.capacity(), Generous());
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    assign(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource & source,
       typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
       Limit)
{
    SEQAN_CHECKPOINT;
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    assign(target, source, limit, Generous());
}
template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
assign(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source,
       typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
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
//append to std::basic_string

///.Function.append.param.target.type:Adaption.std::basic_string
///.Function.append.param.source.type:Adaption.std::basic_string
///.Function.append.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
append(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source,
       Generous)
{
    SEQAN_CHECKPOINT;
    target.append(begin(source, Standard()), end(source, Standard()));
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
append(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source,
       typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
       Generous)
{
    SEQAN_CHECKPOINT;
    typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type target_length = target.length();
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
        target.append(source_begin, source_begin + source_length);
    }
}

//____________________________________________________________________________

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
append(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source,
       Limit)
{
    SEQAN_CHECKPOINT;
    append(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
append(std::basic_string<TChar, TCharTraits, TAlloc> & target,
       TSource const & source,
       typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
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
///.Function.appendValue.param.target.type:Adaption.std::basic_string
///.Function.appendValue.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TValue, typename TTag>
inline void
appendValue(std::basic_string<TChar, TCharTraits, TAlloc> & me,
            TValue const & _value,
            TTag)
{
    SEQAN_CHECKPOINT;
    me.push_back(_value);
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TValue>
inline void
appendValue(std::basic_string<TChar, TCharTraits, TAlloc> & me,
            TValue const & _value,
            Limit)
{
    SEQAN_CHECKPOINT;
    if (capacity(me) > length(me)) me.push_back(_value);
}

//////////////////////////////////////////////////////////////////////////////
//replace to std::basic_string

///.Function.replace.param.target.type:Adaption.std::basic_string
///.Function.replace.param.source.type:Adaption.std::basic_string
///.Function.replace.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
replace(std::basic_string<TChar, TCharTraits, TAlloc> & target,
        typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin,
        typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end,
        TSource const & source,
        Generous)
{
    SEQAN_CHECKPOINT;
    target.replace(target.begin() + pos_begin, target.begin() + pos_end, begin(source, Standard()), end(source, Standard()));
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
replace(std::basic_string<TChar, TCharTraits, TAlloc> & target,
        typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin,
        typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end,
        TSource const & source,
        typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
        Generous)
{
    SEQAN_CHECKPOINT;
    if (pos_begin >= limit)
    {
        target.resize(limit);
    }
    else
    {
        typename Iterator<TSource const, Standard>::Type source_begin = begin(source, Standard());
        typename Size<TSource const>::Type source_length = length(source);
        typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_mid = pos_begin + source_length;
        if (pos_mid > limit)
        {
            target.replace(target.begin() + pos_begin, target.begin() + limit, source_begin, source_begin + limit - pos_begin);
            target.resize(limit);
        }
        else
        {
            target.replace(target.begin() + pos_begin, target.begin() + pos_end, source_begin, end(source, Standard()));
            if (target.length() > limit)
            {
                target.resize(limit);
            }
        }
    }
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
replace(std::basic_string<TChar, TCharTraits, TAlloc> & target,
        typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin,
        typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end,
        TSource const & source,
        Limit)
{
    SEQAN_CHECKPOINT;
    replace(target, pos_begin, pos_end, source, target.capacity(), Generous());
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSource>
inline void
replace(std::basic_string<TChar, TCharTraits, TAlloc> & target,
        typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_begin,
        typename Position< std::basic_string<TChar, TCharTraits, TAlloc> >::Type pos_end,
        TSource const & source,
        typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
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

/*
template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand>
inline void
replace(std::basic_string<TChar, TCharTraits, TAlloc> & target,
        typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_begin,
        typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_end,
        TSource & source,
        Tag<TExpand> tag)
{
    replace(target, position(pos_begin), position(pos_end), source, tag);
}

template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand>
inline void
replace(std::basic_string<TChar, TCharTraits, TAlloc> & target,
        typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_begin,
        typename Iterator< std::basic_string<TChar, TCharTraits, TAlloc>, Rooted>::Type pos_end,
        TSource & source,
        typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type limit,
        Tag<TExpand> tag)
{
    replace(target,  position(pos_begin),  position(pos_end), source, tag);
}
*/

///.Function.reserve.param.object.type:Adaption.std::basic_string
///.Function.reserve.remarks:For @Adaption.std::basic_string|STL Adaptions@, $reserve$ is only guaranteed to have the specified behaviour with $Insist$ and $Generous$.
///.Function.reserve.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize, typename TExpand>
inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type
reserve(
    std::basic_string<TChar, TCharTraits, TAlloc> & seq,
    TSize new_capacity,
    Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    seq.reserve(new_capacity);
    return _capacityReturned(seq, new_capacity, tag);
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize>
inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type
reserve(
    std::basic_string<TChar, TCharTraits, TAlloc> & seq,
    TSize new_capacity,
    Insist const &)
{
    SEQAN_CHECKPOINT;
    // do nothing
    return _capacityReturned(seq, new_capacity, Insist());
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize>
inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type
reserve(
    std::basic_string<TChar, TCharTraits, TAlloc> & seq,
    TSize new_capacity,
    Limit const &)
{
    SEQAN_CHECKPOINT;
    // do nothing
    return _capacityReturned(seq, new_capacity, Limit());
}

///.Function.resize.param.object.type:Adaption.std::basic_string
///.Function.resize.class:Adaption.std::basic_string

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize, typename TExpand>
inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type
resize(
    std::basic_string<TChar, TCharTraits, TAlloc> & me,
    TSize new_length,
    Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    me.resize(new_length);
    return me.length();
}

template <typename TChar, typename TCharTraits, typename TAlloc, typename TSize, typename TExpand>
inline typename Size< std::basic_string<TChar, TCharTraits, TAlloc> >::Type
resize(
    std::basic_string<TChar, TCharTraits, TAlloc> & me,
    TSize new_length,
    TChar const & val,
    Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    me.resize(new_length, val);
    return me.length();
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_ADAPT_STD_STRING_H_
