// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Adaptions for Thrust device_vectors to SeqAn strings.
// ==========================================================================
// NOTE(esiragusa): This adaption is not tested. Anything beside assign might
// not work as expected.

#ifndef SEQAN_SEQUENCE_ADAPT_THRUST_VECTOR_H_
#define SEQAN_SEQUENCE_ADAPT_THRUST_VECTOR_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ===========================================================================
// Metafunctions
// ===========================================================================

// ----------------------------------------------------------------------------
// Metafunction IsSequence
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct IsSequence<thrust::device_vector<TChar, TAlloc> >
{
    typedef True Type;
    enum { VALUE = true };
};

// ----------------------------------------------------------------------------
// Metafunction IsContiguous
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct IsContiguous<thrust::device_vector<TChar, TAlloc> >
{
    enum { VALUE = true };
};

template <typename  TChar, typename TAlloc>
struct IsContiguous<thrust::device_vector<TChar, TAlloc> const>
        : IsContiguous<thrust::device_vector<TChar, TAlloc> > {};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct Value<thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar, TAlloc>::value_type Type;
};

template <typename TChar, typename TAlloc>
struct Value<thrust::device_vector<TChar, TAlloc> const>
        : Value<thrust::device_vector<TChar, TAlloc> > {};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct Reference<thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar, TAlloc>::reference Type;
};

template <typename TChar,  typename TAlloc>
struct Reference<thrust::device_vector<TChar, TAlloc> const>
{
    typedef typename thrust::device_vector<TChar, TAlloc>::const_reference Type;
};

template <typename TValue>
struct Reference<thrust::detail::normal_iterator<thrust::device_ptr<TValue> > >
{
    typedef thrust::device_reference<TValue> Type;
};

template <typename TValue>
struct Reference<thrust::detail::normal_iterator<thrust::device_ptr<TValue> > const>
{
    typedef thrust::device_reference<TValue> const Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct GetValue<thrust::device_vector<TChar, TAlloc> > : Reference<thrust::device_vector<TChar, TAlloc> const> {};

template <typename TChar, typename TAlloc>
struct GetValue<thrust::device_vector<TChar, TAlloc> const> : Reference<thrust::device_vector<TChar, TAlloc> const> {};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct Iterator<thrust::device_vector<TChar, TAlloc>, Rooted>
{
    typedef thrust::device_vector<TChar, TAlloc> TVector_;
    typedef Iter<TVector_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TVector_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar, typename TAlloc>
struct Iterator<thrust::device_vector<TChar, TAlloc> const, Rooted>
{
    typedef thrust::device_vector<TChar, TAlloc> const TVector_;
    typedef Iter<TVector_, StdIteratorAdaptor> TIterator_;
    typedef Iter<TVector_, AdaptorIterator<TIterator_> > Type;
};

template <typename TChar,  typename TAlloc>
struct Iterator<thrust::device_vector<TChar, TAlloc>, Standard>
{
    typedef typename thrust::device_vector<TChar, TAlloc>::iterator Type;
};

template <typename TChar,  typename TAlloc>
struct Iterator<thrust::device_vector<TChar, TAlloc> const, Standard>
{
    typedef typename thrust::device_vector<TChar, TAlloc>::const_iterator Type;
};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc>
struct Position<thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar, TAlloc>::size_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc>
struct Size<thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar, TAlloc>::size_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct Difference<thrust::device_vector<TChar, TAlloc> >
{
    typedef typename thrust::device_vector<TChar, TAlloc>::difference_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct DefaultOverflowImplicit<thrust::device_vector<TChar, TAlloc> >
{
    typedef Generous Type;
};

// ----------------------------------------------------------------------------
// Metafunction StdContainerIterator
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
struct StdContainerIterator<thrust::device_vector<TChar, TAlloc> >
{
    typedef thrust::device_vector<TChar, TAlloc> TContainer_;
    typedef typename TContainer_::iterator Type;
};

template <typename TChar, typename TAlloc>
struct StdContainerIterator<thrust::device_vector<TChar, TAlloc> const>
{
    typedef thrust::device_vector<TChar, TAlloc> TContainer_;
    typedef typename TContainer_::const_iterator Type;
};

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc>
struct View<thrust::device_vector<TValue, TAlloc> >
{
    typedef ContainerView<thrust::device_vector<TValue, TAlloc> >   Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsDevice
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc>
struct IsDevice<thrust::device_vector<TValue, TAlloc> > : public True {};

// ----------------------------------------------------------------------------
// Metafunction Device
// ----------------------------------------------------------------------------

template <typename TValue, typename TAlloc>
struct Device<String<TValue, TAlloc> >
{
    typedef thrust::device_vector<TValue>   Type;
};

template <typename TValue>
struct Device<std::vector<TValue> >
{
    typedef thrust::device_vector<TValue>   Type;
};

template <>
struct Device<std::string>
{
    typedef thrust::device_vector<char>     Type;
};

// ===========================================================================
// Functions
// ===========================================================================

// ----------------------------------------------------------------------------
// Function getObjectId()
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc>
inline void const *
getObjectId(thrust::device_vector<TChar, TAlloc> const & me)
{
    if (me.empty())
        return NULL;
    else
        return (& *(me.end() - 1)) + 1;
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc>
inline typename Iterator<thrust::device_vector<TChar, TAlloc>, Standard>::Type
begin(thrust::device_vector<TChar, TAlloc> & me, Standard)
{
    return typename Iterator<thrust::device_vector<TChar, TAlloc>, Standard>::Type(me.begin());
}

template <typename TChar,  typename TAlloc>
inline typename Iterator<thrust::device_vector<TChar, TAlloc> const, Standard>::Type
begin(thrust::device_vector<TChar, TAlloc> const & me, Standard)
{
    return typename Iterator<thrust::device_vector<TChar, TAlloc> const, Standard>::Type(me.begin());
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
inline typename Iterator<thrust::device_vector<TChar, TAlloc>, Standard>::Type
end(thrust::device_vector<TChar, TAlloc> & me, Standard)
{
    return typename Iterator<thrust::device_vector<TChar, TAlloc>, Standard>::Type(me.end());
}

template <typename TChar,  typename TAlloc>
inline typename Iterator<thrust::device_vector<TChar, TAlloc> const, Standard>::Type
end(thrust::device_vector<TChar, TAlloc> const & me, Standard)
{
    return typename Iterator<thrust::device_vector<TChar, TAlloc> const, Standard>::Type(me.end());
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc, typename TPos>
inline typename Reference<thrust::device_vector<TChar, TAlloc> >::Type
value(thrust::device_vector<TChar, TAlloc> & me, TPos pos)
{
    return me[pos];
}

template <typename TChar,  typename TAlloc, typename TPos>
inline typename Reference<thrust::device_vector<TChar, TAlloc> const>::Type
value(thrust::device_vector<TChar, TAlloc> const & me, TPos pos)
{
    return me[pos];
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type
length(thrust::device_vector<TChar, TAlloc> const & me)
{
    return me.size();
}

// ----------------------------------------------------------------------------
// Function capacity()
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type
capacity(thrust::device_vector<TChar, TAlloc> const & me)
{
    return me.capacity();
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
inline bool
empty(thrust::device_vector<TChar, TAlloc> const & me)
{
    return me.empty();
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc>
inline void
clear(thrust::device_vector<TChar, TAlloc> & me)
{
    me.clear();
}

// ----------------------------------------------------------------------------
// Function front()
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
inline typename Reference<thrust::device_vector<TChar, TAlloc> >::Type
front(thrust::device_vector<TChar, TAlloc> & list)
{
    return list.front();
}

template <typename TChar, typename TAlloc>
inline typename Reference<thrust::device_vector<TChar, TAlloc> const>::Type
front(thrust::device_vector<TChar, TAlloc> const & list)
{
    return list.front();
}

// ----------------------------------------------------------------------------
// Function back()
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc>
inline typename Reference<thrust::device_vector<TChar, TAlloc> >::Type
back(thrust::device_vector<TChar, TAlloc> & list)
{
    return list.back();
}

template <typename TChar, typename TAlloc>
inline typename Reference<thrust::device_vector<TChar, TAlloc> const>::Type
back(thrust::device_vector<TChar, TAlloc> const & list)
{
    return list.back();
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source)
{
    assign(target, source, Generous());
}

template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source)
{
    assign(target, reinterpret_cast<TSource const &>(source), Generous());
}

template <typename TChar,  typename TAlloc, typename TSource, typename TSize>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, TSize limit)
{
    assign(target, source, limit, Generous());
}

template <typename TChar,  typename TAlloc, typename TSource, typename TSize>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source, TSize limit)
{
    assign(target, reinterpret_cast<TSource const &>(source), limit, Generous());
}

// ----------------------------------------------------------------------------
// Function assign(); Generous
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, Generous)
{
    target.assign(begin(source, Standard()), end(source, Standard()));
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source, Generous)
{
    assign(target, reinterpret_cast<TSource const &>(source), Generous());
}

template <typename TChar,  typename TAlloc, typename TSource>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit,
       Generous)
{
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
assign(thrust::device_vector<TChar, TAlloc> & target,
       TSource & source,
       typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit,
       Generous)
{
    assign(target, reinterpret_cast<TSource const &>(source), limit, Generous());
}

// ----------------------------------------------------------------------------
// Function assign(); Limit
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, Limit)
{
    assign(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target, TSource & source, Limit)
{
    assign(target, reinterpret_cast<TSource const &>(source), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit,
       Limit)
{
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    assign(target, source, limit, Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
assign(thrust::device_vector<TChar, TAlloc> & target,
       TSource & source,
       typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit,
       Limit)
{
    assign(target, reinterpret_cast<TSource const &>(source), limit, Generous());
}

// ----------------------------------------------------------------------------
// Function append(); Generous
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TSource>
inline void
append(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, Generous)
{
    target.insert(target.end(), begin(source, Standard()), end(source, Standard()));
}

template <typename TChar,  typename TAlloc, typename TSource>
inline void
append(thrust::device_vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit,
       Generous)
{
    typename Size<thrust::device_vector<TChar, TAlloc> >::Type target_length = target.length();
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

// ----------------------------------------------------------------------------
// Function append(); Limit
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TSource>
inline void
append(thrust::device_vector<TChar, TAlloc> & target, TSource const & source, Limit)
{
    append(target, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
append(thrust::device_vector<TChar, TAlloc> & target,
       TSource const & source,
       typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit,
       Limit)
{
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }

    append(target, source, limit, Generous());
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TValue, typename TTag>
inline void
appendValue(thrust::device_vector<TChar, TAlloc> & me, TValue const & _value, TTag)
{
    me.push_back(_value);
}

template <typename TChar, typename TAlloc, typename TValue>
inline void
appendValue(thrust::device_vector<TChar, TAlloc> & me, TValue const & _value, Limit)
{
    if (capacity(me) > length(me)) me.push_back(_value);
}

// ----------------------------------------------------------------------------
// Function replace()
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc, typename TSource>
inline void
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_begin,
        typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        Generous)
{
    typename Size<thrust::device_vector<TChar, TAlloc> >::Type target_size = pos_end-pos_begin;
    typename Size<thrust::device_vector<TChar, TAlloc> >::Type source_size =length(source);

    if (target_size >= source_size)
    {
        copy(source.begin(),source.end(), target.begin()+pos_begin);
        if (target_size > source_size)
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
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_begin,
        typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit,
        Generous)
{
    if (pos_begin >= limit)
    {
        target.resize(limit);
    }
    else
    {
        typename Size<TSource const>::Type source_length = length(source);
        typename Size<thrust::device_vector<TChar, TAlloc> >::Type pos_mid = pos_begin + source_length;
        typename Size<thrust::device_vector<TChar, TAlloc> >::Type pos_limit(limit);
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
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_begin,
        typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        Limit)
{
    replace(target, pos_begin, pos_end, source, target.capacity(), Generous());
}

template <typename TChar, typename TAlloc, typename TSource>
inline void
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_begin,
        typename Position<thrust::device_vector<TChar, TAlloc> >::Type pos_end,
        TSource const & source,
        typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit,
        Limit)
{
    if (limit > target.capacity())
    {
        limit = target.capacity();
    }
    replace(target, pos_begin, pos_end, source, limit, Generous());
}


// Handling of iterators as begin and end.

template<typename TChar, typename TCharTraits, typename TAlloc, typename TSource, typename TExpand>
inline void
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Iterator<thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_begin,
        typename Iterator<thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_end,
        TSource & source,
        Tag<TExpand> const tag)
{
    replace(target, position(pos_begin), position(pos_end), source, tag);
}

/*
template<typename TChar, typename TAlloc, typename TSource, typename TExpand>
inline void
replace(thrust::device_vector<TChar, TAlloc> & target,
        typename Iterator<thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_begin,
        typename Iterator<thrust::device_vector<TChar, TAlloc>, Rooted>::Type pos_end,
        TSource & source,
        typename Size<thrust::device_vector<TChar, TAlloc> >::Type limit,
        Tag<TExpand> const tag)
{
    replace(target,  position(pos_begin),  position(pos_end), source, tag);
}
*/

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc, typename TSize, typename TExpand>
inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type
reserve(thrust::device_vector<TChar, TAlloc> & seq, TSize new_capacity, Tag<TExpand> const & tag)
{
    seq.reserve(new_capacity);
    return _capacityReturned(seq, new_capacity, tag);
}

template <typename TChar, typename TAlloc, typename TSize>
inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type
reserve(thrust::device_vector<TChar, TAlloc> & seq, TSize new_capacity, Insist const &)
{
    // do nothing
    return _capacityReturned(seq, new_capacity, Insist());
}

template <typename TChar,  typename TAlloc, typename TSize>
inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type
reserve(thrust::device_vector<TChar, TAlloc> & seq, TSize new_capacity, Limit const &)
{
    // do nothing
    return _capacityReturned(seq, new_capacity, Limit());
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TChar,  typename TAlloc, typename TSize, typename TExpand>
inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type
resize(thrust::device_vector<TChar, TAlloc> & me, TSize new_length, Tag<TExpand> const &)
{
    me.resize(new_length);
    return me.size();
}

// ----------------------------------------------------------------------------
// Function fill()
// ----------------------------------------------------------------------------

template <typename TChar, typename TAlloc, typename TSize, typename TExpand>
inline typename Size<thrust::device_vector<TChar, TAlloc> >::Type
fill(thrust::device_vector<TChar, TAlloc> & me, TSize new_length, TChar const & val, Tag<TExpand> const &)
{
    me.resize(new_length, val);
    return me.length();
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TAlloc>
inline typename View<thrust::device_vector<TContainer, TAlloc> >::Type
view(thrust::device_vector<TContainer, TAlloc> & container)
{
    typedef typename View<thrust::device_vector<TContainer, TAlloc> >::Type TView;

    if (empty(container)) return TView();

    return TView(thrust::raw_pointer_cast(&container.front()),
                 thrust::raw_pointer_cast(&container.front()) + container.size());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_ADAPT_THRUST_VECTOR_H_
