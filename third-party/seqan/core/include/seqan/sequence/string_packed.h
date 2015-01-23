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

template <typename TPackedContainer>
struct PackedConsts_;

/**
.Spec.Packed String:
..cat:Strings
..general:Class.String
..summary:A string that stores as many values in one machine word as possible.
..signature:String<TValue, Packed<THostspec> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.THostspec:The specializing type.
...remarks:This is the specialization of the host string that is used for storing the packed values.
...default:@Spec.Alloc String.Alloc<>@
..include:seqan/sequence.h
*/

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

    THost data_host;
    TSize data_length;

    String():
        data_length(0)
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TSource>
    String(TSource & source):
        data_length(0)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source);
    }
    template <typename TSource>
    String(TSource const & source):
        data_length(0)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source);
    }
    String(String const & source):
        data_length(0)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source);
    }

    template <typename TSource>
    String & operator =(TSource const & source)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source);
        return *this;
    }
    String & operator =(String const & source)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source);
        return *this;
    }

    ~String()
    {
        SEQAN_CHECKPOINT;
    }

    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    template <typename TPos>
    inline typename Reference<String>::Type
    operator[](TPos pos)
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<String const>::Type 
    operator[](TPos pos) const
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }
};

// --------------------------------------------------------------------------
// Specialization Packed String Iter
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
class Iter<TContainer, Packed<THostspec> >
{
public:
    typedef typename HostIterator<Iter>::Type THostIterator;
    typedef typename Position<TContainer>::Type TPosition;

    typename Pointer_<TContainer>::Type data_container;
    THostIterator data_iterator;
    unsigned char data_bitpos;

    Iter()
    {
        SEQAN_CHECKPOINT;
    }

    Iter(typename Parameter_<TContainer>::Type container_)
        : data_container(_toPointer(container_)),
          data_iterator(begin(host(container_))),
          data_bitpos(0)
    {
        SEQAN_CHECKPOINT;
    }

    Iter(typename Parameter_<TContainer>::Type container_, TPosition pos_):
        data_container(_toPointer(container_))
    {
        SEQAN_CHECKPOINT;
        setPosition(*this, pos_);
    }

    Iter(Iter const & other_):
        data_container(other_.data_container),
        data_iterator(other_.data_iterator),
        data_bitpos(other_.data_bitpos)
    {
        SEQAN_CHECKPOINT;
    }

    ~Iter()
    {
        SEQAN_CHECKPOINT;
    }

    inline
    Iter const & 
    operator=(Iter const & other_)
    {
        SEQAN_CHECKPOINT;
        data_container = other_.data_container;
        data_iterator = other_.data_iterator;
        data_bitpos = other_.data_bitpos;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

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
// Metafunction Host
// --------------------------------------------------------------------------

///.Metafunction.Host.param.T.type:Spec.Packed String
template <typename TValue, typename THostspec>
struct Host<String<TValue, Packed<THostspec> > >
{
    typedef String<unsigned int, THostspec> Type;
};

template <typename TValue, typename THostspec>
struct Host<String<TValue, Packed<THostspec> > const>
{
    typedef String<unsigned int, THostspec> const Type;
};

// --------------------------------------------------------------------------
// Metafunction GetValue
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
struct GetValue<String<TValue, Packed<THostspec> > >
        : Value<String<TValue, Packed<THostspec> > >
{};

template <typename TValue, typename THostspec>
struct GetValue<String<TValue, Packed<THostspec> > const>
        : Value<String<TValue, Packed<THostspec> > const>
{};

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
struct Reference<String<TValue, Packed<THostspec> > const>
{
    typedef typename Iterator<String<TValue, Packed<THostspec> > const, Standard>::Type TIterator;
    typedef Proxy<IteratorProxy<TIterator> > Type;
};

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

template <typename TValue, typename THostspec, typename TSpec>
struct Iterator<String<TValue, Packed<THostspec> >, TSpec>
{
    typedef Iter<String<TValue, Packed<THostspec> >, Packed<THostspec> > Type;
};

template <typename TValue, typename THostspec, typename TSpec>
struct Iterator<String<TValue, Packed<THostspec> > const, TSpec>
{
    typedef Iter<String<TValue, Packed<THostspec> > const, Packed<THostspec> > Type;
};

// --------------------------------------------------------------------------
// Metafunction HostIterator
// --------------------------------------------------------------------------

// TODO(holtgrew): Actually is internal, mark so here and rename to HostIterator_
template <typename T>
struct HostIterator;

template <typename TContainer, typename THostspec>
struct HostIterator<Iter<TContainer, Packed<THostspec> > >
{
    typedef typename Host<TContainer>::Type THost_;
    typedef typename Iterator<THost_, Standard>::Type Type;
};

template <typename TContainer, typename THostspec>
struct HostIterator<Iter<TContainer, Packed<THostspec> > const>
{
    typedef typename Host<TContainer>::Type THost_;
    typedef typename Iterator<THost_, Standard>::Type const Type;
};

// --------------------------------------------------------------------------
// Metafunction PackedConsts_
// --------------------------------------------------------------------------

template <typename TPackedContainer>
struct PackedConsts_
{
    typedef typename Value<TPackedContainer>::Type TValue;
    typedef typename Host<TPackedContainer>::Type THost;
    typedef typename Value<THost>::Type THostValue;

    enum
    {
        BITS_PER_VALUE = BitsPerValue<TValue>::VALUE,
        BITS_PER_HOST_VALUE = BitsPerValue<THostValue>::VALUE,
        VALUES_PER_WORD = (BITS_PER_VALUE > BITS_PER_HOST_VALUE) ? 1 : (BITS_PER_HOST_VALUE / BITS_PER_VALUE),
        VALUE_MASK = (1 << BITS_PER_VALUE) - 1,
        MAX_BIT_POS = (VALUES_PER_WORD - 1) * BITS_PER_VALUE
    };

    static typename Size<THost>::Type
    toHostLength(typename Size<TPackedContainer>::Type len)
    {
        return (len + VALUES_PER_WORD - 1) / VALUES_PER_WORD;
    }
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

// --------------------------------------------------------------------------
// Function host
// --------------------------------------------------------------------------

///.Function.host.param.object.type:Spec.Packed String
///.Function.host.class:Spec.Packed String

template <typename TValue, typename THostspec>
inline typename Host<String<TValue, Packed<THostspec> > >::Type &
host(String<TValue, Packed<THostspec> > & me)
{
    SEQAN_CHECKPOINT;
    return me.data_host;
}

template <typename TValue, typename THostspec>
inline typename Host<String<TValue, Packed<THostspec> > const>::Type const &
host(String<TValue, Packed<THostspec> > const & me)
{
    SEQAN_CHECKPOINT;
    return me.data_host;
}

// --------------------------------------------------------------------------
// Function length
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline typename Size<String<TValue, Packed<THostspec> > >::Type
length(String<TValue, Packed<THostspec> > & me) 
{
    SEQAN_CHECKPOINT;
    return me.data_length;
}

template <typename TValue, typename THostspec>
inline typename Size<String<TValue, Packed<THostspec> > const>::Type
length(String<TValue, Packed<THostspec> > const & me) 
{
    SEQAN_CHECKPOINT;
    return me.data_length;
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
    SEQAN_CHECKPOINT;
    typedef String<TValue, Packed<THostspec> > TString;
    me.data_length = new_length;
    _setLength(host(me), PackedConsts_<TString>::toHostLength(new_length));
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
    typedef typename Size<TTarget>::Type TSize2;

    assign(host(target), host(source), tag);
    TSize2 new_length_limit = length(host(target)) * PackedConsts_<TTarget>::VALUES_PER_WORD;
    TSize2 new_length = length(source);
    if (new_length > new_length_limit)
    {
        new_length = new_length_limit;
    }
    _setLength(target, new_length);
}

template <typename TTarget, typename TSource, typename TSize, typename TTag>
inline void 
_assignCopyPackedString(TTarget & target,
                        TSource & source,
                        TSize limit,
                        Tag<TTag> const & tag)
{
    typedef typename Size<TTarget>::Type TSize2;

    TSize2 host_limit = PackedConsts_<TTarget>::toHostLength(limit);
    assign(host(target), host(source), host_limit, tag);
    TSize2 new_length_limit = length(host(target)) * PackedConsts_<TTarget>::VALUES_PER_WORD;
    TSize2 new_length = length(source);
    if (new_length > new_length_limit)
    {
        new_length = new_length_limit;
    }
    if (new_length > limit)
    {
        new_length = limit;
    }
    _setLength(target, new_length);
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
    SEQAN_CHECKPOINT;
    return getObjectId(host(me));
}

// --------------------------------------------------------------------------
// Function iter()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TPos, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type 
iter(String<TValue, Packed<THostspec> > & me,
     TPos pos_,
     Tag<TTag> const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type TIterator;
    return TIterator(me, pos_);
}

template <typename TValue, typename THostspec, typename TPos, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type 
iter(String<TValue, Packed<THostspec> > const & me,
     TPos pos_,
     Tag<TTag> const &)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type TIterator;
    return TIterator(me, pos_);
}

// --------------------------------------------------------------------------
// Function begin()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type 
begin(String<TValue, Packed<THostspec> > & me,
      Tag<TTag> const & tag_)
{
    SEQAN_CHECKPOINT;
    return iter(me, 0, tag_);
}

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type 
begin(String<TValue, Packed<THostspec> > const & me,
      Tag<TTag> const & tag_)
{
    SEQAN_CHECKPOINT;
    return iter(me, 0, tag_);
}

// --------------------------------------------------------------------------
// Function end()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type 
end(String<TValue, Packed<THostspec> > & me,
    Tag<TTag> const & tag_)
{
    SEQAN_CHECKPOINT;
    return iter(me, length(me), tag_);
}

template <typename TValue, typename THostspec, typename TTag>
inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type 
end(String<TValue, Packed<THostspec> > const & me,
    Tag<TTag> const & tag_)
{
    SEQAN_CHECKPOINT;
    return iter(me, length(me), tag_);
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec, typename TPos>
inline typename Reference<String<TValue, Packed<THostspec> > >::Type
value(String<TValue, Packed<THostspec> > & me, 
      TPos pos)
{
    SEQAN_CHECKPOINT;
    
    return *iter(me, pos, Standard());
} 

template <typename TValue, typename THostspec, typename TPos>
inline typename Reference<String<TValue, Packed<THostspec> > const>::Type
value(String<TValue, Packed<THostspec> > const & me, 
      TPos pos)
{
    SEQAN_CHECKPOINT;
    
    return *iter(me, pos, Standard());
} 

// --------------------------------------------------------------------------
// Function capacity()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline typename Size<String<TValue, Packed<THostspec> > const>::Type
capacity(String<TValue, Packed<THostspec> > const & me)
{
    SEQAN_CHECKPOINT;
    typedef typename Size<String<TValue, Packed<THostspec> > const>::Type TSize;
    TSize len = capacity(host(me));
    len *= PackedConsts_<String<TValue, Packed<THostspec> > >::VALUES_PER_WORD;
    return len;
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename TValue, typename THostspec>
inline void 
clear(String<TValue, Packed<THostspec> > & me)
{
    SEQAN_CHECKPOINT;
    clear(host(me));
    _setLength(me, 0);
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
    SEQAN_CHECKPOINT;
        typedef typename Size<T>::Type TSize;
        TSize wanted_host_length = PackedConsts_<T>::toHostLength(size);
        TSize new_host_length = resize(host(seq), wanted_host_length, TExpand());
        if (new_host_length < wanted_host_length)
        {
            size = new_host_length * PackedConsts_<T>::VALUES_PER_WORD;
        }
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
        if (limit < size)
        {
    SEQAN_CHECKPOINT;
            size = limit;
        }
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
    SEQAN_CHECKPOINT;
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
    SEQAN_CHECKPOINT;
//??? TODO: This function can be accelerated this way: 
//              - move values in host
//              - avoid double moving of the rest-part if "resize" allocates a new block

        typedef typename Size<T>::Type TSize;

        TSize old_length = length(seq);
        TSize old_size = end - start;
        TSize wanted_new_length = old_length + size - old_size;

        if (wanted_new_length > limit)
        {
            wanted_new_length = limit;
        }

        TSize wanted_host_length = PackedConsts_<T>::toHostLength(wanted_new_length);
        TSize new_host_length = resize(host(seq), wanted_host_length, TExpand());

        TSize new_length;
        if (new_host_length < wanted_host_length)
        {
            new_length = new_host_length * PackedConsts_<T>::VALUES_PER_WORD;
            if (new_length <= start + size)
            {
                goto FINISH;
            }
            old_length = new_length - size + old_size;
        }
        else
        {
            new_length = wanted_new_length;
        }
/*
        //move [end:right_end] to [start + size:..]
        if (old_size > size)
        {//move rest to left
            ::std::copy(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq, end + size - old_size, Standard()));
        }
        else
        {//move rest to right
            ::std::copy_backward(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq,  new_length, Standard()));
        }
*/
        if (old_size > size)
        {
            arrayMoveForward(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq, end + size - old_size, Standard()));
        }
        else
        {
            arrayMoveBackward(iter(seq, end, Standard()), iter(seq, old_length, Standard()), iter(seq, end + size - old_size, Standard()));
        }
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
    SEQAN_CHECKPOINT;
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
    SEQAN_CHECKPOINT;
        typename Iterator<T>::Type seq_begin = begin(seq);
        return _clearSpace(seq, size, start - seq_begin, end - seq_begin, limit, Insist());
    }
*/
};

template<typename TValue, typename THostspec, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type 
_clearSpace(String<TValue, Packed<THostspec> > & me, 
        typename Size< String<TValue, Packed<THostspec> > >::Type size, 
        Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size);
}

template<typename TValue, typename THostspec, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type 
_clearSpace(String<TValue, Packed<THostspec> > & me, 
        typename Size< String<TValue, Packed<THostspec> > >::Type size, 
        typename Size< String<TValue, Packed<THostspec> > >::Type limit, 
        Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size, limit);
}

template<typename TValue, typename THostspec, typename TPosition, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type 
_clearSpace(String<TValue, Packed<THostspec> > & me, 
            typename Size< String<TValue, Packed<THostspec> > >::Type size, 
            TPosition pos_begin, 
            TPosition pos_end, 
            Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end);
}

template<typename TValue, typename THostspec, typename TPosition, typename TExpand>
inline typename Size< String<TValue, Packed<THostspec> > >::Type 
_clearSpace(String<TValue, Packed<THostspec> > & me, 
            typename Size< String<TValue, Packed<THostspec> > >::Type size, 
            TPosition pos_begin, 
            TPosition pos_end, 
            typename Size< String<TValue, Packed<THostspec> > >::Type limit, 
            Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    return ClearSpaceStringPacked_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end, limit);
}

// --------------------------------------------------------------------------
// Function reserve()
// --------------------------------------------------------------------------

///.Function.reserve.param.object.type:Spec.Packed String

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline typename Size< String<TValue, Packed<TSpec> > >::Type
reserve(
    String<TValue, Packed<TSpec> > & seq, 
    TSize_ new_capacity,
    Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;

    typedef String<TValue, Packed<TSpec> > TString;
    typedef typename Size<TString>::Type TSize;
    TSize ret_value = reserve(host(seq), PackedConsts_<TString>::toHostLength(new_capacity), tag);
    return ret_value * PackedConsts_<TString>::VALUES_PER_WORD;
}

// ****************************************************************************
// Functions for Packed String Iter
// ****************************************************************************

// --------------------------------------------------------------------------
// Function container()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline typename Parameter_<TContainer>::Type 
container(Iter<TContainer, Packed<THostspec> > & me)
{
    SEQAN_CHECKPOINT;
    return _toParameter<TContainer>(me.data_container);
}

template <typename TContainer, typename THostspec>
inline typename Parameter_<TContainer>::Type 
container(Iter<TContainer, Packed<THostspec> > const & me)
{
    SEQAN_CHECKPOINT;
    return _toParameter<TContainer>(me.data_container);
}

// --------------------------------------------------------------------------
// Function setContainer()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec, typename TContainer2>
inline void
setContainer(Iter<TContainer, Packed<THostspec> > & me,
             TContainer2 container_)
{
    SEQAN_CHECKPOINT;
   typedef Iter<TContainer, Packed<THostspec> > TIter;
    typename Position<TIter>::Type pos = position(me);
    me.data_container = _toPointer(container_);
    setPosition(me, pos);
}

// --------------------------------------------------------------------------
// Function hostIterator()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline typename HostIterator<Iter<TContainer, Packed<THostspec> > >::Type &
hostIterator(Iter<TContainer, Packed<THostspec> > & me)
{
    SEQAN_CHECKPOINT;
    return me.data_iterator;
}

template <typename TContainer, typename THostspec>
inline typename HostIterator<Iter<TContainer, Packed<THostspec> > const>::Type  &
hostIterator(Iter<TContainer, Packed<THostspec> > const & me)
{
    SEQAN_CHECKPOINT;
    return me.data_iterator;
}

// --------------------------------------------------------------------------
// Helper Function _bitpos()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline unsigned char &
_bitpos(Iter<TContainer, Packed<THostspec> > & me)
{
    SEQAN_CHECKPOINT;
    return me.data_bitpos;
}
template <typename TContainer, typename THostspec>
inline unsigned char
_bitpos(Iter<TContainer, Packed<THostspec> > const & me)
{
    SEQAN_CHECKPOINT;
    return me.data_bitpos;
}

// --------------------------------------------------------------------------
// Function position()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline typename Position<Iter<TContainer, Packed<THostspec> > const>::Type 
position(Iter<TContainer, Packed<THostspec> > const & me)
{
    SEQAN_CHECKPOINT;
    typedef typename Host<TContainer>::Type THost;
    THost const & host_ = host(container(me));
    return (hostIterator(me) - begin(host_)) * PackedConsts_<TContainer>::VALUES_PER_WORD + _bitpos(me) / PackedConsts_<TContainer>::BITS_PER_VALUE;
}

// --------------------------------------------------------------------------
// Function setPosition()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec, typename TPosition>
inline void 
setPosition(Iter<TContainer, Packed<THostspec> > & me,
            TPosition pos_)
{
    SEQAN_CHECKPOINT;
    hostIterator(me) = begin(host(container(me))) + pos_ / PackedConsts_<TContainer>::VALUES_PER_WORD;
    _bitpos(me) = (pos_ % PackedConsts_<TContainer>::VALUES_PER_WORD) * PackedConsts_<TContainer>::BITS_PER_VALUE;
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline typename Reference<Iter<TContainer, Packed<THostspec> > >::Type 
value(Iter<TContainer, Packed<THostspec> > & me)
{
    SEQAN_CHECKPOINT;
    return typename Reference<Iter<TContainer, Packed<THostspec> > >::Type(me);
}

template <typename TContainer, typename THostspec>
inline typename Reference<Iter<TContainer, Packed<THostspec> > const>::Type 
value(Iter<TContainer, Packed<THostspec> > const & me)
{
    SEQAN_CHECKPOINT;
    return typename Reference<Iter<TContainer, Packed<THostspec> > const>::Type(me);
}

// --------------------------------------------------------------------------
// Function getValue()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline typename GetValue<Iter<TContainer, Packed<THostspec> > >::Type 
getValue(Iter<TContainer, Packed<THostspec> > & me)
{
    SEQAN_CHECKPOINT;
    return (value(hostIterator(me)) >> _bitpos(me)) & PackedConsts_<TContainer>::VALUE_MASK;
}

template <typename TContainer, typename THostspec>
inline typename GetValue<Iter<TContainer, Packed<THostspec> > const>::Type 
getValue(Iter<TContainer, Packed<THostspec> > const & me)
{
    SEQAN_CHECKPOINT;
    return (value(hostIterator(me)) >> _bitpos(me)) & PackedConsts_<TContainer>::VALUE_MASK;
}

// --------------------------------------------------------------------------
// Function assignValue()
// --------------------------------------------------------------------------

template <typename TIter, typename TValue>
inline void
_assignValuePackedStringIterator(TIter & me,
                                    TValue & _value)
{
    typedef typename Container<TIter>::Type TContainer;
    typedef typename Host<TContainer>::Type THost;
    typedef typename Value<THost>::Type THostValue;
    THostValue mask_ = PackedConsts_<TContainer>::VALUE_MASK << _bitpos(me);
    THostValue val_;
    assign(val_, _value);
    val_ <<= _bitpos(me);

    assignValue(hostIterator(me), (getValue(hostIterator(me)) & ~(mask_)) | val_);
}

template <typename TContainer, typename THostspec, typename TValue>
inline void
assignValue(Iter<TContainer, Packed<THostspec> > const & me,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    typedef Iter<TContainer, Packed<THostspec> > const TIterator;
    typename Value<TIterator>::Type _temp_value = _value; //conversion
    _assignValuePackedStringIterator(me, _temp_value);
}

template <typename TContainer, typename THostspec, typename TValue>
inline void
assignValue(Iter<TContainer, Packed<THostspec> > & me,
            TValue const & _value)
{
    SEQAN_CHECKPOINT;
    typedef Iter<TContainer, Packed<THostspec> > TIter;
    assignValue(static_cast<TIter const &>(me), _value);
}

// --------------------------------------------------------------------------
// Function moveValue()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec, typename TValue>
inline void
moveValue(Iter<TContainer, Packed<THostspec> > & me,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    assignValue(me, _value);
}

template <typename TContainer, typename THostspec, typename TValue>
inline void
moveValue(Iter<TContainer, Packed<THostspec> > const & me,
          TValue const & _value)
{
    SEQAN_CHECKPOINT;
    assignValue(me, _value);
}

// --------------------------------------------------------------------------
// Function valueConstruct()
// --------------------------------------------------------------------------

//emulate construction and destruction 

template <typename TContainer, typename THostspec>
inline void
valueConstruct(Iter<TContainer, Packed<THostspec> > const & /*it*/)
{
    // TODO(holtgrew): Why not assign default-constructed?
}

template <typename TContainer, typename THostspec, typename TParam>
inline void
valueConstruct(Iter<TContainer, Packed<THostspec> > const & it,
               TParam const & param_)
{
    assignValue(it, param_);
}

template <typename TContainer, typename THostspec, typename TParam>
inline void
valueConstruct(Iter<TContainer, Packed<THostspec> > const & it,
               TParam const & param_,
               Move const & /*tag*/)
{
    moveValue(it, param_);
}

// --------------------------------------------------------------------------
// Function valueDestruct()
// --------------------------------------------------------------------------

// Packed strings cannot contain non-POD data types.

template <typename TContainer, typename THostspec>
inline void
valueDestruct(Iter<TContainer, Packed<THostspec> > const & /*it*/)
{
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline bool 
operator==(Iter<TContainer, Packed<THostspec> > const & left,
           Iter<TContainer, Packed<THostspec> > const & right)
{
    SEQAN_CHECKPOINT;
    return (hostIterator(left) == hostIterator(right)) && (_bitpos(left) == _bitpos(right));
}

// --------------------------------------------------------------------------
// Function operator!=()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline bool 
operator!=(Iter<TContainer, Packed<THostspec> > const & left,
           Iter<TContainer, Packed<THostspec> > const & right)
{
    SEQAN_CHECKPOINT;
    return (hostIterator(left) != hostIterator(right)) || (_bitpos(left) != _bitpos(right));
}

// --------------------------------------------------------------------------
// Function operator>()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline bool 
operator>(Iter<TContainer, Packed<THostspec> > const & left,
          Iter<TContainer, Packed<THostspec> > const & right)
{
    SEQAN_CHECKPOINT;
    return (hostIterator(left) > hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (_bitpos(left) > _bitpos(right)));
}

// --------------------------------------------------------------------------
// Function operator>=()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline bool 
operator>=(Iter<TContainer, Packed<THostspec> > const & left,
           Iter<TContainer, Packed<THostspec> > const & right)
{
    SEQAN_CHECKPOINT;
    return (hostIterator(left) > hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (_bitpos(left) >= _bitpos(right)));
}

// --------------------------------------------------------------------------
// Function operator<()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline bool 
operator<(Iter<TContainer, Packed<THostspec> > const & left,
          Iter<TContainer, Packed<THostspec> > const & right)
{
    SEQAN_CHECKPOINT;
    return (hostIterator(left) < hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (_bitpos(left) < _bitpos(right)));
}

// --------------------------------------------------------------------------
// Function operator<=()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline bool 
operator <= (Iter<TContainer, Packed<THostspec> > const & left,
             Iter<TContainer, Packed<THostspec> > const & right)
{
    SEQAN_CHECKPOINT;
    return (hostIterator(left) < hostIterator(right)) || ((hostIterator(left) == hostIterator(right)) && (_bitpos(left) <= _bitpos(right)));
}

// --------------------------------------------------------------------------
// Function goNext()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline void
goNext(Iter<TContainer, Packed<THostspec> > & me)
{
    SEQAN_CHECKPOINT;
    int new_bitpos = _bitpos(me) + PackedConsts_<TContainer>::BITS_PER_VALUE;
    if (new_bitpos <= PackedConsts_<TContainer>::MAX_BIT_POS)
    {
        _bitpos(me) = (unsigned char) new_bitpos;
    }
    else
    {
        _bitpos(me) = 0;
        goNext(hostIterator(me));
    }
}

// --------------------------------------------------------------------------
// Function goPrevious()
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline void
goPrevious(Iter<TContainer, Packed<THostspec> > & me)
{
    SEQAN_CHECKPOINT;
    int new_bitpos = _bitpos(me) - PackedConsts_<TContainer>::BITS_PER_VALUE;
    if (new_bitpos >= 0)
    {
        _bitpos(me) = (unsigned char) new_bitpos;
    }
    else
    {
        _bitpos(me) = PackedConsts_<TContainer>::MAX_BIT_POS;
        goPrevious(hostIterator(me));
    }
}

// --------------------------------------------------------------------------
// Function operator+() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> >  
operator+(Iter<TContainer, Packed<THostspec> > const & left,
          TIntegral const & right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, Packed<THostspec> >(container(left), position(left) + right);
}

template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> >  
operator+(TIntegral const & left,
          Iter<TContainer, Packed<THostspec> > const & right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, Packed<THostspec> >(container(right), position(right) + left);
}

// --------------------------------------------------------------------------
// Function operator+=() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> > &
operator+=(Iter<TContainer, Packed<THostspec> > & left,
           TIntegral const & right)
{
    SEQAN_CHECKPOINT;
    setPosition(left, position(left) + right);
    return left;
}

// --------------------------------------------------------------------------
// Function operator-() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> >  
operator-(Iter<TContainer, Packed<THostspec> > const & left,
          TIntegral const & right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, Packed<THostspec> >(container(left), position(left) - right);
}

// --------------------------------------------------------------------------
// Function operator-=() for (iter, integral)
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec, typename TIntegral>
inline Iter<TContainer, Packed<THostspec> > &
operator-=(Iter<TContainer, Packed<THostspec> > & left,
           TIntegral const & right)
{
    SEQAN_CHECKPOINT;
    setPosition(left, position(left) - right);
    return left;
}

// --------------------------------------------------------------------------
// Function operator-() for (iter, iter)
// --------------------------------------------------------------------------

template <typename TContainer, typename THostspec>
inline typename Difference<Iter<TContainer, Packed<THostspec> > >::Type  
operator-(Iter<TContainer, Packed<THostspec> > const & left,
          Iter<TContainer, Packed<THostspec> > const & right)
{
    SEQAN_CHECKPOINT;
    return position(left) - position(right);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_PACKED_H_
