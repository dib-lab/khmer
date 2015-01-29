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
// ==========================================================================
// Implementation of the Block String which allows constant-sized, e.g.
// linear expansion
// ==========================================================================

#ifndef SEQAN_HEADER_GRAPH_STACK_H
#define SEQAN_HEADER_GRAPH_STACK_H

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class BlockString Block String
 * @extends String
 * @headerfile <seqan/sequence.h>
 * @brief String optimized for push_back, top, and pop (Stack behaviour).
 * 
 * @signature template <typename TValue, unsigned SPACE = 4096>
 *            class String<TValue, Block<SIZE> >;
 * 
 * @tparam TValue The value type, that is the type of the items/characters stored in the string.  Use
 *                @link String#Value @endlink to get the value type for a given class.
 * @tparam SIZE A positive integer that specifies the number of values in each
 *              allocated block.  Size should be a power of 2, e.g., 1024.
 * 
 */

template<unsigned int SPACE = 4096>
struct Block;

/**
.Spec.Block String:
..cat:Strings
..general:Class.String
..summary:String optimized for push_back, top, and pop (Stack behaviour).
..signature:String<TValue, Block<size> >
..param.TValue:The value type, that is the type of the items/characters stored in the string.
...remarks:Use @Metafunction.Value@ to get the value type for a given class.
..param.size:A positive integer that specifies the number of values in each allocated block.
...remarks: Size should be a power of 2, e.g., 1024.
..include:seqan/sequence.h
*/
template<typename TValue, unsigned int SPACE>
class String<TValue, Block<SPACE> >
{
    typedef String<TValue, Array<SPACE> >               TBlock;
    typedef TBlock*                                     PBlock;
    typedef Allocator< SinglePool<sizeof(TBlock)> >     TAllocator;

public:
    typedef typename Iterator<TBlock, Standard>::Type   TBlockIter;
    typedef String<PBlock>                              TBlockTable;

    TBlockTable     blocks;
    TBlockIter      blockFirst, blockLast;  // current block boundaries
    TBlockIter      lastValue;              // pointer to top value
    TAllocator      alloc;

    String():
        blockFirst(TBlockIter()),
        blockLast(TBlockIter()),
        lastValue(TBlockIter()) {}

    template<typename TSource>
    String(TSource const& source):
        blockFirst(TBlockIter()),
        blockLast(TBlockIter()),
        lastValue(TBlockIter())
    {
    SEQAN_CHECKPOINT
        assign(*this, source);
    }

    String(String const & source):
        blockFirst(TBlockIter()),
        blockLast(TBlockIter()),
        lastValue(TBlockIter())
    {
    SEQAN_CHECKPOINT
        assign(*this, source);
    }

    template<typename TSource>
    String & operator =(TSource const& source)
    {
    SEQAN_CHECKPOINT
        assign(*this, source);
        return *this;
    }

    String & operator =(String const& _other)
    {
    SEQAN_CHECKPOINT
        if (this == &_other) return *this;
        assign(*this, _other);
        return *this;
    }

    ~String()
    {
        clear(*this);
    }

    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    template<typename TPos>
    inline typename Reference<String>::Type
        operator[] (TPos pos)
    {
    SEQAN_CHECKPOINT
        return value(*this, pos);
    }

    template<typename TPos>
    inline typename Reference<String const>::Type
        operator[] (TPos pos) const
    {
    SEQAN_CHECKPOINT
        return value(*this, pos);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE>
struct DefaultOverflowImplicit< String<TValue, Block<SPACE> > >
{
    typedef Generous Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

///.Metafunction.Iterator.param.T.type:Spec.Block String
///.Metafunction.Iterator.class:Spec.Block String

template<typename TValue, unsigned int SPACE>
struct Iterator<String<TValue, Block<SPACE> >, Standard>
{
    typedef Iter<String<TValue, Block<SPACE> >, PositionIterator> Type;
};

template<typename TValue, unsigned int SPACE>
struct Iterator<String<TValue, Block<SPACE> > const, Standard>
{
    typedef Iter<String<TValue, Block<SPACE> > const, PositionIterator> Type;
};

template<typename TValue, unsigned int SPACE>
struct Iterator<String<TValue, Block<SPACE> >, Rooted>
{
    typedef Iter<String<TValue, Block<SPACE> >, PositionIterator> Type;
};

template<typename TValue, unsigned int SPACE>
struct Iterator<String<TValue, Block<SPACE> > const, Rooted>
{
    typedef Iter<String<TValue, Block<SPACE> > const, PositionIterator> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE, typename TSpec>
inline typename Iterator<String<TValue, Block<SPACE> >, Tag<TSpec> const >::Type
begin(String<TValue, Block<SPACE> > & me, Tag<TSpec> const)
{
    SEQAN_CHECKPOINT;
    return Iter<String<TValue, Block<SPACE> >, PositionIterator>(me, 0);
}

template<typename TValue, unsigned int SPACE, typename TSpec>
inline typename Iterator<String<TValue, Block<SPACE> > const, Tag<TSpec> const>::Type
begin(String<TValue, Block<SPACE> > const & me, Tag<TSpec> const)
{
    SEQAN_CHECKPOINT;
    return Iter<String<TValue, Block<SPACE> > const, PositionIterator>(me, 0);
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE, typename TSpec>
inline typename Iterator<String<TValue, Block<SPACE> >, Tag<TSpec> const >::Type
end(String<TValue, Block<SPACE> > & me, Tag<TSpec> const)
{
    SEQAN_CHECKPOINT;
    return Iter<String<TValue, Block<SPACE> >, PositionIterator>(me, length(me));
}

template<typename TValue, unsigned int SPACE, typename TSpec>
inline typename Iterator<String<TValue, Block<SPACE> > const, Tag<TSpec> const>::Type
end(String<TValue, Block<SPACE> > const & me, Tag<TSpec> const)
{
    SEQAN_CHECKPOINT;
    return Iter<String<TValue, Block<SPACE> > const, PositionIterator>(me, length(me));
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE, typename TSource>
inline void
assign(
    String<TValue, Block<SPACE> >& target,
    TSource const& source)
{
    SEQAN_CHECKPOINT;
    clear(target);
    typedef typename Iterator<TSource const, Standard>::Type TIter;
    for(TIter it = begin(source, Standard()); !atEnd(it, source); goNext(it))
        push(target, *it);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE, typename TPos>
inline typename Reference<String<TValue, Block<SPACE> > >::Type
value(
    String<TValue, Block<SPACE> >& stack,
    TPos const pos)
{
    SEQAN_CHECKPOINT;
    return value(*(stack.blocks[pos / SPACE]), pos % SPACE);
}

template<typename TValue, unsigned int SPACE, typename TPos>
inline typename Reference<String<TValue, Block<SPACE> > >::Type
value(
    String<TValue, Block<SPACE> > const& stack,
    TPos const pos)
{
    SEQAN_CHECKPOINT;
    return value(*(stack.blocks[pos / SPACE]), pos % SPACE);
}

// ----------------------------------------------------------------------------
// Function atEnd()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Probably wrong place?
template<typename TValue, unsigned int SPACE, typename TIteratorSpec>
inline bool
atEnd(
    Iter<String<TValue, Block<SPACE> >, TIteratorSpec>& it,
    String<TValue, Block<SPACE> >& container)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<String<TValue, Block<SPACE> >, Standard>::Type TIter;
    TIter endIt = end(container, Standard());
    return (it == endIt);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE>
inline void
clear(String<TValue, Block<SPACE> >& me)
{
    SEQAN_CHECKPOINT;
    typedef String<TValue, Block<SPACE> >           TBlockString;
    typedef typename TBlockString::TBlockTable      TBlockTable;
    typedef typename Iterator<TBlockTable, Standard>::Type  TIter;

    TIter it = begin(me.blocks), itEnd = end(me.blocks);
    while (it != itEnd) {
        valueDestruct(it);
        deallocate(me.alloc, *it, 1);
        ++it;
    }
    clear(me.blocks);
    me.lastValue = me.blockLast = typename TBlockString::TBlockIter();
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE, typename TSize2, typename TExpand>
inline typename Size< String<TValue, Block<SPACE> > >::Type
resize(String<TValue, Block<SPACE> > & me,
    TSize2 new_length,
    Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    typedef String<TValue, Block<SPACE> >           TBlockString;
    typedef typename Size<TBlockString>::Type       TSize;
    TSize len = length(me);

    if ((TSize)new_length > len)
    {
        for (; len < (TSize)new_length; ++len) push(me);
    }
    else if ((TSize)new_length < len)
    {
        for (; len > (TSize)new_length; --len) pop(me);
    }
    return new_length;
}

template<typename TValue, unsigned int SPACE, typename TSize2>
inline typename Size< String<TValue, Block<SPACE> > >::Type
resize(String<TValue, Block<SPACE> > & me,
    TSize2 new_length,
    Limit)
{
    SEQAN_CHECKPOINT;
    typedef String<TValue, Block<SPACE> >           TBlockString;
    typedef typename Size<TBlockString>::Type       TSize;
    TSize len = length(me);

    if (new_length > capacity(me)) new_length = capacity(me);

    if (new_length > len)
    {
        TValue val;
        for (; len < new_length; ++len) push(me, val);
    }
    else if (new_length < len)
    {
        for (; len > new_length; --len) pop(me);
    }
    return new_length;
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Why is this only a dummy implementation?
///.Function.reserve.param.object.type:Spec.Block String
///.Function.reserve.class:Spec.Block String
/*
template <typename TValue, unsigned int SPACE, typename TSize, typename TExpand>
inline typename Size< String<TValue, Block<SPACE> > >::Type
reserve(
    String<TValue, Block<SPACE> >& me,
    TSize new_capacity,
    Tag<TExpand> tag)
{
SEQAN_CHECKPOINT
    reserve(me.blocks, (new_capacity + SPACE - 1) / SPACE, tag);
    return capacity(me.blocks) * SPACE;
}
*/

// dummy implementation
template<typename TValue, unsigned int SPACE, typename TSize, typename TExpand>
inline typename Size< String<TValue, Block<SPACE> > >::Type
reserve(String<TValue, Block<SPACE> > & /*me*/,
    TSize new_capacity,
    Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    return new_capacity;
}

// ----------------------------------------------------------------------------
// Function append()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE, typename TSource, typename TExpand>
inline void
append(
    String<TValue, Block<SPACE> >& me,
    TSource const& source,
    Tag<TExpand> /*tag*/)
{
    SEQAN_CHECKPOINT;
    typedef typename Iterator<TSource const, Standard>::Type TIter;
    for(TIter it = begin(source, Standard()); !atEnd(it, source); goNext(it))
        appendValue(me, *it);
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------

///.Function.appendValue.param.target.type:Spec.Block String

template<typename TValue, unsigned int SPACE, typename TVal, typename TExpand>
inline void
appendValue(
    String<TValue, Block<SPACE> >& me,
    TVal const& source,
    Tag<TExpand> tag)
{
    // TODO(holtgrew): Why does this operate on raw memory instead of using appendValue(me.blocks[last], X)?
    SEQAN_CHECKPOINT;
    if (me.lastValue == me.blockLast) {
        typename Size< String<TValue, Block<SPACE> > >::Type last = length(me.blocks);

        resize(me.blocks, last + 1, tag);
        allocate(me.alloc, me.blocks[last], 1);
        valueConstruct(me.blocks[last]);
        me.lastValue = me.blockFirst = begin(*me.blocks[last]);
        me.blockLast = (me.blockFirst + (SPACE - 1));
        back(me.blocks)->data_end += 1;
    } else {
        ++me.lastValue;
        back(me.blocks)->data_end += 1;
    }
    valueConstruct(me.lastValue, source);
}

// ----------------------------------------------------------------------------
// Function push()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE, typename TVal>
inline void
push(
    String<TValue, Block<SPACE> >& me,
    TVal const& source)
{
    appendValue(me, source);
}

template<typename TValue, unsigned int SPACE>
inline void
push(String<TValue, Block<SPACE> >& me)
{
    SEQAN_CHECKPOINT;
    if (me.lastValue == me.blockLast) {
        typename Size< String<TValue, Block<SPACE> > >::Type last = length(me.blocks);

        resize(me.blocks, last + 1, typename DefaultOverflowImplicit<String<TValue, Block<SPACE> > >::Type());
        allocate(me.alloc, me.blocks[last], 1);
        me.lastValue = me.blockFirst = begin(*me.blocks[last]);
        me.blockLast = (me.blockFirst + (SPACE - 1));
        back(me.blocks)->data_end += 1;
    } else {
        ++me.lastValue;
        back(me.blocks)->data_end += 1;
    }
    valueConstruct(me.lastValue);
}

// ----------------------------------------------------------------------------
// Function push_back()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Breaks naming-conventions.
template<typename TValue, unsigned int SPACE, typename TVal>
inline void
push_back(
    String<TValue, Block<SPACE> >& me,
    TVal const& source)
{
    appendValue(me, source);
}

// ----------------------------------------------------------------------------
// Function top()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE>
inline TValue &
top(String<TValue, Block<SPACE> > & me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT_MSG(empty(me), "top() called on an empty string.");

    return *me.lastValue;
}

template<typename TValue, unsigned int SPACE>
inline TValue const &
top(String<TValue, Block<SPACE> > const & me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT_MSG(empty(me), "top() called on an empty string.");

    return *me.lastValue;
}

template<typename TValue, unsigned int SPACE>
inline TValue &
back(String<TValue, Block<SPACE> > & me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT_MSG(empty(me), "back() called on an empty string.");

    return *me.lastValue;
}

template<typename TValue, unsigned int SPACE>
inline TValue const &
back(String<TValue, Block<SPACE> > const & me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT_MSG(empty(me), "back() called on an empty string.");

    return *me.lastValue;
}

// ----------------------------------------------------------------------------
// Function topPrev()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE>
inline TValue &
topPrev(String<TValue, Block<SPACE> > & me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ_MSG(length(me), 2u, "topPrev() called on a string with less than 2 elements.");

    if (me.lastValue != me.blockFirst)
        return *(me.lastValue - 1);
    else
        return *(begin(*me.blocks[length(me.blocks) - 1]) + (SPACE - 1));
}

template<typename TValue, unsigned int SPACE>
inline TValue const &
topPrev(String<TValue, Block<SPACE> > const& me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GEQ_MSG(length(me), 2u, "topPrev() called on a string with less than 2 elements.");

   if (me.lastValue != me.blockFirst)
        return *(me.lastValue - 1);
    else
        return *(begin(*me.blocks[length(me.blocks) - 1]) + (SPACE - 1));
}

// ----------------------------------------------------------------------------
// Function pop()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void
pop(String<TValue, TSpec> & me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT_MSG(empty(me), "pop() called on an empty string.");

    resize(me, length(me) - 1);
}

template<typename TValue, unsigned int SPACE>
inline void
pop(String<TValue, Block<SPACE> >& me)
{
    typedef typename String<TValue, Block<SPACE> >::TBlockIter TBlockIter;
    
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT_MSG(empty(me), "pop() called on an empty string.");
    
    if (me.lastValue == me.blockFirst) {
        typename Size< String<TValue, Block<SPACE> > >::Type last = length(me.blocks);

        if (last) {
            back(me.blocks)->data_end -= 1;
            valueDestruct(me.lastValue);
            valueDestruct(me.blocks[--last]);
            deallocate(me.alloc, me.blocks[last], 1);
            resize(me.blocks, last);
            if (last) {
                me.blockFirst = begin(*me.blocks[--last]);
                me.lastValue = me.blockLast = (me.blockFirst + (SPACE - 1));
            } else
            {
                me.lastValue = me.blockFirst = me.blockLast = TBlockIter();
            }
        }
    } else {
        back(me.blocks)->data_end -= 1;
        valueDestruct(me.lastValue);
        --me.lastValue;
    }
}

// ----------------------------------------------------------------------------
// Function pop_back()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Breaks naming-conventions.
template<typename TValue, unsigned int SPACE>
inline void
pop_back(String<TValue, Block<SPACE> >& me)
{
    SEQAN_CHECKPOINT;
    pop(me);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE>
inline bool
empty(String<TValue, Block<SPACE> > const& me)
{
    SEQAN_CHECKPOINT;
    return length(me.blocks) == 0;
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE>
inline typename Size<String<TValue, Block<SPACE> > >::Type
length(String<TValue, Block<SPACE> > const & me)
{
    SEQAN_CHECKPOINT;
    if (length(me.blocks))
        return (length(me.blocks) - 1) * SPACE + (me.lastValue - me.blockFirst) + 1;
    else
        return 0;
}

// ----------------------------------------------------------------------------
// Function capacity()
// ----------------------------------------------------------------------------

template<typename TValue, unsigned int SPACE>
inline typename Size<String<TValue, Block<SPACE> > >::Type
capacity(String<TValue, Block<SPACE> > const & me)
{
    SEQAN_CHECKPOINT;
    if (length(me.blocks))
        return length(me.blocks) * SPACE;
    else
        return 0;
}

} // namespace seqan

#endif //#ifndef SEQAN_HEADER_...
