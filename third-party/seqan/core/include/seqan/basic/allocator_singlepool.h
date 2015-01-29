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
// Allocator that pools blocks of a given size;  Different-sized blocks are
// not pooled.
// ==========================================================================

#ifndef SEQAN_BASIC_BASIC_ALLOCATOR_SINGLE_POOL_H_
#define SEQAN_BASIC_BASIC_ALLOCATOR_SINGLE_POOL_H_

#include <seqan/basic/allocator_interface.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Single Pool Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:Allocator that pools memory blocks of specific size.
..signature:Allocator< SinglePool<SIZE, ParentAllocator> >
..param.SIZE:Size of memory blocks that are pooled.
...value:An unsigned integer with $SIZE >= sizeof(void *)$.
..param.ParentAllocator:An allocator that is by the pool allocator used to allocate memory.
...default:@Spec.Simple Allocator@
...note:The single pool allocator only supports @Function.clear@ if this function is also implemented for $ParentAllocator$.
..remarks:A pool allocator allocates several memory blocks at once. 
Freed blocks are not immediately deallocated but recycled in subsequential allocations.
This way, the number of calls to the heap manager is reduced, and that speeds up memory management.
...text:The single pool allocator only pools memory blocks of size $SIZE$.
Blocks of other sizes are allocated and deallocated using an allocator of type $ParentAllocator$.
...text:Using the single pool allocator for blocksizes larger than some KB is not advised.
..include:seqan/basic.h
*/

template <size_t SIZE, typename TParentAllocator = SimpleAllocator>
struct SinglePool;

template <size_t SIZE, typename TParentAllocator>
struct Allocator<SinglePool<SIZE, TParentAllocator> >
{
    enum
    {
        SIZE_PER_ITEM = SIZE,
        ITEMS_PER_BLOCK = (SIZE_PER_ITEM < 0x0100) ? 0x01000 / SIZE_PER_ITEM : 16,
        STORAGE_SIZE = SIZE * ITEMS_PER_BLOCK,

        STORAGE_SIZE_MIN = SIZE
    };

    char * data_recycled_blocks;
    char * data_current_begin;
    char * data_current_end;
    char * data_current_free;
    Holder<TParentAllocator, Tristate> data_parent_allocator;

    Allocator()
    {
        SEQAN_CHECKPOINT;
        data_recycled_blocks = data_current_end = data_current_free = 0;
        //dont need to initialize data_current_begin
    }

    Allocator(size_t reserve_item_count)
    {
        SEQAN_CHECKPOINT;
        data_recycled_blocks = 0;

        size_t storage_size = (reserve_item_count * SIZE > STORAGE_SIZE_MIN) ? reserve_item_count * SIZE : STORAGE_SIZE_MIN;
        allocate( parentAllocator( *this ), data_current_begin, storage_size );
        data_current_end = data_current_begin + storage_size;
        data_current_free = data_current_begin;
    }

    Allocator(TParentAllocator & parent_alloc)
    {
        SEQAN_CHECKPOINT;
        setValue(data_parent_allocator, parent_alloc);

        data_recycled_blocks = data_current_end = data_current_free = 0;
        //dont need to initialize data_current_begin
    }

    Allocator(size_t reserve_item_count, TParentAllocator & parent_alloc)
    {
        SEQAN_CHECKPOINT;
        data_recycled_blocks = 0;

        setValue(data_parent_allocator, parent_alloc);

        size_t storage_size = (reserve_item_count * SIZE > STORAGE_SIZE_MIN) ? reserve_item_count * SIZE : STORAGE_SIZE_MIN;
        allocate( parentAllocator( *this ), data_current_begin, storage_size );
        data_current_end = data_current_begin + storage_size;
        data_current_free = data_current_begin;
    }

    //Dummy copy
    Allocator(Allocator const &)
    {
        data_recycled_blocks = data_current_end = data_current_free = 0;
        //dont need to initialize data_current_begin
    }

    inline Allocator &
    operator=(Allocator const &)
    {
        clear(*this);
        return *this;
    }

    ~Allocator()
    {
        SEQAN_CHECKPOINT;
        clear(*this);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function parentAllocator()
// ----------------------------------------------------------------------------

template <size_t SIZE, typename TParentAllocator>
inline TParentAllocator &
parentAllocator(Allocator<SinglePool<SIZE, TParentAllocator> > & me)
{
    SEQAN_CHECKPOINT;
    return value(me.data_parent_allocator);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <size_t SIZE, typename TParentAllocator>
void
clear(Allocator<SinglePool<SIZE, TParentAllocator> > & me)
{
    SEQAN_CHECKPOINT;

    me.data_recycled_blocks = me.data_current_end = me.data_current_free = 0;

    clear(parentAllocator(me));
}

// ----------------------------------------------------------------------------
// Function allocate()
// ----------------------------------------------------------------------------

template <size_t SIZE, typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void
allocate(Allocator<SinglePool<SIZE, TParentAllocator> > & me, 
         TValue * & data,
         TSize count,
         Tag<TUsage> const tag_)
{
    SEQAN_CHECKPOINT;
    typedef Allocator<SinglePool<SIZE, TParentAllocator> > TAllocator;
    size_t bytes_needed = count * sizeof(TValue);

    if (bytes_needed != TAllocator::SIZE_PER_ITEM)
    {//no blocking
        allocate(parentAllocator(me), data, count, tag_);
        return;
    }

    char * ptr;
    if (me.data_recycled_blocks)
    {//use recycled
        ptr = me.data_recycled_blocks;
        me.data_recycled_blocks = * reinterpret_cast<char **>(ptr);
    }
    else
    {//use new
        ptr = me.data_current_free;
        if (ptr + bytes_needed > me.data_current_end)
        {//not enough free space in current storage: allocate new
            allocate(parentAllocator(me), ptr, (size_t) TAllocator::STORAGE_SIZE, tag_);
            me.data_current_begin = ptr;
            me.data_current_end = ptr + TAllocator::STORAGE_SIZE;
        }
        me.data_current_free = ptr + bytes_needed;
    }

    data = reinterpret_cast<TValue *>(ptr);
}

// ----------------------------------------------------------------------------
// Function deallocate()
// ----------------------------------------------------------------------------

template <size_t SIZE, typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(Allocator<SinglePool<SIZE, TParentAllocator> > & me,
           TValue * data, 
           TSize count,
           Tag<TUsage> const tag_)
{
    SEQAN_CHECKPOINT;
    typedef Allocator<SinglePool<SIZE, TParentAllocator> > TAllocator;

    size_t bytes_needed = count * sizeof(TValue);

    if (bytes_needed != TAllocator::SIZE_PER_ITEM)
    {//no blocking
        deallocate(parentAllocator(me), data, count, tag_);
        return;
    }

    //link in recycling list
    *reinterpret_cast<char **>(data) = me.data_recycled_blocks;
    me.data_recycled_blocks = reinterpret_cast<char *>(data);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_ALLOCATOR_SINGLE_POOL_H_
