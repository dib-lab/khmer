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
// Allocator that pools one or more consecutive memory blocks of a specific
// size.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALLOCATOR_CHUNKPOOL_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALLOCATOR_CHUNKPOOL_H_

#include <seqan/basic/allocator_interface.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Chunk Pool Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:Allocator that pools one or more consecutive memory blocks of a specific size.
..signature:Allocator< ChunkPool<SIZE, MAX_COUNT, ParentAllocator> >
..param.SIZE:Size of memory blocks that are pooled.
...value:An unsigned integer with $SIZE >= sizeof(void *)$.
..param.MAX_COUNT:Maximum number of consecutive memory blocks that are pooled.
...default:26
...remarks:Longer "chunks" are allocated and deallocated without pooling.
..param.ParentAllocator:An allocator that is by the pool allocator used to allocate memory.
...default:@Spec.Simple Allocator@
...note:The multi pool allocator only supports @Function.clear@ if this function is also implemented for $ParentAllocator$.
..remarks:A pool allocator allocates several memory blocks at once. 
Freed blocks are not immediately deallocated but recycled in subsequential allocations.
This way, the number of calls to the heap manager is reduced, and that speeds up memory management.
...text:Note that memory blocks of size different than $SIZE$, $2*SIZE$, $3*SIZE$, ..., $MAX_COUNT * SIZE$ 
are not pooled but immediately allocated and deallocated using $ParentAllocator$.
..include:seqan/basic.h
*/

template <
    size_t SIZE, 
    size_t MAX_COUNT = 26, 
    typename TParentAllocator = Allocator<SimpleAlloc<Default> > >
struct ChunkPool;

template <size_t SIZE, size_t MAX_COUNT, typename TParentAllocator>
struct Allocator<ChunkPool<SIZE, MAX_COUNT, TParentAllocator> >
{
    enum
    {
        STORAGE_SIZE_1 = 0x1000UL,
        STORAGE_SIZE_2 = SIZE * MAX_COUNT * 8,
        STORAGE_SIZE_UPPER = (STORAGE_SIZE_1 > STORAGE_SIZE_2) ? STORAGE_SIZE_1 : STORAGE_SIZE_2,
        ITEMS_PER_STORAGE = STORAGE_SIZE_UPPER / SIZE,
        STORAGE_SIZE = ITEMS_PER_STORAGE * SIZE,

        STORAGE_SIZE_MIN = SIZE * MAX_COUNT //minimal storage size
    };

    char * data_recycled_blocks [MAX_COUNT];
    char * data_current_begin;
    char * data_current_end;
    char * data_current_free;
    Holder<TParentAllocator, Tristate> data_parent_allocator;

    Allocator()
    {
        SEQAN_CHECKPOINT;
        ::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
        data_current_end = data_current_free = 0;
        //dont need to initialize data_current_begin
    }

    Allocator(size_t reserve_item_count)
    {
        SEQAN_CHECKPOINT;
        ::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));

        size_t storage_size = (reserve_item_count * SIZE > STORAGE_SIZE_MIN) ? reserve_item_count * SIZE : STORAGE_SIZE_MIN;
        allocate( parentAllocator( *this ), data_current_begin, storage_size );
        data_current_end = data_current_begin + storage_size;
        data_current_free = data_current_begin;
    }

    Allocator(TParentAllocator & parent_alloc)
    {
        SEQAN_CHECKPOINT;
        ::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
        data_current_end = data_current_free = 0;
        //dont need to initialize data_current_begin

        setValue(data_parent_allocator, parent_alloc);
    }

    Allocator(size_t reserve_item_count, TParentAllocator & parent_alloc)
    {
        SEQAN_CHECKPOINT;
        ::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));

        setValue(data_parent_allocator, parent_alloc);

        size_t storage_size = (reserve_item_count * SIZE > STORAGE_SIZE_MIN) ? reserve_item_count * SIZE : STORAGE_SIZE_MIN;
        allocate( parentAllocator( *this ), data_current_begin, storage_size );
        data_current_end = data_current_begin + storage_size;
        data_current_free = data_current_begin;
    }

    //Dummy copy
    Allocator(Allocator const &)
    {
        ::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
        data_current_end = data_current_free = 0;
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

template <size_t SIZE, size_t MAX_COUNT, typename TParentAllocator>
inline TParentAllocator &
parentAllocator(Allocator<ChunkPool<SIZE, MAX_COUNT, TParentAllocator> > & me)
{
    SEQAN_CHECKPOINT;
    return value(me.data_parent_allocator);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <size_t SIZE, size_t MAX_COUNT, typename TParentAllocator>
void
clear(Allocator<ChunkPool<SIZE, MAX_COUNT, TParentAllocator> > & me)
{
    SEQAN_CHECKPOINT;
    ::std::memset(me.data_recycled_blocks, 0, sizeof(me.data_recycled_blocks));
    me.data_current_end = me.data_current_free = 0;

    clear(parentAllocator(me));
}

// ----------------------------------------------------------------------------
// Function allocate()
// ----------------------------------------------------------------------------

template <size_t SIZE, size_t MAX_COUNT, typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void
allocate(Allocator<ChunkPool<SIZE, MAX_COUNT, TParentAllocator> > & me, 
         TValue * & data,
         TSize count,
         Tag<TUsage> const tag_)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GT(count, static_cast<TSize>(0));

    typedef Allocator<ChunkPool<SIZE, MAX_COUNT, TParentAllocator> > TAllocator;

    char * ptr;

    if ((sizeof(TValue) != SIZE) || ((size_t) count > MAX_COUNT))
    {//no blocking
        return allocate(parentAllocator(me), data, count, tag_);
    }

    size_t bytes_needed = count * SIZE;
    if (me.data_recycled_blocks[count - 1])
    {//use recycled
        ptr = me.data_recycled_blocks[count - 1];
        me.data_recycled_blocks[count - 1] = * reinterpret_cast<char **>(ptr);
    }
    else
    {//use new
        ptr = me.data_current_free;
        if (ptr + bytes_needed > me.data_current_end)
        {//not enough free space in current storage: allocate new
            size_t rest_block_number = (me.data_current_end - me.data_current_free) / SIZE;
            if (ptr && rest_block_number)
            {//link rest to recycle list
                *reinterpret_cast<char **>(ptr) = me.data_recycled_blocks[rest_block_number - 1];
                me.data_recycled_blocks[rest_block_number - 1] = reinterpret_cast<char *>(ptr);
            }

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

template <size_t SIZE, size_t MAX_COUNT, typename TParentAllocator, typename TValue, typename TSize, typename TUsage>
inline void
deallocate(Allocator<ChunkPool<SIZE, MAX_COUNT, TParentAllocator> > & me,
           TValue * data,
           TSize count,
           Tag<TUsage> const tag_)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_GT(count, 0);

    if ((sizeof(TValue) != SIZE) || (static_cast<size_t>(count) > MAX_COUNT))
    {//no blocking
        return deallocate(parentAllocator(me), data, count, tag_);
    }

    //link in recycling list
    *reinterpret_cast<char **>(data) = me.data_recycled_blocks[count - 1];
    me.data_recycled_blocks[count - 1] = reinterpret_cast<char *>(data);
}

}  // namespace seqan

#endif  // SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALLOCATOR_CHUNKPOOL_H_
