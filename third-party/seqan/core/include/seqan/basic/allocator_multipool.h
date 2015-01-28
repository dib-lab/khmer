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
// Allocator that pools blocks up to a given size.
// ==========================================================================

#ifndef SEQAN_BASIC_BASIC_ALLOCATOR_MULTIPOOL_H_
#define SEQAN_BASIC_BASIC_ALLOCATOR_MULTIPOOL_H_

#include <seqan/basic/allocator_interface.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Multi Pool Allocator:
..cat:Allocators
..general:Class.Allocator
..summary:Allocator that pools memory blocks.
..signature:Allocator MultiPool<ParentAllocator, BLOCKING_LIMIT> >
..param.ParentAllocator:An allocator that is by the pool allocator used to allocate memory.
...default:@Spec.Simple Allocator@
...note:The multi pool allocator only supports @Function.clear@ if this function is also implemented for $ParentAllocator$.
..remarks:A pool allocator allocates several memory blocks at once. 
..param.BLOCKING_LIMIT:The maximum size for memory blocks to be pooled.
...default:256
Freed blocks are not immediately deallocated but recycled in subsequential allocations.
This way, the number of calls to the heap manager is reduced, and that speeds up memory management.
...text:Note that memory blocks larger than $BLOCKING_LIMIT$ are not pooled 
but immediately allocated and deallocated using $ParentAllocator$.
..include:seqan/basic.h
*/

template <typename TParentAllocator = Allocator<SimpleAlloc<Default> >, unsigned int BLOCKING_LIMIT = 0x100>
struct MultiPool;

typedef Allocator<MultiPool<Allocator<SimpleAlloc<Default> >, 0x100> > PoolAllocator;

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT_>
struct Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT_> >
{
    enum
    {
        BLOCKING_LIMIT = BLOCKING_LIMIT_,
        GRANULARITY_BITS = 2,
        BLOCKING_COUNT = BLOCKING_LIMIT >> GRANULARITY_BITS,
        STORAGE_SIZE = 0xf80
    };

    char * data_recycled_blocks [BLOCKING_COUNT];
    char * data_current_begin [BLOCKING_COUNT];
    char * data_current_free [BLOCKING_COUNT];
    Holder<TParentAllocator, Tristate> data_parent_allocator;

    Allocator()
    {
        SEQAN_CHECKPOINT;
        // TODO(holtrew): Why not SeqAn's memset? or use using?
        ::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
        ::std::memset(data_current_begin, 0, sizeof(data_current_begin));
        ::std::memset(data_current_free, 0, sizeof(data_current_free));
    }

    Allocator(TParentAllocator & parent_alloc)
    {
        SEQAN_CHECKPOINT;
        // TODO(holtrew): Why not SeqAn's memset? or use using?
        ::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
        ::std::memset(data_current_begin, 0, sizeof(data_current_begin));
        ::std::memset(data_current_free, 0, sizeof(data_current_free));

        setValue(data_parent_allocator, parent_alloc);
    }

    //Dummy copy
    Allocator(Allocator const &)
    {
        // TODO(holtrew): Why not SeqAn's memset? or use using?
        ::std::memset(data_recycled_blocks, 0, sizeof(data_recycled_blocks));
        ::std::memset(data_current_begin, 0, sizeof(data_current_begin));
        ::std::memset(data_current_free, 0, sizeof(data_current_free));
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

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT>
inline TParentAllocator &
parentAllocator(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > & me)
{
    SEQAN_CHECKPOINT;
    return value(me.data_parent_allocator);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT>
void
clear(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > & me)
{
    SEQAN_CHECKPOINT;
    ::std::memset(me.data_recycled_blocks, 0, sizeof(me.data_recycled_blocks));
    ::std::memset(me.data_current_begin, 0, sizeof(me.data_current_begin));
    ::std::memset(me.data_current_free, 0, sizeof(me.data_current_free));

    clear(parentAllocator(me));
}

// ----------------------------------------------------------------------------
// Helper function _allocatorBlockNumber()
// ----------------------------------------------------------------------------

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT>
inline unsigned int
_allocatorBlockNumber(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > &,
                      size_t size_)
{
    SEQAN_CHECKPOINT;
    typedef Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > TAllocator;

    SEQAN_ASSERT_GT(size_, 0u);

    if (size_ < BLOCKING_LIMIT)
    {//blocks
        return size_ >> TAllocator::GRANULARITY_BITS;
    }
    else
    {//no blocking
        return TAllocator::BLOCKING_COUNT;
    }
}

// ----------------------------------------------------------------------------
// Function allocate()
// ----------------------------------------------------------------------------

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT, typename TValue, typename TSize, typename TUsage>
inline void
allocate(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > & me, 
         TValue * & data,
         TSize count,
         Tag<TUsage> const & tag_)
{
    SEQAN_CHECKPOINT;
    typedef Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > TAllocator;

    size_t bytes_needed = count * sizeof(TValue);
    char * ptr;

    unsigned int block_number =  _allocatorBlockNumber(me, bytes_needed);
    if (block_number == TAllocator::BLOCKING_COUNT)
    {//no blocking
        return allocate(parentAllocator(me), data, count, tag_);
    }

    bytes_needed = (block_number + 1) << TAllocator::GRANULARITY_BITS;

    if (me.data_recycled_blocks[block_number])
    {//use recycled
        ptr = me.data_recycled_blocks[block_number];
        me.data_recycled_blocks[block_number] = * reinterpret_cast<char **>(ptr);
    }
    else
    {//use new
        ptr = me.data_current_free[block_number];
        if (!ptr || (ptr + bytes_needed > me.data_current_begin[block_number] + TAllocator::STORAGE_SIZE))
        {//not enough free space in current storage: allocate new
            allocate(parentAllocator(me), ptr, (size_t) TAllocator::STORAGE_SIZE, tag_);
            me.data_current_begin[block_number] = ptr;
        }
        me.data_current_free[block_number] = ptr + bytes_needed;
    }

    data = reinterpret_cast<TValue *>(ptr);
}

// ----------------------------------------------------------------------------
// Function deallocate()
// ----------------------------------------------------------------------------

template <typename TParentAllocator, unsigned int BLOCKING_LIMIT, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > & me,
           TValue * data, 
           TSize count,
           Tag<TUsage> const tag_)
{
    SEQAN_CHECKPOINT;
    typedef Allocator<MultiPool<TParentAllocator, BLOCKING_LIMIT> > TAllocator;

    size_t bytes_needed = count * sizeof(TValue);

    unsigned int block_number = _allocatorBlockNumber(me, bytes_needed);
    if (block_number == TAllocator::BLOCKING_COUNT)
    {//no blocking
        return deallocate(parentAllocator(me), data, count, tag_);
    }

    bytes_needed = (block_number + 1) << TAllocator::GRANULARITY_BITS;

    //link in recycling list
    *reinterpret_cast<char **>(data) = me.data_recycled_blocks[block_number];
    me.data_recycled_blocks[block_number] = reinterpret_cast<char *>(data);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_ALLOCATOR_MULTIPOOL_H_
