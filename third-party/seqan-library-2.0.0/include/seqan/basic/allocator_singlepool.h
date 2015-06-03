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

/*!
 * @class SinglePoolAllocator
 * @extends Allocator
 * @headerfile <seqan/basic.h>
 * @brief Allocator that pools memory blocks of a specific size.
 *
 * @signature template <unsigned SIZE, typename TParentAllocator>
 *            class Allocator;
 *
 * @tparam SIZE             The size of the blocks.
 * @tparam TParentAllocator The parent allocator to use.
 *
 * A pool allocator allocates several memory blocks at once.  Freed blocks are not immediately deallocated but
 * recycled in subsequential allocations.  This way, the number of calls to the heap manager is reduced, and that
 * might speed up memory management.
 *
 * The single pool allocator only pools memory blocks of size at most $SIZE$.  Blocks of other sizes are allocated and
 * deallocated using an allocator of type $ParentAllocator$. Using the single pool allocator for blocksizes larger
 * than a few KB is not advised.
 */

template <size_t SIZE, typename TParentAllocator = SimpleAllocator>
struct SinglePool;

template <size_t SIZE, typename TParentAllocator>
struct Allocator<SinglePool<SIZE, TParentAllocator> >
{
    enum
    {
        // item must be large enough to keep a pointer to the next free item
        SIZE_PER_ITEM = SIZE < sizeof(void*)? sizeof(void*) : SIZE,
        ITEMS_PER_BLOCK = (SIZE_PER_ITEM < 0x0100) ? 0x01000 / SIZE_PER_ITEM : 16,
        STORAGE_SIZE = SIZE_PER_ITEM * ITEMS_PER_BLOCK,
        STORAGE_SIZE_MIN = SIZE_PER_ITEM
    };

    char * data_recycled_blocks;
    char * data_current_begin;
    char * data_current_end;
    char * data_current_free;
    Holder<TParentAllocator, Tristate> data_parent_allocator;

    Allocator() : data_recycled_blocks(), data_current_begin(), data_current_end(), data_current_free()
    {}

    Allocator(size_t reserve_item_count) : data_recycled_blocks()
    {
        size_t storage_size = std::max(reserve_item_count * SIZE_PER_ITEM, STORAGE_SIZE_MIN);
        allocate(parentAllocator(*this), data_current_begin, storage_size);
        data_current_end = data_current_begin + storage_size;
        data_current_free = data_current_begin;
    }

    Allocator(TParentAllocator & parent_alloc)
    {
        setValue(data_parent_allocator, parent_alloc);

        data_recycled_blocks = data_current_end = data_current_free = 0;
        //dont need to initialize data_current_begin
    }

    Allocator(size_t reserve_item_count, TParentAllocator & parent_alloc)
    {
        data_recycled_blocks = 0;

        setValue(data_parent_allocator, parent_alloc);

        size_t storage_size = std::max(reserve_item_count * SIZE_PER_ITEM, STORAGE_SIZE_MIN);
        allocate(parentAllocator(*this), data_current_begin, storage_size);
        data_current_end = data_current_begin + storage_size;
        data_current_free = data_current_begin;
    }

    // Dummy copy
    Allocator(Allocator const &) :
        data_recycled_blocks(), data_current_begin(), data_current_end(),
        data_current_free()
    {
        data_recycled_blocks = data_current_end = data_current_free = 0;
    }

    inline Allocator &
    operator=(Allocator const &)
    {
        clear(*this);
        return *this;
    }

    ~Allocator()
    {
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
    return value(me.data_parent_allocator);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <size_t SIZE, typename TParentAllocator>
void
clear(Allocator<SinglePool<SIZE, TParentAllocator> > & me)
{
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
    typedef Allocator<SinglePool<SIZE, TParentAllocator> > TAllocator;
    size_t bytes_needed = count * sizeof(TValue);

    if (bytes_needed > TAllocator::SIZE_PER_ITEM)
    {//no blocking
        allocate(parentAllocator(me), data, count, tag_);
        return;
    }

    if (bytes_needed < TAllocator::SIZE_PER_ITEM)
        bytes_needed = TAllocator::SIZE_PER_ITEM;

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
    typedef Allocator<SinglePool<SIZE, TParentAllocator> > TAllocator;

    size_t bytes_needed = count * sizeof(TValue);

    if (bytes_needed > TAllocator::SIZE_PER_ITEM)
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
