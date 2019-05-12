/***************************************************************************
 *  include/stxxl/bits/mng/block_manager.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2007 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007, 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2008-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_BLOCK_MANAGER_HEADER
#define STXXL_MNG_BLOCK_MANAGER_HEADER

#include <stxxl/bits/config.h>

#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <cstdlib>

#if STXXL_MSVC
#include <memory.h>
#endif

#include <stxxl/bits/defines.h>
#include <stxxl/bits/deprecated.h>
#include <stxxl/bits/io/request.h>
#include <stxxl/bits/io/file.h>
#include <stxxl/bits/io/create_file.h>
#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/singleton.h>
#include <stxxl/bits/mng/bid.h>
#include <stxxl/bits/mng/disk_allocator.h>
#include <stxxl/bits/mng/block_alloc.h>
#include <stxxl/bits/mng/config.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/common/simple_vector.h>

STXXL_BEGIN_NAMESPACE

#ifndef STXXL_MNG_COUNT_ALLOCATION
#define STXXL_MNG_COUNT_ALLOCATION 1
#endif // STXXL_MNG_COUNT_ALLOCATION

//! \addtogroup mnglayer
//! \{

//! Block manager class.
//!
//! Manages allocation and deallocation of blocks in multiple/single disk setting
//! \remarks is a singleton
class block_manager : public singleton<block_manager>
{
    friend class singleton<block_manager>;

    disk_allocator** disk_allocators;
    file** disk_files;

    size_t ndisks;
    block_manager();

#if STXXL_MNG_COUNT_ALLOCATION
    //! total requested allocation in bytes
    uint64 m_total_allocation;

    //! currently allocated bytes
    uint64 m_current_allocation;

    //! maximum number of bytes allocated during program run.
    uint64 m_maximum_allocation;
#endif // STXXL_MNG_COUNT_ALLOCATION

protected:
    template <class BIDType, class DiskAssignFunctor, class BIDIteratorClass>
    void new_blocks_int(
        const unsigned_type nblocks,
        const DiskAssignFunctor& functor,
        unsigned_type offset,
        BIDIteratorClass out);

public:
    //! return total number of bytes available in all disks
    uint64 get_total_bytes() const;

    //! Return total number of free disk allocations
    uint64 get_free_bytes() const;

    //! Allocates new blocks.
    //!
    //! Allocates new blocks according to the strategy
    //! given by \b functor and stores block identifiers
    //! to the range [ \b bidbegin, \b bidend)
    //! Allocation will be lined up with previous partial allocations
    //! of \b offset blocks.
    //! \param functor object of model of \b allocation_strategy concept
    //! \param bidbegin bidirectional BID iterator object
    //! \param bidend bidirectional BID iterator object
    //! \param offset advance for \b functor to line up partial allocations
    template <class DiskAssignFunctor, class BIDIteratorClass>
    void new_blocks(
        const DiskAssignFunctor& functor,
        BIDIteratorClass bidbegin,
        BIDIteratorClass bidend,
        unsigned_type offset = 0)
    {
        typedef typename std::iterator_traits<BIDIteratorClass>::value_type bid_type;
        new_blocks_int<bid_type>(std::distance(bidbegin, bidend), functor, offset, bidbegin);
    }

    //! Allocates new blocks according to the strategy
    //! given by \b functor and stores block identifiers
    //! to the output iterator \b out
    //! Allocation will be lined up with previous partial allocations
    //! of \b offset blocks.
    //! \param nblocks the number of blocks to allocate
    //! \param functor object of model of \b allocation_strategy concept
    //! \param out iterator object of OutputIterator concept
    //! \param offset advance for \b functor to line up partial allocations
    //!
    //! The \c BlockType template parameter defines the type of block to allocate
    template <class BlockType, class DiskAssignFunctor, class BIDIteratorClass>
    void new_blocks(
        const unsigned_type nblocks,
        const DiskAssignFunctor& functor,
        BIDIteratorClass out,
        unsigned_type offset = 0)
    {
        typedef typename BlockType::bid_type bid_type;
        new_blocks_int<bid_type>(nblocks, functor, offset, out);
    }

    //! Allocates a new block according to the strategy
    //! given by \b functor and stores the block identifier
    //! to bid.
    //! Allocation will be lined up with previous partial allocations
    //! of \b offset blocks.
    //! \param functor object of model of \b allocation_strategy concept
    //! \param bid BID to store the block identifier
    //! \param offset advance for \b functor to line up partial allocations
    template <typename DiskAssignFunctor, unsigned BLK_SIZE>
    void new_block(const DiskAssignFunctor& functor, BID<BLK_SIZE>& bid, unsigned_type offset = 0)
    {
        new_blocks_int<BID<BLK_SIZE> >(1, functor, offset, &bid);
    }

    //! Deallocates blocks.
    //!
    //! Deallocates blocks in the range [ \b bidbegin, \b bidend)
    //! \param bidbegin iterator object of \b bid_iterator concept
    //! \param bidend iterator object of \b bid_iterator concept
    template <class BIDIteratorClass>
    void delete_blocks(const BIDIteratorClass& bidbegin, const BIDIteratorClass& bidend);

    //! Deallocates a block.
    //! \param bid block identifier
    template <unsigned BLK_SIZE>
    void delete_block(const BID<BLK_SIZE>& bid);

    ~block_manager();

#if STXXL_MNG_COUNT_ALLOCATION
    //! return total requested allocation in bytes
    uint64 get_total_allocation() const
    { return m_total_allocation; }

    //! return currently allocated bytes
    uint64 get_current_allocation() const
    { return m_current_allocation; }

    //! return maximum number of bytes allocated during program run.
    uint64 get_maximum_allocation() const
    { return m_maximum_allocation; }
#endif // STXXL_MNG_COUNT_ALLOCATION
};

template <class BIDType, class DiskAssignFunctor, class OutputIterator>
void block_manager::new_blocks_int(
    const unsigned_type nblocks,
    const DiskAssignFunctor& functor,
    unsigned_type offset,
    OutputIterator out)
{
    typedef BIDType bid_type;
    typedef BIDArray<bid_type::t_size> bid_array_type;

    simple_vector<int_type> bl(ndisks);
    simple_vector<bid_array_type> disk_bids(ndisks);
    simple_vector<file*> disk_ptrs(nblocks);

    // choose disks by calling DiskAssignFunctor

    bl.memzero();
    for (unsigned_type i = 0; i < nblocks; ++i)
    {
        unsigned_type disk = functor(offset + i);
        disk_ptrs[i] = disk_files[disk];
        bl[disk]++;
    }

    // allocate blocks on disks

    for (unsigned_type i = 0; i < ndisks; ++i)
    {
        if (bl[i])
        {
            disk_bids[i].resize(bl[i]);
            disk_allocators[i]->new_blocks(disk_bids[i]);
        }
    }

    bl.memzero();

    OutputIterator it = out;
    for (unsigned_type i = 0; i != nblocks; ++it, ++i)
    {
        const int disk = disk_ptrs[i]->get_allocator_id();
        bid_type bid(disk_ptrs[i], disk_bids[disk][bl[disk]++].offset);
        *it = bid;
        STXXL_VERBOSE_BLOCK_LIFE_CYCLE("BLC:new    " << FMT_BID(bid));
    }

#if STXXL_MNG_COUNT_ALLOCATION
    m_total_allocation += nblocks * BIDType::size;
    m_current_allocation += nblocks * BIDType::size;
    m_maximum_allocation = STXXL_MAX(m_maximum_allocation, m_current_allocation);
#endif // STXXL_MNG_COUNT_ALLOCATION
}

template <unsigned BlockSize>
void block_manager::delete_block(const BID<BlockSize>& bid)
{
    if (!bid.valid()) {
        //STXXL_MSG("Warning: invalid block to be deleted.");
        return;
    }
    if (!bid.is_managed())
        return;  // self managed disk
    STXXL_VERBOSE_BLOCK_LIFE_CYCLE("BLC:delete " << FMT_BID(bid));
    assert(bid.storage->get_allocator_id() >= 0);
    disk_allocators[bid.storage->get_allocator_id()]->delete_block(bid);
    disk_files[bid.storage->get_allocator_id()]->discard(bid.offset, bid.size);

#if STXXL_MNG_COUNT_ALLOCATION
    m_current_allocation -= BlockSize;
#endif // STXXL_MNG_COUNT_ALLOCATION
}

template <class BIDIteratorClass>
void block_manager::delete_blocks(
    const BIDIteratorClass& bidbegin,
    const BIDIteratorClass& bidend)
{
    for (BIDIteratorClass it = bidbegin; it != bidend; it++)
    {
        delete_block(*it);
    }
}

// in bytes
#ifndef STXXL_DEFAULT_BLOCK_SIZE
    #define STXXL_DEFAULT_BLOCK_SIZE(type) (2 * 1024 * 1024) // use traits
#endif

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_MNG_BLOCK_MANAGER_HEADER
// vim: et:ts=4:sw=4
