/***************************************************************************
 *  include/stxxl/bits/mng/read_write_pool.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_READ_WRITE_POOL_HEADER
#define STXXL_MNG_READ_WRITE_POOL_HEADER

#include <stxxl/bits/mng/write_pool.h>
#include <stxxl/bits/mng/prefetch_pool.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup schedlayer
//! \{

//! Implements dynamically resizable buffered writing and prefetched reading pool.
template <typename BlockType>
class read_write_pool : private noncopyable
{
public:
    typedef BlockType block_type;
    typedef typename block_type::bid_type bid_type;
    typedef unsigned_type size_type;

protected:
    typedef write_pool<block_type> write_pool_type;
    typedef prefetch_pool<block_type> prefetch_pool_type;

    write_pool_type* w_pool;
    prefetch_pool_type* p_pool;
    bool delete_pools;

public:
    //! Constructs pool.
    //! \param init_size_prefetch initial number of blocks in the prefetch pool
    //! \param init_size_write initial number of blocks in the write pool
    explicit read_write_pool(size_type init_size_prefetch = 1, size_type init_size_write = 1)
        : delete_pools(true)
    {
        w_pool = new write_pool_type(init_size_write);
        p_pool = new prefetch_pool_type(init_size_prefetch);
    }

    STXXL_DEPRECATED(read_write_pool(prefetch_pool_type& p_pool, write_pool_type& w_pool))
        : w_pool(&w_pool), p_pool(&p_pool), delete_pools(false)
    { }

    void swap(read_write_pool& obj)
    {
        std::swap(w_pool, obj.w_pool);
        std::swap(p_pool, obj.p_pool);
        std::swap(delete_pools, obj.delete_pools);
    }

    //! Waits for completion of all ongoing requests and frees memory.
    ~read_write_pool()
    {
        if (delete_pools) {
            delete w_pool;
            delete p_pool;
        }
    }

    //! Returns number of blocks owned by the write_pool.
    size_type size_write() const { return w_pool->size(); }

    //! Returns number of blocks owned by the prefetch_pool.
    size_type size_prefetch() const { return p_pool->size(); }

    //! Resizes size of the pool.
    //! \param new_size new size of the pool after the call
    void resize_write(size_type new_size)
    {
        w_pool->resize(new_size);
    }

    //! Resizes size of the pool.
    //! \param new_size new size of the pool after the call
    void resize_prefetch(size_type new_size)
    {
        p_pool->resize(new_size);
    }

    // WRITE POOL METHODS

    //! Passes a block to the pool for writing.
    //! \param block block to write. Ownership of the block goes to the pool.
    //! \c block must be allocated dynamically with using \c new .
    //! \param bid location, where to write
    //! \warning \c block must be allocated dynamically with using \c new .
    //! \return request object of the write operation
    request_ptr write(block_type*& block, bid_type bid)
    {
        request_ptr result = w_pool->write(block, bid);

        // if there is a copy of this block in the prefetch pool,
        // it is now a stale copy, so invalidate it and re-hint the block
        if (p_pool->invalidate(bid))
            p_pool->hint(bid, *w_pool);

        return result;
    }

    //! Take out a block from the pool.
    //! \return pointer to the block. Ownership of the block goes to the caller.
    block_type * steal()
    {
        return w_pool->steal();
    }

    //! Add block to write pool
    void add(block_type*& block)
    {
        w_pool->add(block);
    }

    // PREFETCH POOL METHODS

    //! Gives a hint for prefetching a block.
    //! \param bid address of a block to be prefetched
    //! \return \c true if there was a free block to do prefetch and prefetching
    //! was scheduled, \c false otherwise
    //! \note If there are no free blocks available (all blocks
    //! are already in reading or read but not retrieved by user calling \c read
    //! method) calling \c hint function has no effect
    bool hint(bid_type bid)
    {
        return p_pool->hint(bid, *w_pool);
    }

    //! Cancel a hint request in case the block is no longer desired.
    bool invalidate(bid_type bid)
    {
        return p_pool->invalidate(bid);
    }

    /*!
     * Reads block. If this block is cached block is not read but passed from
     * the cache.
     *
     * \param block block object, where data to be read to. If block was cached
     * \c block 's ownership goes to the pool and block from cache is returned
     * in \c block value.
     *
     * \param bid address of the block
     * \warning \c block parameter must be allocated dynamically using \c new .
     * \return request pointer object of read operation
     */
    request_ptr read(block_type*& block, bid_type bid)
    {
        return p_pool->read(block, bid, *w_pool);
    }

    //! Returns the request pointer for a hinted block, or an invalid NULL
    //! request in case it was not requested due to lack of prefetch buffers.
    request_ptr find_hint(bid_type bid)
    {
        return p_pool->find(bid);
    }

    //! Returns true if the blocks was hinted and the request is finished.
    bool poll_hint(bid_type bid)
    {
        return p_pool->poll(bid);
    }

    //! Add block to prefetch pool
    void add_prefetch(block_type*& block)
    {
        p_pool->add(block);
    }

    //! Take out a block from the prefetch pool, one unhinted free block must
    //! be available.
    //! \return pointer to the block. Ownership of the block goes to the caller.
    block_type * steal_prefetch()
    {
        return p_pool->steal();
    }

    //! Checks if a block is in the hinted block set.
    bool in_prefetching(bid_type bid)
    {
        return p_pool->in_prefetching(bid);
    }

    //! Returns the number of free prefetching blocks.
    unsigned_type free_size_prefetch() const
    {
        return p_pool->free_size();
    }

    //! Returns the number of busy prefetching blocks.
    unsigned_type busy_size_prefetch() const
    {
        return p_pool->busy_size();
    }
};

//! \}

STXXL_END_NAMESPACE

namespace std {

template <class BlockType>
void swap(stxxl::read_write_pool<BlockType>& a,
          stxxl::read_write_pool<BlockType>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_MNG_READ_WRITE_POOL_HEADER
// vim: et:ts=4:sw=4
