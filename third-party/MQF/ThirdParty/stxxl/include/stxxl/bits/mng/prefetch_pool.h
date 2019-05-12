/***************************************************************************
 *  include/stxxl/bits/mng/prefetch_pool.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003-2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_PREFETCH_POOL_HEADER
#define STXXL_MNG_PREFETCH_POOL_HEADER

#include <list>
#include <stxxl/bits/config.h>
#include <stxxl/bits/mng/write_pool.h>
#include <stxxl/bits/compat/hash_map.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup schedlayer
//! \{

//! Implements dynamically resizable prefetching pool.
template <class BlockType>
class prefetch_pool : private noncopyable
{
public:
    typedef BlockType block_type;
    typedef typename block_type::bid_type bid_type;

protected:
    struct bid_hash
    {
        size_t operator () (const bid_type& bid) const
        {
            size_t result = size_t(bid.storage) +
                            size_t(bid.offset & 0xffffffff) +
                            size_t(bid.offset >> 32);
            return result;
        }
#if STXXL_MSVC
        bool operator () (const bid_type& a, const bid_type& b) const
        {
            return (a.storage < b.storage) || (a.storage == b.storage && a.offset < b.offset);
        }
        enum
        {                               // parameters for hash table
            bucket_size = 4,            // 0 < bucket_size
            min_buckets = 8             // min_buckets = 2 ^^ N, 0 < N
        };
#endif
    };
    typedef std::pair<block_type*, request_ptr> busy_entry;
    typedef typename compat_hash_map<bid_type, busy_entry, bid_hash>::result hash_map_type;
    typedef typename std::list<block_type*>::iterator free_blocks_iterator;
    typedef typename hash_map_type::iterator busy_blocks_iterator;

    //! contains free prefetch blocks
    std::list<block_type*> free_blocks;

    //! blocks that are in reading or already read but not retrieved by user
    hash_map_type busy_blocks;

    //! count number of free blocks, since traversing the std::list is slow.
    unsigned_type free_blocks_size;

public:
    //! Constructs pool.
    //! \param init_size initial number of blocks in the pool
    explicit prefetch_pool(unsigned_type init_size = 1)
        : free_blocks_size(init_size)
    {
        unsigned_type i = 0;
        for ( ; i < init_size; ++i)
            free_blocks.push_back(new block_type);
    }

    void swap(prefetch_pool& obj)
    {
        std::swap(free_blocks, obj.free_blocks);
        std::swap(busy_blocks, obj.busy_blocks);
        std::swap(free_blocks_size, obj.free_blocks_size);
    }

    //! Waits for completion of all ongoing read requests and frees memory.
    virtual ~prefetch_pool()
    {
        while (!free_blocks.empty())
        {
            delete free_blocks.back();
            free_blocks.pop_back();
        }

        try
        {
            busy_blocks_iterator i2 = busy_blocks.begin();
            for ( ; i2 != busy_blocks.end(); ++i2)
            {
                i2->second.second->wait();
                delete i2->second.first;
            }
        }
        catch (...)
        { }
    }

    //! Returns number of owned blocks.
    unsigned_type size() const
    {
        return free_blocks_size + busy_blocks.size();
    }

    //! Returns the number of free prefetching blocks.
    unsigned_type free_size() const
    {
        return free_blocks_size;
    }

    //! Returns the number of busy prefetching blocks.
    unsigned_type busy_size() const
    {
        return busy_blocks.size();
    }

    //! Add a new block to prefetch pool, enlarges size of pool.
    void add(block_type*& block)
    {
        free_blocks.push_back(block);
        ++free_blocks_size;
        block = NULL; // prevent caller from using the block any further
    }

    //! Take out a block from the pool, one unhinted free block must be
    //! available.
    //! \return pointer to the block. Ownership of the block goes to the caller.
    block_type * steal()
    {
        STXXL_CHECK(!free_blocks.empty());

        block_type* p = free_blocks.back();
        free_blocks.pop_back();
        --free_blocks_size;
        return p;
    }

    /*!
     * Gives a hint for prefetching a block, the block may or may not be read
     * into a prefetch buffer.
     *
     * \param bid address of a block to be prefetched
     * \return \c true if there was a free block to do prefetch and
     * prefetching was scheduled, \c false otherwise
     *
     * \note If there are no free blocks available (all blocks are already in
     * reading or read but not retrieved by user calling \c read method)
     * calling \c hint function has no effect
     */
    bool hint(bid_type bid)
    {
        // if block is already hinted, no need to hint it again
        if (in_prefetching(bid)) {
            STXXL_VERBOSE2("prefetch_pool::hint2 bid=" << bid << " was already cached");
            return true;
        }

        if (free_blocks_size) //  only if we have a free block
        {
            --free_blocks_size;
            block_type* block = free_blocks.back();
            free_blocks.pop_back();
            STXXL_VERBOSE2("prefetch_pool::hint bid=" << bid << " => prefetching");
            request_ptr req = block->read(bid);
            busy_blocks[bid] = busy_entry(block, req);
            return true;
        }
        STXXL_VERBOSE2("prefetch_pool::hint bid=" << bid << " => no free blocks for prefetching");
        return false;
    }

    /*!
     * Gives a hint for prefetching a block, the block may or may not be read
     * into a prefetch buffer. This variant checks if the write pool is
     * currently writing said block.
     *
     * \param bid address of a block to be prefetched
     * \return \c true if there was a free block to do prefetch and
     * prefetching was scheduled, \c false otherwise
     *
     * \note If there are no free blocks available (all blocks are already in
     * reading or read but not retrieved by user calling \c read method)
     * calling \c hint function has no effect
     */
    bool hint(bid_type bid, write_pool<block_type>& w_pool)
    {
        // if block is already hinted, no need to hint it again
        if (in_prefetching(bid)) {
            STXXL_VERBOSE2("prefetch_pool::hint2 bid=" << bid << " was already cached");
            return true;
        }

        if (free_blocks_size) //  only if we have a free block
        {
            --free_blocks_size;
            block_type* block = free_blocks.back();
            free_blocks.pop_back();
            if (w_pool.has_request(bid))
            {
                busy_entry wp_request = w_pool.steal_request(bid);
                STXXL_VERBOSE1("prefetch_pool::hint2 bid=" << bid << " was in write cache at " << wp_request.first);
                assert(wp_request.first != 0);
                w_pool.add(block);  //in exchange
                busy_blocks[bid] = wp_request;
                return true;
            }
            STXXL_VERBOSE2("prefetch_pool::hint2 bid=" << bid << " => prefetching");
            request_ptr req = block->read(bid);
            busy_blocks[bid] = busy_entry(block, req);
            return true;
        }
        STXXL_VERBOSE2("prefetch_pool::hint2 bid=" << bid << " => no free blocks for prefetching");
        return false;
    }

    //! Cancel a hint request in case the block is no longer desired.
    bool invalidate(bid_type bid)
    {
        busy_blocks_iterator cache_el = busy_blocks.find(bid);
        if (cache_el == busy_blocks.end())
            return false;

        // cancel request if it is a read request, there might be
        // write requests 'stolen' from a write_pool that may not be canceled
        if (cache_el->second.second->get_type() == request::READ)
            cache_el->second.second->cancel();
        // finish the request
        cache_el->second.second->wait();
        ++free_blocks_size;
        free_blocks.push_back(cache_el->second.first);
        busy_blocks.erase(cache_el);
        return true;
    }

    //! Checks if a block is in the hinted block set.
    bool in_prefetching(bid_type bid)
    {
        return (busy_blocks.find(bid) != busy_blocks.end());
    }

    //! Returns the request pointer for a hinted block, or an invalid NULL
    //! request in case it was not requested due to lack of prefetch buffers.
    request_ptr find(bid_type bid)
    {
        busy_blocks_iterator cache_el = busy_blocks.find(bid);

        if (cache_el == busy_blocks.end())
            return request_ptr(); // invalid pointer
        else
            return cache_el->second.second;
    }

    //! Returns true if the blocks was hinted and the request is finished.
    bool poll(bid_type bid)
    {
        request_ptr req = find(bid);
        return req.valid() ? req->poll() : false;
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
     *
     * \warning \c block parameter must be allocated dynamically using \c new .
     *
     * \return request pointer object of read operation
     */
    request_ptr read(block_type*& block, bid_type bid)
    {
        busy_blocks_iterator cache_el = busy_blocks.find(bid);
        if (cache_el == busy_blocks.end())
        {
            // not cached
            STXXL_VERBOSE1("prefetch_pool::read bid=" << bid << " => no copy in cache, retrieving to " << block);
            return block->read(bid);
        }

        // cached
        STXXL_VERBOSE1("prefetch_pool::read bid=" << bid << " => copy in cache exists");
        ++free_blocks_size;
        free_blocks.push_back(block);
        block = cache_el->second.first;
        request_ptr result = cache_el->second.second;
        busy_blocks.erase(cache_el);
        return result;
    }

    request_ptr read(block_type*& block, bid_type bid, write_pool<block_type>& w_pool)
    {
        // try cache
        busy_blocks_iterator cache_el = busy_blocks.find(bid);
        if (cache_el != busy_blocks.end())
        {
            // cached
            STXXL_VERBOSE1("prefetch_pool::read bid=" << bid << " => copy in cache exists");
            ++free_blocks_size;
            free_blocks.push_back(block);
            block = cache_el->second.first;
            request_ptr result = cache_el->second.second;
            busy_blocks.erase(cache_el);
            return result;
        }

        // try w_pool cache
        if (w_pool.has_request(bid))
        {
            busy_entry wp_request = w_pool.steal_request(bid);
            STXXL_VERBOSE1("prefetch_pool::read bid=" << bid << " was in write cache at " << wp_request.first);
            assert(wp_request.first != 0);
            w_pool.add(block);  //in exchange
            block = wp_request.first;
            return wp_request.second;
        }

        // not cached
        STXXL_VERBOSE1("prefetch_pool::read bid=" << bid << " => no copy in cache, retrieving to " << block);
        return block->read(bid);
    }

    //! Resizes size of the pool.
    //! \param new_size desired size of the pool. If some
    //! blocks are used for prefetching, these blocks can't be freed.
    //! Only free blocks (not in prefetching) can be freed by reducing
    //! the size of the pool calling this method.
    //! \return new size of the pool
    unsigned_type resize(unsigned_type new_size)
    {
        int_type diff = int_type(new_size) - int_type(size());
        if (diff > 0)
        {
            free_blocks_size += diff;
            while (--diff >= 0)
                free_blocks.push_back(new block_type);

            return size();
        }

        while (diff < 0 && free_blocks_size > 0)
        {
            ++diff;
            --free_blocks_size;
            delete free_blocks.back();
            free_blocks.pop_back();
        }
        return size();
    }
};

//! \}

STXXL_END_NAMESPACE

namespace std {

template <class BlockType>
void swap(stxxl::prefetch_pool<BlockType>& a,
          stxxl::prefetch_pool<BlockType>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_MNG_PREFETCH_POOL_HEADER
// vim: et:ts=4:sw=4
