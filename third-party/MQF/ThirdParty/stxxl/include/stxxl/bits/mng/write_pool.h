/***************************************************************************
 *  include/stxxl/bits/mng/write_pool.h
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

#ifndef STXXL_MNG_WRITE_POOL_HEADER
#define STXXL_MNG_WRITE_POOL_HEADER

#include <list>
#include <stxxl/bits/config.h>
#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/deprecated.h>
#include <stxxl/bits/io/request_operations.h>

#define STXXL_VERBOSE_WPOOL(msg) STXXL_VERBOSE1("write_pool[" << static_cast<void*>(this) << "]" << msg)

STXXL_BEGIN_NAMESPACE

//! \addtogroup schedlayer
//! \{

//! Implements dynamically resizable buffered writing pool.
template <class BlockType>
class write_pool : private noncopyable
{
public:
    typedef BlockType block_type;
    typedef typename block_type::bid_type bid_type;

    // a hack to make wait_any work with busy_entry type
    struct busy_entry
    {
        block_type* block;
        request_ptr req;
        bid_type bid;

        busy_entry() : block(NULL) { }
        busy_entry(const busy_entry& a) : block(a.block), req(a.req), bid(a.bid) { }
        busy_entry(block_type*& bl, request_ptr& r, bid_type& bi)
            : block(bl), req(r), bid(bi) { }

        operator request_ptr () { return req; }
    };
    typedef typename std::list<block_type*>::iterator free_blocks_iterator;
    typedef typename std::list<busy_entry>::iterator busy_blocks_iterator;

protected:
    // contains free write blocks
    std::list<block_type*> free_blocks;
    // blocks that are in writing
    std::list<busy_entry> busy_blocks;

public:
    //! Constructs pool.
    //! \param init_size initial number of blocks in the pool
    explicit write_pool(unsigned_type init_size = 1)
    {
        for (unsigned_type i = 0; i < init_size; ++i)
        {
            free_blocks.push_back(new block_type);
            STXXL_VERBOSE_WPOOL("  create block=" << free_blocks.back());
        }
    }

    void swap(write_pool& obj)
    {
        std::swap(free_blocks, obj.free_blocks);
        std::swap(busy_blocks, obj.busy_blocks);
    }

    //! Waits for completion of all ongoing write requests and frees memory.
    ~write_pool()
    {
        STXXL_VERBOSE_WPOOL("::~write_pool free_blocks.size()=" << free_blocks.size() <<
                            " busy_blocks.size()=" << busy_blocks.size());
        while (!free_blocks.empty())
        {
            STXXL_VERBOSE_WPOOL("  delete free block=" << free_blocks.back());
            delete free_blocks.back();
            free_blocks.pop_back();
        }

        try
        {
            for (busy_blocks_iterator i2 = busy_blocks.begin(); i2 != busy_blocks.end(); ++i2)
            {
                i2->req->wait();
                if (free_blocks.empty())
                    STXXL_VERBOSE_WPOOL("  delete busy block=(empty)");
                else
                    STXXL_VERBOSE_WPOOL("  delete busy block=" << free_blocks.back());
                delete i2->block;
            }
        }
        catch (...)
        { }
    }

    //! Returns number of owned blocks.
    unsigned_type size() const { return free_blocks.size() + busy_blocks.size(); }

    //! Passes a block to the pool for writing.
    //! \param block block to write. Ownership of the block goes to the pool.
    //! \c block must be allocated dynamically with using \c new .
    //! \param bid location, where to write
    //! \warning \c block must be allocated dynamically with using \c new .
    //! \return request object of the write operation
    request_ptr write(block_type*& block, bid_type bid)
    {
        STXXL_VERBOSE_WPOOL("::write: " << block << " @ " << bid);
        for (busy_blocks_iterator i2 = busy_blocks.begin(); i2 != busy_blocks.end(); ++i2)
        {
            if (i2->bid == bid) {
                assert(i2->block != block);
                STXXL_VERBOSE_WPOOL("WAW dependency");
                // try to cancel the obsolete request
                i2->req->cancel();
                // invalidate the bid of the stale write request,
                // prevents prefetch_pool from stealing a stale block
                i2->bid.storage = 0;
            }
        }
        request_ptr result = block->write(bid);
        busy_blocks.push_back(busy_entry(block, result, bid));
        block = NULL; // prevent caller from using the block any further
        return result;
    }

    //! Take out a block from the pool.
    //! \return pointer to the block. Ownership of the block goes to the caller.
    block_type * steal()
    {
        STXXL_ASSERT(size() > 0);
        if (!free_blocks.empty())
        {
            block_type* p = free_blocks.back();
            STXXL_VERBOSE_WPOOL("::steal : " << free_blocks.size() << " free blocks available, serve block=" << p);
            free_blocks.pop_back();
            return p;
        }
        STXXL_VERBOSE_WPOOL("::steal : all " << busy_blocks.size() << " are busy");
        busy_blocks_iterator completed = wait_any(busy_blocks.begin(), busy_blocks.end());
        assert(completed != busy_blocks.end()); // we got something reasonable from wait_any
        assert(completed->req->poll());         // and it is *really* completed
        block_type* p = completed->block;
        busy_blocks.erase(completed);
        check_all_busy();                       // for debug
        STXXL_VERBOSE_WPOOL("  serve block=" << p);
        return p;
    }

    // deprecated name for the steal()
    STXXL_DEPRECATED(block_type * get())
    {
        return steal();
    }

    //! Resizes size of the pool.
    //! \param new_size new size of the pool after the call
    void resize(unsigned_type new_size)
    {
        int_type diff = int_type(new_size) - int_type(size());
        if (diff > 0)
        {
            while (--diff >= 0)
            {
                free_blocks.push_back(new block_type);
                STXXL_VERBOSE_WPOOL("  create block=" << free_blocks.back());
            }

            return;
        }

        while (++diff <= 0)
            delete steal();
    }

    STXXL_DEPRECATED(request_ptr get_request(bid_type bid))
    {
        busy_blocks_iterator i2 = busy_blocks.begin();
        for ( ; i2 != busy_blocks.end(); ++i2)
        {
            if (i2->bid == bid)
                return i2->req;
        }
        return request_ptr();
    }

    bool has_request(bid_type bid)
    {
        for (busy_blocks_iterator i2 = busy_blocks.begin(); i2 != busy_blocks.end(); ++i2)
        {
            if (i2->bid == bid)
                return true;
        }
        return false;
    }

    STXXL_DEPRECATED(block_type * steal(bid_type bid))
    {
        busy_blocks_iterator i2 = busy_blocks.begin();
        for ( ; i2 != busy_blocks.end(); ++i2)
        {
            if (i2->bid == bid)
            {
                block_type* p = i2->block;
                i2->req->wait();
                busy_blocks.erase(i2);
                return p;
            }
        }
        return NULL;
    }

    // returns a block and a (potentially unfinished) I/O request associated with it
    std::pair<block_type*, request_ptr> steal_request(bid_type bid)
    {
        for (busy_blocks_iterator i2 = busy_blocks.begin(); i2 != busy_blocks.end(); ++i2)
        {
            if (i2->bid == bid)
            {
                // remove busy block from list, request has not yet been waited for!
                block_type* blk = i2->block;
                request_ptr req = i2->req;
                busy_blocks.erase(i2);

                STXXL_VERBOSE_WPOOL("::steal_request block=" << blk);
                // hand over block and (unfinished) request to caller
                return std::pair<block_type*, request_ptr>(blk, req);
            }
        }
        STXXL_VERBOSE_WPOOL("::steal_request NOT FOUND");
        // not matching request found, return a dummy
        return std::pair<block_type*, request_ptr>((block_type*)NULL, request_ptr());
    }

    void add(block_type*& block)
    {
        STXXL_VERBOSE_WPOOL("::add " << block);
        free_blocks.push_back(block);
        block = NULL; // prevent caller from using the block any further
    }

protected:
    void check_all_busy()
    {
        busy_blocks_iterator cur = busy_blocks.begin();
        int_type cnt = 0;
        while (cur != busy_blocks.end())
        {
            if (cur->req->poll())
            {
                free_blocks.push_back(cur->block);
                cur = busy_blocks.erase(cur);
                ++cnt;
                continue;
            }
            ++cur;
        }
        STXXL_VERBOSE_WPOOL("::check_all_busy : " << cnt <<
                            " are completed out of " << busy_blocks.size() + cnt << " busy blocks");
    }
};

//! \}

STXXL_END_NAMESPACE

namespace std {

template <class BlockType>
void swap(stxxl::write_pool<BlockType>& a,
          stxxl::write_pool<BlockType>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_MNG_WRITE_POOL_HEADER
// vim: et:ts=4:sw=4
