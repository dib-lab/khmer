/***************************************************************************
 *  include/stxxl/bits/containers/queue.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2005 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2009, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_QUEUE_HEADER
#define STXXL_CONTAINERS_QUEUE_HEADER

#include <vector>
#include <queue>
#include <deque>

#include <stxxl/bits/deprecated.h>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/common/simple_vector.h>
#include <stxxl/bits/common/tmeta.h>
#include <stxxl/bits/mng/read_write_pool.h>
#include <stxxl/bits/mng/write_pool.h>
#include <stxxl/bits/mng/prefetch_pool.h>

STXXL_BEGIN_NAMESPACE

#ifndef STXXL_VERBOSE_QUEUE
#define STXXL_VERBOSE_QUEUE STXXL_VERBOSE2
#endif

//! \addtogroup stlcont
//! \{

//! External FIFO queue container. \n
//! <b> Introduction </b> to queue container: see \ref tutorial_queue tutorial\n
//! <b> Design and Internals </b> of queue container: see \ref design_queue.
//!
//! \tparam ValueType type of the contained objects (POD with no references to internal memory)
//! \tparam BlockSize size of the external memory block in bytes, default is \c STXXL_DEFAULT_BLOCK_SIZE(ValueType)
//! \tparam AllocStr parallel disk allocation strategy, default is \c STXXL_DEFAULT_ALLOC_STRATEGY
//! \tparam SizeType size data type, default is \c stxxl::uint64
template <class ValueType,
          unsigned BlockSize = STXXL_DEFAULT_BLOCK_SIZE(ValueType),
          class AllocStr = STXXL_DEFAULT_ALLOC_STRATEGY,
          class SizeType = stxxl::uint64>
class queue : private noncopyable
{
public:
    typedef ValueType value_type;
    typedef AllocStr alloc_strategy_type;
    typedef SizeType size_type;
    enum {
        block_size = BlockSize
    };

    typedef typed_block<block_size, value_type> block_type;
    typedef BID<block_size> bid_type;

private:
    typedef read_write_pool<block_type> pool_type;

    size_type m_size;
    bool delete_pool;
    pool_type* pool;
    block_type* front_block;
    block_type* back_block;
    value_type* front_element;
    value_type* back_element;
    alloc_strategy_type alloc_strategy;
    unsigned_type alloc_count;
    std::deque<bid_type> bids;
    block_manager* bm;
    unsigned_type blocks2prefetch;

public:
    //! \name Constructors/Destructors
    //! \{

    //! Constructs empty queue with own write and prefetch block pool.
    //!
    //! \param D  number of parallel disks, defaulting to the configured number of scratch disks,
    //!           memory consumption will be 2 * D + 2 blocks
    //!           (first and last block, D blocks as write cache, D block for prefetching)
    explicit queue(int_type D = -1)
        : m_size(0),
          delete_pool(true),
          alloc_count(0),
          bm(block_manager::get_instance())
    {
        if (D < 1)
            D = config::get_instance()->disks_number();
        STXXL_VERBOSE_QUEUE("queue[" << this << "]::queue(D)");
        pool = new pool_type(D, D + 2);
        init();
    }

    //! Constructs empty queue with own write and prefetch block pool.
    //!
    //! \param w_pool_size  number of blocks in the write pool, must be at least 2, recommended at least 3
    //! \param p_pool_size  number of blocks in the prefetch pool, recommended at least 1
    //! \param blocks2prefetch_  defines the number of blocks to prefetch (\c front side),
    //!                          default is number of block in the prefetch pool
    explicit queue(unsigned_type w_pool_size, unsigned_type p_pool_size, int blocks2prefetch_ = -1)
        : m_size(0),
          delete_pool(true),
          alloc_count(0),
          bm(block_manager::get_instance())
    {
        STXXL_VERBOSE_QUEUE("queue[" << this << "]::queue(sizes)");
        pool = new pool_type(p_pool_size, w_pool_size);
        init(blocks2prefetch_);
    }

    //! Constructs empty queue.
    //!
    //! \param w_pool write pool
    //! \param p_pool prefetch pool
    //! \param blocks2prefetch_  defines the number of blocks to prefetch (\c front side),
    //!                          default is number of blocks in the prefetch pool
    //!  \warning Number of blocks in the write pool must be at least 2, recommended at least 3
    //!  \warning Number of blocks in the prefetch pool recommended at least 1
    STXXL_DEPRECATED(
        queue(write_pool<block_type>& w_pool, prefetch_pool<block_type>& p_pool, int blocks2prefetch_ = -1))
        : m_size(0),
          delete_pool(true),
          alloc_count(0),
          bm(block_manager::get_instance())
    {
        STXXL_VERBOSE_QUEUE("queue[" << this << "]::queue(pools)");
        pool = new pool_type(p_pool, w_pool);
        init(blocks2prefetch_);
    }

    //! Constructs empty queue.
    //!
    //! \param pool_ block write/prefetch pool
    //! \param blocks2prefetch_  defines the number of blocks to prefetch (\c front side),
    //!                          default is number of blocks in the prefetch pool
    //!  \warning Number of blocks in the write pool must be at least 2, recommended at least 3
    //!  \warning Number of blocks in the prefetch pool recommended at least 1
    queue(pool_type& pool_, int blocks2prefetch_ = -1)
        : m_size(0),
          delete_pool(false),
          pool(&pool_),
          alloc_count(0),
          bm(block_manager::get_instance())
    {
        STXXL_VERBOSE_QUEUE("queue[" << this << "]::queue(pool)");
        init(blocks2prefetch_);
    }

    //! \}

    //! \name Modifiers
    //! \{

    void swap(queue& obj)
    {
        std::swap(m_size, obj.m_size);
        std::swap(delete_pool, obj.delete_pool);
        std::swap(pool, obj.pool);
        std::swap(front_block, obj.front_block);
        std::swap(back_block, obj.back_block);
        std::swap(front_element, obj.front_element);
        std::swap(back_element, obj.back_element);
        std::swap(alloc_strategy, obj.alloc_strategy);
        std::swap(alloc_count, obj.alloc_count);
        std::swap(bids, obj.bids);
        std::swap(bm, obj.bm);
        std::swap(blocks2prefetch, obj.blocks2prefetch);
    }

    //! \}

private:
    void init(int blocks2prefetch_ = -1)
    {
        if (pool->size_write() < 2) {
            STXXL_ERRMSG("queue: invalid configuration, not enough blocks (" << pool->size_write() <<
                         ") in write pool, at least 2 are needed, resizing to 3");
            pool->resize_write(3);
        }

        if (pool->size_write() < 3) {
            STXXL_MSG("queue: inefficient configuration, no blocks for buffered writing available");
        }

        if (pool->size_prefetch() < 1) {
            STXXL_MSG("queue: inefficient configuration, no blocks for prefetching available");
        }

        front_block = back_block = pool->steal();
        back_element = back_block->begin() - 1;
        front_element = back_block->begin();
        set_prefetch_aggr(blocks2prefetch_);
    }

public:
    //! \name Miscellaneous
    //! \{

    //! Defines the number of blocks to prefetch (\c front side).
    //! This method should be called whenever the prefetch pool is resized
    //! \param blocks2prefetch_  defines the number of blocks to prefetch (\c front side),
    //!                          a negative value means to use the number of blocks in the prefetch pool
    void set_prefetch_aggr(int_type blocks2prefetch_)
    {
        if (blocks2prefetch_ < 0)
            blocks2prefetch = pool->size_prefetch();
        else
            blocks2prefetch = blocks2prefetch_;
    }

    //! Returns the number of blocks prefetched from the \c front side.
    unsigned_type get_prefetch_aggr() const
    {
        return blocks2prefetch;
    }
    //! \}

    //! \name Modifiers
    //! \{

    //! Adds an element in the queue.
    void push(const value_type& val)
    {
        if (UNLIKELY(back_element == back_block->begin() + (block_type::size - 1)))
        {
            // back block is filled
            if (front_block == back_block)
            {             // can not write the back block because it
                // is the same as the front block, must keep it memory
                STXXL_VERBOSE1("queue::push Case 1");
            }
            else if (size() < 2 * block_type::size)
            {
                STXXL_VERBOSE1("queue::push Case 1.5");
                // only two blocks with a gap in the beginning, move elements within memory
                assert(bids.empty());
                size_t gap = front_element - front_block->begin();
                assert(gap > 0);
                std::copy(front_element, front_block->end(), front_block->begin());
                std::copy(back_block->begin(), back_block->begin() + gap, front_block->begin() + (block_type::size - gap));
                std::copy(back_block->begin() + gap, back_block->end(), back_block->begin());
                front_element -= gap;
                back_element -= gap;

                ++back_element;
                *back_element = val;
                ++m_size;
                return;
            }
            else
            {
                STXXL_VERBOSE1("queue::push Case 2");
                // write the back block
                // need to allocate new block
                bid_type newbid;

                bm->new_block(alloc_strategy, newbid, alloc_count++);

                STXXL_VERBOSE_QUEUE("queue[" << this << "]: push block " << back_block << " @ " << FMT_BID(newbid));
                bids.push_back(newbid);
                pool->write(back_block, newbid);
                if (bids.size() <= blocks2prefetch) {
                    STXXL_VERBOSE1("queue::push Case Hints");
                    pool->hint(newbid);
                }
            }
            back_block = pool->steal();

            back_element = back_block->begin();
            *back_element = val;
            ++m_size;
            return;
        }
        ++back_element;
        *back_element = val;
        ++m_size;
    }

    //! Removes element from the queue.
    void pop()
    {
        assert(!empty());

        if (UNLIKELY(front_element == front_block->begin() + (block_type::size - 1)))
        {
            // if there is only one block, it implies ...
            if (back_block == front_block)
            {
                STXXL_VERBOSE1("queue::pop Case 3");
                assert(size() == 1);
                assert(back_element == front_element);
                assert(bids.empty());
                // reset everything
                back_element = back_block->begin() - 1;
                front_element = back_block->begin();
                m_size = 0;
                return;
            }

            --m_size;
            if (m_size <= block_type::size)
            {
                STXXL_VERBOSE1("queue::pop Case 4");
                assert(bids.empty());
                // the back_block is the next block
                pool->add(front_block);
                front_block = back_block;
                front_element = back_block->begin();
                return;
            }
            STXXL_VERBOSE1("queue::pop Case 5");

            assert(!bids.empty());
            request_ptr req = pool->read(front_block, bids.front());
            STXXL_VERBOSE_QUEUE("queue[" << this << "]: pop block  " << front_block << " @ " << FMT_BID(bids.front()));

            // give prefetching hints
            for (unsigned_type i = 0; i < blocks2prefetch && i < bids.size() - 1; ++i)
            {
                STXXL_VERBOSE1("queue::pop Case Hints");
                pool->hint(bids[i + 1]);
            }

            front_element = front_block->begin();
            req->wait();

            bm->delete_block(bids.front());
            bids.pop_front();
            return;
        }

        ++front_element;
        --m_size;
    }
    //! \}

    //! \name Operators
    //! \{

    //! Returns a mutable reference at the back of the queue.
    value_type & back()
    {
        assert(!empty());
        return *back_element;
    }

    //! Returns a const reference at the back of the queue.
    const value_type & back() const
    {
        assert(!empty());
        return *back_element;
    }

    //! Returns a mutable reference at the front of the queue.
    value_type & front()
    {
        assert(!empty());
        return *front_element;
    }

    //! Returns a const reference at the front of the queue.
    const value_type & front() const
    {
        assert(!empty());
        return *front_element;
    }

    //! \}

    //! \name Constructors/Destructors
    //! \{

    ~queue()
    {
        if (front_block != back_block)
            pool->add(back_block);
        pool->add(front_block);

        if (delete_pool)
        {
            delete pool;
        }

        if (!bids.empty())
            bm->delete_blocks(bids.begin(), bids.end());
    }

    //! \}

    //! \name Capacity
    //! \{

    //! Returns the size of the queue.
    size_type size() const
    {
        return m_size;
    }

    //! Returns \c true if queue is empty.
    bool empty() const
    {
        return (m_size == 0);
    }

    //! \}
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_QUEUE_HEADER
// vim: et:ts=4:sw=4
