/***************************************************************************
 *  include/stxxl/bits/containers/sequence.h
 *
 *  based on include/stxxl/bits/containers/queue.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_SEQUENCE_HEADER
#define STXXL_CONTAINERS_SEQUENCE_HEADER

#include <deque>

#include <stxxl/bits/deprecated.h>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/common/tmeta.h>
#include <stxxl/bits/mng/read_write_pool.h>
#include <stxxl/bits/mng/write_pool.h>
#include <stxxl/bits/mng/prefetch_pool.h>

STXXL_BEGIN_NAMESPACE

#ifndef STXXL_VERBOSE_SEQUENCE
#define STXXL_VERBOSE_SEQUENCE STXXL_VERBOSE2
#endif

//! \addtogroup stlcont
//! \{

//! External sequence or deque container without random access. \n
//! <b> Introduction </b> to sequence container: see \ref tutorial_sequence tutorial. \n
//! <b> Design and Internals </b> of sequence container: see \ref design_queue

/**
 * Sequence is a primitive container consisting of only a sequence of blocks in
 * external memory. The sequence provides appending methods similar to a deque:
 * push_back and push_front; and also the corresponding pop functions. However,
 * different from stxxl::deque (which is a vector in disguise), the sequence
 * does not allow random access. Instead, the sequence can only be iterated
 * using streams: either from front to back or in reverse.
 *
 * As with queue and stack, sequences of pushes and pops are made efficient
 * using overlapping or read-ahead via block pools. The stream access likewise
 * uses overlapped I/O, just like stream::vector_iterator2stream.
 *
 * \tparam ValueType type of the contained objects (POD with no references to internal memory)
 * \tparam BlockSize size of the external memory block in bytes, default is \c STXXL_DEFAULT_BLOCK_SIZE(ValTp)
 * \tparam AllocStr parallel disk allocation strategy, default is \c STXXL_DEFAULT_ALLOC_STRATEGY
 * \tparam SizeType size data type, default is \c stxxl::uint64
 */
template <class ValueType,
          unsigned BlockSize = STXXL_DEFAULT_BLOCK_SIZE(ValueType),
          class AllocStr = STXXL_DEFAULT_ALLOC_STRATEGY,
          class SizeType = stxxl::uint64>
class sequence : private noncopyable
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

    typedef std::deque<bid_type> bid_deque_type;

private:
    typedef read_write_pool<block_type> pool_type;

    /// current number of items in the sequence
    size_type m_size;

    /// whether the m_pool object is own and should be deleted.
    bool m_owns_pool;

    /// read_write_pool of blocks
    pool_type* m_pool;

    /// current front block of sequence
    block_type* m_front_block;

    /// current back block of sequence
    block_type* m_back_block;

    /// pointer to current front element in m_front_block
    value_type* m_front_element;

    /// pointer to current back element in m_back_block
    value_type* m_back_element;

    /// block allocation strategy
    alloc_strategy_type m_alloc_strategy;

    /// block allocation counter
    unsigned_type m_alloc_count;

    /// allocated block identifiers
    bid_deque_type m_bids;

    /// block manager used
    block_manager* m_bm;

    /// number of blocks to prefetch
    unsigned_type m_blocks2prefetch;

public:
    //! \name Constructors/Destructors
    //! \{

    //! Constructs empty sequence with own write and prefetch block pool
    //!
    //! \param D  number of parallel disks, defaulting to the configured number of scratch disks,
    //!           memory consumption will be 2 * D + 2 blocks
    //!           (first and last block, D blocks as write cache, D block for prefetching)
    explicit sequence(int_type D = -1)
        : m_size(0),
          m_owns_pool(true),
          m_alloc_count(0),
          m_bm(block_manager::get_instance())
    {
        if (D < 1) D = config::get_instance()->disks_number();
        STXXL_VERBOSE_SEQUENCE("sequence[" << this << "]::sequence(D)");
        m_pool = new pool_type(D, D + 2);
        init();
    }

    //! Constructs empty sequence with own write and prefetch block pool
    //!
    //! \param w_pool_size  number of blocks in the write pool, must be at least 2, recommended at least 3
    //! \param p_pool_size  number of blocks in the prefetch pool, recommended at least 1
    //! \param blocks2prefetch  defines the number of blocks to prefetch (\c front side),
    //!                          default is number of block in the prefetch pool
    explicit sequence(unsigned_type w_pool_size, unsigned_type p_pool_size, int blocks2prefetch = -1)
        : m_size(0),
          m_owns_pool(true),
          m_alloc_count(0),
          m_bm(block_manager::get_instance())
    {
        STXXL_VERBOSE_SEQUENCE("sequence[" << this << "]::sequence(sizes)");
        m_pool = new pool_type(p_pool_size, w_pool_size);
        init(blocks2prefetch);
    }

    //! Constructs empty sequence
    //!
    //! \param pool block write/prefetch pool
    //! \param blocks2prefetch  defines the number of blocks to prefetch (\c front side), default is number of blocks in the prefetch pool
    //!  \warning Number of blocks in the write pool must be at least 2, recommended at least 3
    //!  \warning Number of blocks in the prefetch pool recommended at least 1
    sequence(pool_type& pool, int blocks2prefetch = -1)
        : m_size(0),
          m_owns_pool(false),
          m_pool(&pool),
          m_alloc_count(0),
          m_bm(block_manager::get_instance())
    {
        STXXL_VERBOSE_SEQUENCE("sequence[" << this << "]::sequence(pool)");
        init(blocks2prefetch);
    }

    //! \}

    //! \name Modifiers
    //! \{

    void swap(sequence& obj)
    {
        std::swap(m_size, obj.m_size);
        std::swap(m_owns_pool, obj.m_owns_pool);
        std::swap(m_pool, obj.m_pool);
        std::swap(m_front_block, obj.m_front_block);
        std::swap(m_back_block, obj.m_back_block);
        std::swap(m_front_element, obj.m_front_element);
        std::swap(m_back_element, obj.m_back_element);
        std::swap(m_alloc_strategy, obj.m_alloc_strategy);
        std::swap(m_alloc_count, obj.m_alloc_count);
        std::swap(m_bids, obj.m_bids);
        std::swap(m_bm, obj.m_bm);
        std::swap(m_blocks2prefetch, obj.m_blocks2prefetch);
    }

    //! \}

private:
    void init(int blocks2prefetch = -1)
    {
        if (m_pool->size_write() < 2) {
            STXXL_ERRMSG("sequence: invalid configuration, not enough blocks (" << m_pool->size_write() <<
                         ") in write pool, at least 2 are needed, resizing to 3");
            m_pool->resize_write(3);
        }

        if (m_pool->size_write() < 3) {
            STXXL_MSG("sequence: inefficient configuration, no blocks for buffered writing available");
        }

        if (m_pool->size_prefetch() < 1) {
            STXXL_MSG("sequence: inefficient configuration, no blocks for prefetching available");
        }

        /// initialize empty sequence
        m_front_block = m_back_block = m_pool->steal();
        m_back_element = m_back_block->begin() - 1;
        m_front_element = m_back_block->begin();
        set_prefetch_aggr(blocks2prefetch);
    }

public:
    //! \name Miscellaneous
    //! \{

    //! Defines the number of blocks to prefetch (\c front side).
    //! This method should be called whenever the prefetch pool is resized
    //! \param blocks2prefetch  defines the number of blocks to prefetch (\c front side),
    //!                         a negative value means to use the number of blocks in the prefetch pool
    void set_prefetch_aggr(int_type blocks2prefetch)
    {
        if (blocks2prefetch < 0)
            m_blocks2prefetch = m_pool->size_prefetch();
        else
            m_blocks2prefetch = blocks2prefetch;
    }

    //! Returns the number of blocks prefetched from the \c front side
    unsigned_type get_prefetch_aggr() const
    {
        return m_blocks2prefetch;
    }

    //! \}

    //! \name Modifiers
    //! \{

    //! Adds an element to the front of the sequence
    void push_front(const value_type& val)
    {
        if (UNLIKELY(m_front_element == m_front_block->begin()))
        {
            if (m_size == 0)
            {
                STXXL_VERBOSE1("sequence::push_front Case 0");
                assert(m_front_block == m_back_block);
                m_front_element = m_back_element = m_front_block->end() - 1;
                *m_front_element = val;
                ++m_size;
                return;
            }
            // front block is completely filled
            else if (m_front_block == m_back_block)
            {
                // can not write the front block because it
                // is the same as the back block, must keep it memory
                STXXL_VERBOSE1("sequence::push_front Case 1");
            }
            else if (size() < 2 * block_type::size)
            {
                STXXL_VERBOSE1("sequence::push_front Case 1.5");
                // only two blocks with a gap at the end, move elements within memory
                assert(m_bids.empty());
                size_t gap = m_back_block->end() - (m_back_element + 1);
                assert(gap > 0);
                std::copy_backward(m_back_block->begin(), m_back_element + 1, m_back_block->end());
                std::copy_backward(m_front_block->end() - gap, m_front_block->end(), m_back_block->begin() + gap);
                std::copy_backward(m_front_block->begin(), m_front_block->end() - gap, m_front_block->end());
                m_front_element += gap;
                m_back_element += gap;

                --m_front_element;
                *m_front_element = val;
                ++m_size;
                return;
            }
            else
            {
                STXXL_VERBOSE1("sequence::push_front Case 2");
                // write the front block
                // need to allocate new block
                bid_type newbid;

                m_bm->new_block(m_alloc_strategy, newbid, m_alloc_count++);

                STXXL_VERBOSE_SEQUENCE("sequence[" << this << "]: push_front block " << m_front_block << " @ " << FMT_BID(newbid));
                m_bids.push_front(newbid);
                m_pool->write(m_front_block, newbid);
                if (m_bids.size() <= m_blocks2prefetch) {
                    STXXL_VERBOSE1("sequence::push Case Hints");
                    m_pool->hint(newbid);
                }
            }

            m_front_block = m_pool->steal();
            m_front_element = m_front_block->end() - 1;
            *m_front_element = val;
            ++m_size;
        }
        else // not at beginning of a block
        {
            --m_front_element;
            *m_front_element = val;
            ++m_size;
        }
    }

    //! Adds an element to the end of the sequence
    void push_back(const value_type& val)
    {
        if (UNLIKELY(m_back_element == m_back_block->begin() + (block_type::size - 1)))
        {
            // back block is completely  filled
            if (m_front_block == m_back_block)
            {
                // can not write the back block because it
                // is the same as the front block, must keep it memory
                STXXL_VERBOSE1("sequence::push_back Case 1");
            }
            else if (size() < 2 * block_type::size)
            {
                STXXL_VERBOSE1("sequence::push_back Case 1.5");
                // only two blocks with a gap in the beginning, move elements within memory
                assert(m_bids.empty());
                size_t gap = m_front_element - m_front_block->begin();
                assert(gap > 0);
                std::copy(m_front_element, m_front_block->end(), m_front_block->begin());
                std::copy(m_back_block->begin(), m_back_block->begin() + gap, m_front_block->begin() + (block_type::size - gap));
                std::copy(m_back_block->begin() + gap, m_back_block->end(), m_back_block->begin());
                m_front_element -= gap;
                m_back_element -= gap;

                ++m_back_element;
                *m_back_element = val;
                ++m_size;
                return;
            }
            else
            {
                STXXL_VERBOSE1("sequence::push_back Case 2");
                // write the back block
                // need to allocate new block
                bid_type newbid;

                m_bm->new_block(m_alloc_strategy, newbid, m_alloc_count++);

                STXXL_VERBOSE_SEQUENCE("sequence[" << this << "]: push_back block " << m_back_block << " @ " << FMT_BID(newbid));
                m_bids.push_back(newbid);
                m_pool->write(m_back_block, newbid);
                if (m_bids.size() <= m_blocks2prefetch) {
                    STXXL_VERBOSE1("sequence::push_back Case Hints");
                    m_pool->hint(newbid);
                }
            }
            m_back_block = m_pool->steal();

            m_back_element = m_back_block->begin();
            *m_back_element = val;
            ++m_size;
        }
        else // not at end of a block
        {
            ++m_back_element;
            *m_back_element = val;
            ++m_size;
        }
    }

    //! Removes element from the front of the sequence
    void pop_front()
    {
        assert(!empty());

        if (UNLIKELY(m_front_element == m_front_block->begin() + (block_type::size - 1)))
        {
            // if there is only one block, it implies ...
            if (m_back_block == m_front_block)
            {
                STXXL_VERBOSE1("sequence::pop_front Case 1");
                assert(size() == 1);
                assert(m_back_element == m_front_element);
                assert(m_bids.empty());
                // reset everything
                m_back_element = m_back_block->begin() - 1;
                m_front_element = m_back_block->begin();
                m_size = 0;
                return;
            }

            --m_size;
            if (m_size <= block_type::size)
            {
                STXXL_VERBOSE1("sequence::pop_front Case 2");
                assert(m_bids.empty());
                // the m_back_block is the next block
                m_pool->add(m_front_block);
                m_front_block = m_back_block;
                m_front_element = m_back_block->begin();
                return;
            }
            STXXL_VERBOSE1("sequence::pop_front Case 3");

            assert(!m_bids.empty());
            request_ptr req = m_pool->read(m_front_block, m_bids.front());
            STXXL_VERBOSE_SEQUENCE("sequence[" << this << "]: pop_front block  " << m_front_block << " @ " << FMT_BID(m_bids.front()));

            // give prefetching hints
            for (unsigned_type i = 0; i < m_blocks2prefetch && i < m_bids.size() - 1; ++i)
            {
                STXXL_VERBOSE1("sequence::pop_front Case Hints");
                m_pool->hint(m_bids[i + 1]);
            }

            m_front_element = m_front_block->begin();
            req->wait();

            m_bm->delete_block(m_bids.front());
            m_bids.pop_front();
        }
        else
        {
            ++m_front_element;
            --m_size;
        }
    }

    //! Removes element from the back of the sequence
    void pop_back()
    {
        assert(!empty());

        if (UNLIKELY(m_back_element == m_back_block->begin()))
        {
            // if there is only one block, it implies ...
            if (m_back_block == m_front_block)
            {
                STXXL_VERBOSE1("sequence::pop_back Case 1");
                assert(size() == 1);
                assert(m_back_element == m_front_element);
                assert(m_bids.empty());
                // reset everything
                m_back_element = m_back_block->begin() - 1;
                m_front_element = m_back_block->begin();
                m_size = 0;
                return;
            }

            --m_size;
            if (m_size <= block_type::size)
            {
                STXXL_VERBOSE1("sequence::pop_back Case 2");
                assert(m_bids.empty());
                // the m_front_block is the next block
                m_pool->add(m_back_block);
                m_back_block = m_front_block;
                m_back_element = m_back_block->end() - 1;
                return;
            }

            STXXL_VERBOSE1("sequence::pop_back Case 3");

            assert(!m_bids.empty());
            request_ptr req = m_pool->read(m_back_block, m_bids.back());
            STXXL_VERBOSE_SEQUENCE("sequence[" << this << "]: pop_back block  " << m_back_block << " @ " << FMT_BID(m_bids.back()));

            // give prefetching hints
            for (unsigned_type i = 1; i < m_blocks2prefetch && i < m_bids.size() - 1; ++i)
            {
                STXXL_VERBOSE1("sequence::pop_front Case Hints");
                m_pool->hint(m_bids[m_bids.size() - 1 - i]);
            }

            m_back_element = m_back_block->end() - 1;
            req->wait();

            m_bm->delete_block(m_bids.back());
            m_bids.pop_back();
        }
        else
        {
            --m_back_element;
            --m_size;
        }
    }

    //! \}

    //! \name Capacity
    //! \{

    //! Returns the size of the sequence
    size_type size() const
    {
        return m_size;
    }

    //! Returns \c true if sequence is empty
    bool empty() const
    {
        return (m_size == 0);
    }

    //! \}

    //! \name Operators
    //! \{

    //! Returns a mutable reference at the back of the sequence
    value_type & back()
    {
        assert(!empty());
        return *m_back_element;
    }

    //! Returns a const reference at the back of the sequence
    const value_type & back() const
    {
        assert(!empty());
        return *m_back_element;
    }

    //! Returns a mutable reference at the front of the sequence
    value_type & front()
    {
        assert(!empty());
        return *m_front_element;
    }

    //! Returns a const reference at the front of the sequence
    const value_type & front() const
    {
        assert(!empty());
        return *m_front_element;
    }

    //! \}

    //! \name Constructors/Destructors
    //! \{

    ~sequence()
    {
        if (m_front_block != m_back_block)
            m_pool->add(m_back_block);

        m_pool->add(m_front_block);

        if (m_owns_pool)
            delete m_pool;

        if (!m_bids.empty())
            m_bm->delete_blocks(m_bids.begin(), m_bids.end());
    }

    //! \}

    /**************************************************************************/

    class stream
    {
    public:
        typedef typename sequence::value_type value_type;

        typedef typename bid_deque_type::const_iterator bid_iter_type;

    protected:
        const sequence& m_sequence;

        size_type m_size;

        value_type* m_current_element;

        block_type* m_current_block;

        bid_iter_type m_next_bid;

    public:
        stream(const sequence& sequence)
            : m_sequence(sequence),
              m_size(sequence.size())
        {
            m_current_block = sequence.m_front_block;
            m_current_element = sequence.m_front_element;
            m_next_bid = sequence.m_bids.begin();
        }

        ~stream()
        {
            if (m_current_block != m_sequence.m_front_block &&
                m_current_block != m_sequence.m_back_block)   // give m_current_block back to pool
                m_sequence.m_pool->add(m_current_block);
        }

        //! return number of element left till end-of-stream.
        size_type size() const
        {
            return m_size;
        }

        //! standard stream method
        bool empty() const
        {
            return (m_size == 0);
        }

        //! standard stream method
        const value_type& operator * () const
        {
            assert(!empty());
            return *m_current_element;
        }

        //! standard stream method
        stream& operator ++ ()
        {
            assert(!empty());

            if (UNLIKELY(m_current_element == m_current_block->begin() + (block_type::size - 1)))
            {
                // next item position is beyond end of current block, find next block
                --m_size;

                if (m_size == 0)
                {
                    STXXL_VERBOSE1("sequence::stream::operator++ last block finished clean at block end");
                    assert(m_next_bid == m_sequence.m_bids.end());
                    assert(m_current_block == m_sequence.m_back_block);
                    // nothing to give back to sequence pool
                    m_current_element = NULL;
                    return *this;
                }
                else if (m_size <= block_type::size)    // still items left in last partial block
                {
                    STXXL_VERBOSE1("sequence::stream::operator++ reached last block");
                    assert(m_next_bid == m_sequence.m_bids.end());
                    // the m_back_block is the next block
                    if (m_current_block != m_sequence.m_front_block) // give current_block back to pool
                        m_sequence.m_pool->add(m_current_block);
                    m_current_block = m_sequence.m_back_block;
                    m_current_element = m_current_block->begin();
                    return *this;
                }
                else if (m_current_block == m_sequence.m_front_block)
                {
                    STXXL_VERBOSE1("sequence::stream::operator++ first own-block case: steal block from sequence's pool");
                    m_current_block = m_sequence.m_pool->steal();
                }

                STXXL_VERBOSE1("sequence::stream::operator++ default case: fetch next block");

                assert(m_next_bid != m_sequence.m_bids.end());
                request_ptr req = m_sequence.m_pool->read(m_current_block, *m_next_bid);
                STXXL_VERBOSE_SEQUENCE("sequence[" << this << "]::stream::operator++ read block " << m_current_block << " @ " << FMT_BID(*m_next_bid));

                // give prefetching hints
                bid_iter_type bid = m_next_bid + 1;
                for (unsigned_type i = 0; i < m_sequence.m_blocks2prefetch && bid != m_sequence.m_bids.end(); ++i, ++bid)
                {
                    STXXL_VERBOSE1("sequence::stream::operator++ giving prefetch hints");
                    m_sequence.m_pool->hint(*bid);
                }

                m_current_element = m_current_block->begin();
                req->wait();

                ++m_next_bid;
            }
            else
            {
                --m_size;
                ++m_current_element;
            }
            return *this;
        }
    };

    //! \name Miscellaneous
    //! \{

    //! construct a forward stream from this sequence
    stream get_stream()
    {
        return stream(*this);
    }

    //! \}

    /**************************************************************************/

    class reverse_stream
    {
    public:
        typedef typename sequence::value_type value_type;

        typedef typename bid_deque_type::const_reverse_iterator bid_iter_type;

    protected:
        const sequence& m_sequence;

        size_type m_size;

        value_type* m_current_element;

        block_type* m_current_block;

        bid_iter_type m_next_bid;

    public:
        reverse_stream(const sequence& sequence)
            : m_sequence(sequence),
              m_size(sequence.size())
        {
            m_current_block = sequence.m_back_block;
            m_current_element = sequence.m_back_element;
            m_next_bid = sequence.m_bids.rbegin();
        }

        ~reverse_stream()
        {
            if (m_current_block != m_sequence.m_front_block &&
                m_current_block != m_sequence.m_back_block)   // give m_current_block back to pool
                m_sequence.m_pool->add(m_current_block);
        }

        //! return number of element left till end-of-stream.
        size_type size() const
        {
            return m_size;
        }

        //! standard stream method
        bool empty() const
        {
            return (m_size == 0);
        }

        //! standard stream method
        const value_type& operator * () const
        {
            assert(!empty());
            return *m_current_element;
        }

        //! standard stream method
        reverse_stream& operator ++ ()
        {
            assert(!empty());

            if (UNLIKELY(m_current_element == m_current_block->begin()))
            {
                // next item position is beyond begin of current block, find next block
                --m_size;

                if (m_size == 0)
                {
                    STXXL_VERBOSE1("sequence::reverse_stream::operator++ last block finished clean at block begin");
                    assert(m_next_bid == m_sequence.m_bids.rend());
                    assert(m_current_block == m_sequence.m_front_block);
                    // nothing to give back to sequence pool
                    m_current_element = NULL;
                    return *this;
                }
                else if (m_size <= block_type::size)
                {
                    STXXL_VERBOSE1("sequence::reverse_stream::operator++ reached first block");
                    assert(m_next_bid == m_sequence.m_bids.rend());
                    // the m_back_block is the next block
                    if (m_current_block != m_sequence.m_back_block) // give current_block back to pool
                        m_sequence.m_pool->add(m_current_block);
                    m_current_block = m_sequence.m_front_block;
                    m_current_element = m_current_block->begin() + (block_type::size - 1);
                    return *this;
                }
                else if (m_current_block == m_sequence.m_back_block)
                {
                    STXXL_VERBOSE1("sequence::reverse_stream::operator++ first own-block case: steal block from sequence's pool");
                    m_current_block = m_sequence.m_pool->steal();
                }

                STXXL_VERBOSE1("sequence::reverse_stream::operator++ default case: fetch previous block");

                assert(m_next_bid != m_sequence.m_bids.rend());
                request_ptr req = m_sequence.m_pool->read(m_current_block, *m_next_bid);
                STXXL_VERBOSE_SEQUENCE("sequence[" << this << "]::reverse_stream::operator++ read block " << m_current_block << " @ " << FMT_BID(*m_next_bid));

                // give prefetching hints
                bid_iter_type bid = m_next_bid + 1;
                for (unsigned_type i = 0; i < m_sequence.m_blocks2prefetch && bid != m_sequence.m_bids.rend(); ++i, ++bid)
                {
                    STXXL_VERBOSE1("sequence::reverse_stream::operator++ giving prefetch hints");
                    m_sequence.m_pool->hint(*bid);
                }

                m_current_element = m_current_block->begin() + (block_type::size - 1);
                req->wait();

                ++m_next_bid;
            }
            else
            {
                --m_size;
                --m_current_element;
            }
            return *this;
        }
    };

    //! \name Miscellaneous
    //! \{

    //! construct a reverse stream from this sequence
    reverse_stream get_reverse_stream()
    {
        return reverse_stream(*this);
    }

    //! \}
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_SEQUENCE_HEADER
// vim: et:ts=4:sw=4
