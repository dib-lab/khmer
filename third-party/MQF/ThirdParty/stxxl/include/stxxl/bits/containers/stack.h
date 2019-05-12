/***************************************************************************
 *  include/stxxl/bits/containers/stack.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003-2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_STACK_HEADER
#define STXXL_CONTAINERS_STACK_HEADER

#include <stack>
#include <vector>

#include <stxxl/bits/deprecated.h>
#include <stxxl/bits/io/request_operations.h>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/common/simple_vector.h>
#include <stxxl/bits/common/tmeta.h>
#include <stxxl/bits/common/error_handling.h>
#include <stxxl/bits/mng/read_write_pool.h>
#include <stxxl/bits/mng/write_pool.h>
#include <stxxl/bits/mng/prefetch_pool.h>

STXXL_BEGIN_NAMESPACE

//! \defgroup stlcont_stack stack
//! \ingroup stlcont
//! External stack implementations
//! \{

template <class ValueType,
          unsigned BlocksPerPage = 4,
          unsigned BlockSize = STXXL_DEFAULT_BLOCK_SIZE(ValueType),
          class AllocStr = STXXL_DEFAULT_ALLOC_STRATEGY,
          class SizeType = stxxl::uint64>
struct stack_config_generator
{
    typedef ValueType value_type;
    enum { blocks_per_page = BlocksPerPage };
    typedef AllocStr alloc_strategy;
    enum { block_size = BlockSize };
    typedef SizeType size_type;
};

//! External stack container.
//! <b> Introduction </b> to stack container: see \ref tutorial_stack tutorial. \n
//! <b> Design and Internals </b> of stack container: see \ref design_stack
//! Conservative implementation. Fits best if your access pattern consists of irregularly mixed
//! push'es and pop's.
//! For semantics of the methods see documentation of the STL \c std::stack. <BR>
//! To gain full bandwidth of disks \c StackConfig::BlocksPerPage must >= number of disks <BR>
//! \internal
template <class StackConfig>
class normal_stack : private noncopyable
{
public:
    typedef StackConfig cfg;
    //! type of the elements stored in the stack
    typedef typename cfg::value_type value_type;
    typedef typename cfg::alloc_strategy alloc_strategy_type;
    //! type for sizes (64-bit)
    typedef typename cfg::size_type size_type;
    enum {
        blocks_per_page = cfg::blocks_per_page,
        block_size = cfg::block_size
    };

    //! type of block used in disk-memory transfers
    typedef typed_block<block_size, value_type> block_type;
    typedef BID<block_size> bid_type;

private:
    size_type m_size;
    unsigned_type cache_offset;
    value_type* current_element;
    simple_vector<block_type> cache;
    typename simple_vector<block_type>::iterator front_page;
    typename simple_vector<block_type>::iterator back_page;
    std::vector<bid_type> bids;
    alloc_strategy_type alloc_strategy;

public:
    //! \name Constructors/Destructors
    //! \{

    //! Default constructor: creates empty stack.
    normal_stack()
        : m_size(0),
          cache_offset(0),
          current_element(NULL),
          cache(blocks_per_page * 2),
          front_page(cache.begin() + blocks_per_page),
          back_page(cache.begin()),
          bids(0)
    {
        bids.reserve(blocks_per_page);
    }

    //! \name Accessor Functions
    //! \{

    void swap(normal_stack& obj)
    {
        std::swap(m_size, obj.m_size);
        std::swap(cache_offset, obj.cache_offset);
        std::swap(current_element, obj.current_element);
        std::swap(cache, obj.cache);
        std::swap(front_page, obj.front_page);
        std::swap(back_page, obj.back_page);
        std::swap(bids, obj.bids);
        std::swap(alloc_strategy, obj.alloc_strategy);
    }

    //! \}

    //! \name Constructors/Destructors
    //! \{

    //! Copy-construction from a another stack of any type.
    //! \param stack_ stack object (could be external or internal, important is that it must
    //! have a copy constructor, \c top() and \c pop() methods )
    template <class StackType>
    normal_stack(const StackType& stack_)
        : m_size(0),
          cache_offset(0),
          current_element(NULL),
          cache(blocks_per_page * 2),
          front_page(cache.begin() + blocks_per_page),
          back_page(cache.begin()),
          bids(0)
    {
        bids.reserve(blocks_per_page);

        StackType stack_copy = stack_;
        size_t sz = stack_copy.size();

        std::vector<value_type> tmp(sz);

        for (size_t i = 0; i < sz; ++i) {
            tmp[sz - i - 1] = stack_copy.top();
            stack_copy.pop();
        }

        for (size_t i = 0; i < sz; ++i)
            push(tmp[i]);
    }

    virtual ~normal_stack()
    {
        STXXL_VERBOSE(STXXL_PRETTY_FUNCTION_NAME);
        block_manager::get_instance()->delete_blocks(bids.begin(), bids.end());
    }

    //! \}

    //! \name Capacity
    //! \{

    //! Returns the number of elements contained in the stack
    size_type size() const
    {
        return m_size;
    }

    //! Returns true if the stack is empty.
    bool empty() const
    {
        return (!m_size);
    }

    //! \}

    //! \name Accessor Functions
    //! \{

    //! Return mutable reference to the element at the top of the
    //! stack. Precondition: stack is not empty().
    value_type & top()
    {
        assert(m_size > 0);
        return (*current_element);
    }

    //! Return constant reference to the element at the top of the
    //! stack. Precondition: stack is not empty().
    const value_type & top() const
    {
        assert(m_size > 0);
        return (*current_element);
    }

    //! Inserts an element at the top of the stack. Postconditions: size() is
    //! incremented by 1, and top() is the inserted element.
    void push(const value_type& val)
    {
        assert(cache_offset <= 2 * blocks_per_page * block_type::size);
        //assert(cache_offset >= 0);

        if (UNLIKELY(cache_offset == 2 * blocks_per_page * block_type::size))  // cache overflow
        {
            STXXL_VERBOSE2("growing, size: " << m_size);

            bids.resize(bids.size() + blocks_per_page);
            typename std::vector<bid_type>::iterator cur_bid = bids.end() - blocks_per_page;
            block_manager::get_instance()->new_blocks(alloc_strategy, cur_bid, bids.end(), cur_bid - bids.begin());

            simple_vector<request_ptr> requests(blocks_per_page);

            for (int i = 0; i < blocks_per_page; ++i, ++cur_bid)
            {
                requests[i] = (back_page + i)->write(*cur_bid);
            }

            std::swap(back_page, front_page);

            bids.reserve(bids.size() + blocks_per_page);

            cache_offset = blocks_per_page * block_type::size + 1;
            current_element = &((*front_page)[0]);
            ++m_size;

            wait_all(requests.begin(), blocks_per_page);

            *current_element = val;

            return;
        }

        current_element = element(cache_offset);
        *current_element = val;
        ++m_size;
        ++cache_offset;
    }

    //! Removes the element at the top of the stack. Precondition: stack is not
    //! empty(). Postcondition: size() is decremented.
    void pop()
    {
        assert(cache_offset <= 2 * blocks_per_page * block_type::size);
        assert(cache_offset > 0);
        assert(m_size > 0);

        if (UNLIKELY(cache_offset == 1 && bids.size() >= blocks_per_page))
        {
            STXXL_VERBOSE2("shrinking, size: " << m_size);

            simple_vector<request_ptr> requests(blocks_per_page);

            {
                typename std::vector<bid_type>::const_iterator cur_bid = bids.end();
                for (int i = blocks_per_page - 1; i >= 0; --i)
                {
                    requests[i] = (front_page + i)->read(*(--cur_bid));
                }
            }

            std::swap(front_page, back_page);

            cache_offset = blocks_per_page * block_type::size;
            --m_size;
            current_element = &((*(back_page + (blocks_per_page - 1)))[block_type::size - 1]);

            wait_all(requests.begin(), blocks_per_page);

            block_manager::get_instance()->delete_blocks(bids.end() - blocks_per_page, bids.end());
            bids.resize(bids.size() - blocks_per_page);

            return;
        }

        --m_size;

        current_element = element((--cache_offset) - 1);
    }

    //! \}

private:
    value_type * element(unsigned_type offset)
    {
        if (offset < blocks_per_page * block_type::size)
            return &((*(back_page + offset / block_type::size))[offset % block_type::size]);

        unsigned_type unbiased_offset = offset - blocks_per_page * block_type::size;
        return &((*(front_page + unbiased_offset / block_type::size))[unbiased_offset % block_type::size]);
    }
};

//! Efficient implementation that uses prefetching and overlapping using internal buffers.
//!
//! Use it if your access pattern consists of many repeated push'es and pop's
//! For semantics of the methods see documentation of the STL \c std::stack.
//! \warning The amortized complexity of operation is not O(1/DB), rather O(DB)
template <class StackConfig>
class grow_shrink_stack : private noncopyable
{
public:
    typedef StackConfig cfg;
    //! type of the elements stored in the stack
    typedef typename cfg::value_type value_type;
    typedef typename cfg::alloc_strategy alloc_strategy_type;
    //! type for sizes (64-bit)
    typedef typename cfg::size_type size_type;
    enum {
        blocks_per_page = cfg::blocks_per_page,
        block_size = cfg::block_size
    };

    //! type of block used in disk-memory transfers
    typedef typed_block<block_size, value_type> block_type;
    typedef BID<block_size> bid_type;

private:
    size_type m_size;
    unsigned_type cache_offset;
    value_type* current_element;
    simple_vector<block_type> cache;
    typename simple_vector<block_type>::iterator cache_buffers;
    typename simple_vector<block_type>::iterator overlap_buffers;
    simple_vector<request_ptr> requests;
    std::vector<bid_type> bids;
    alloc_strategy_type alloc_strategy;

public:
    //! \name Constructors/Destructor
    //! \{

    //! Default constructor: creates empty stack.
    grow_shrink_stack()
        : m_size(0),
          cache_offset(0),
          current_element(NULL),
          cache(blocks_per_page * 2),
          cache_buffers(cache.begin()),
          overlap_buffers(cache.begin() + blocks_per_page),
          requests(blocks_per_page),
          bids(0)
    {
        bids.reserve(blocks_per_page);
    }

    //! \}

    //! \name Accessor Functions
    //! \{

    void swap(grow_shrink_stack& obj)
    {
        std::swap(m_size, obj.m_size);
        std::swap(cache_offset, obj.cache_offset);
        std::swap(current_element, obj.current_element);
        std::swap(cache, obj.cache);
        std::swap(cache_buffers, obj.cache_buffers);
        std::swap(overlap_buffers, obj.overlap_buffers);
        std::swap(requests, obj.requests);
        std::swap(bids, obj.bids);
        std::swap(alloc_strategy, obj.alloc_strategy);
    }

    //! \}

    //! \name Constructors/Destructors
    //! \{

    //! Copy-construction from a another stack of any type.
    //! \param stack_ stack object (could be external or internal, important is that it must
    //! have a copy constructor, \c top() and \c pop() methods )
    template <class StackType>
    grow_shrink_stack(const StackType& stack_)
        : m_size(0),
          cache_offset(0),
          current_element(NULL),
          cache(blocks_per_page * 2),
          cache_buffers(cache.begin()),
          overlap_buffers(cache.begin() + blocks_per_page),
          requests(blocks_per_page),
          bids(0)
    {
        bids.reserve(blocks_per_page);

        StackType stack_copy = stack_;
        size_t sz = stack_copy.size();

        std::vector<value_type> tmp(sz);

        for (size_t i = 0; i < sz; ++i)
        {
            tmp[sz - i - 1] = stack_copy.top();
            stack_copy.pop();
        }

        for (size_t i = 0; i < sz; ++i)
            push(tmp[i]);
    }
    virtual ~grow_shrink_stack()
    {
        STXXL_VERBOSE(STXXL_PRETTY_FUNCTION_NAME);
        try
        {
            if (requests[0].get())
                wait_all(requests.begin(), blocks_per_page);
        }
        catch (const io_error&)
        { }
        block_manager::get_instance()->delete_blocks(bids.begin(), bids.end());
    }

    //! \}

    //! \name Capacity
    //! \{

    //! Returns the number of elements contained in the stack
    size_type size() const
    {
        return m_size;
    }
    //! Returns true if the stack is empty.
    bool empty() const
    {
        return (!m_size);
    }

    //! \}

    //! \name Accessor Functions
    //! \{

    //! Return mutable reference to the element at the top of the
    //! stack. Precondition: stack is not empty().
    value_type & top()
    {
        assert(m_size > 0);
        return (*current_element);
    }

    //! Return constant reference to the element at the top of the
    //! stack. Precondition: stack is not empty().
    const value_type & top() const
    {
        assert(m_size > 0);
        return (*current_element);
    }

    //! Inserts an element at the top of the stack. Postconditions: size() is
    //! incremented by 1, and top() is the inserted element.
    void push(const value_type& val)
    {
        assert(cache_offset <= blocks_per_page * block_type::size);
        //assert(cache_offset >= 0);

        if (UNLIKELY(cache_offset == blocks_per_page * block_type::size))  // cache overflow
        {
            STXXL_VERBOSE2("growing, size: " << m_size);

            bids.resize(bids.size() + blocks_per_page);
            typename std::vector<bid_type>::iterator cur_bid = bids.end() - blocks_per_page;
            block_manager::get_instance()->new_blocks(alloc_strategy, cur_bid, bids.end(), cur_bid - bids.begin());

            for (int i = 0; i < blocks_per_page; ++i, ++cur_bid)
            {
                if (requests[i].get())
                    requests[i]->wait();

                requests[i] = (cache_buffers + i)->write(*cur_bid);
            }

            std::swap(cache_buffers, overlap_buffers);

            bids.reserve(bids.size() + blocks_per_page);

            cache_offset = 1;
            current_element = &((*cache_buffers)[0]);
            ++m_size;

            *current_element = val;

            return;
        }

        current_element = &((*(cache_buffers + cache_offset / block_type::size))[cache_offset % block_type::size]);
        *current_element = val;
        ++m_size;
        ++cache_offset;
    }

    //! Removes the element at the top of the stack. Precondition: stack is not
    //! empty(). Postcondition: size() is decremented.
    void pop()
    {
        assert(cache_offset <= blocks_per_page * block_type::size);
        assert(cache_offset > 0);
        assert(m_size > 0);

        if (UNLIKELY(cache_offset == 1 && bids.size() >= blocks_per_page))
        {
            STXXL_VERBOSE2("shrinking, size: " << m_size);

            if (requests[0].get())
                wait_all(requests.begin(), blocks_per_page);

            std::swap(cache_buffers, overlap_buffers);

            if (bids.size() > blocks_per_page)
            {
                STXXL_VERBOSE2("prefetching, size: " << m_size);
                typename std::vector<bid_type>::const_iterator cur_bid = bids.end() - blocks_per_page;
                for (int i = blocks_per_page - 1; i >= 0; --i)
                    requests[i] = (overlap_buffers + i)->read(*(--cur_bid));
            }

            block_manager::get_instance()->delete_blocks(bids.end() - blocks_per_page, bids.end());
            bids.resize(bids.size() - blocks_per_page);

            cache_offset = blocks_per_page * block_type::size;
            --m_size;
            current_element = &((*(cache_buffers + (blocks_per_page - 1)))[block_type::size - 1]);

            return;
        }

        --m_size;
        unsigned_type cur_offset = (--cache_offset) - 1;
        current_element = &((*(cache_buffers + cur_offset / block_type::size))[cur_offset % block_type::size]);
    }

    //! \}
};

//! Efficient implementation that uses prefetching and overlapping using (shared) buffers pools.
//! \warning This is a single buffer stack! Each direction change (push() followed by pop() or vice versa) may cause one I/O.
template <class StackConfig>
class grow_shrink_stack2 : private noncopyable
{
public:
    typedef StackConfig cfg;
    //! type of the elements stored in the stack
    typedef typename cfg::value_type value_type;
    typedef typename cfg::alloc_strategy alloc_strategy_type;
    //! type for sizes (64-bit)
    typedef typename cfg::size_type size_type;
    enum {
        blocks_per_page = cfg::blocks_per_page,     // stack of this type has only one page
        block_size = cfg::block_size
    };

    //! type of block used in disk-memory transfers
    typedef typed_block<block_size, value_type> block_type;
    typedef BID<block_size> bid_type;

private:
    typedef read_write_pool<block_type> pool_type;

    size_type m_size;
    unsigned_type cache_offset;
    block_type* cache;
    std::vector<bid_type> bids;
    alloc_strategy_type alloc_strategy;
    unsigned_type pref_aggr;
    pool_type* owned_pool;
    pool_type* pool;

public:
    //! \name Constructors/Destructors
    //! \{

    //! Default constructor: creates empty stack. The stack will use the
    //! read_write_pool for prefetching and buffered writing.
    //! \param pool_ block write/prefetch pool
    //! \param prefetch_aggressiveness number of blocks that will be used from prefetch pool
    grow_shrink_stack2(pool_type& pool_,
                       unsigned_type prefetch_aggressiveness = 0)
        : m_size(0),
          cache_offset(0),
          cache(new block_type),
          pref_aggr(prefetch_aggressiveness),
          owned_pool(NULL),
          pool(&pool_)
    {
        STXXL_VERBOSE2("grow_shrink_stack2::grow_shrink_stack2(...)");
    }

    //! Default constructor: creates empty stack. The stack will use the pair
    //! of prefetch_pool and write_pool for prefetching and buffered writing.
    //! This constructor is deprecated in favor of the read_write_pool
    //! constructor.
    //!
    //! \param p_pool_ prefetch pool, that will be used for block prefetching
    //! \param w_pool_ write pool, that will be used for block writing
    //! \param prefetch_aggressiveness number of blocks that will be used from prefetch pool
    STXXL_DEPRECATED(
        grow_shrink_stack2(prefetch_pool<block_type>& p_pool_,
                           write_pool<block_type>& w_pool_,
                           unsigned_type prefetch_aggressiveness = 0)
        )
        : m_size(0),
          cache_offset(0),
          cache(new block_type),
          pref_aggr(prefetch_aggressiveness),
          owned_pool(new pool_type(p_pool_, w_pool_)),
          pool(owned_pool)
    {
        STXXL_VERBOSE2("grow_shrink_stack2::grow_shrink_stack2(...)");
    }

    //! \}

    //! \name Accessor Functions
    //! \{

    void swap(grow_shrink_stack2& obj)
    {
        std::swap(m_size, obj.m_size);
        std::swap(cache_offset, obj.cache_offset);
        std::swap(cache, obj.cache);
        std::swap(bids, obj.bids);
        std::swap(alloc_strategy, obj.alloc_strategy);
        std::swap(pref_aggr, obj.pref_aggr);
        std::swap(owned_pool, obj.owned_pool);
        std::swap(pool, obj.pool);
    }

    //! \}

    //! \name Constructors/Destructors
    //! \{

    virtual ~grow_shrink_stack2()
    {
        try
        {
            STXXL_VERBOSE2("grow_shrink_stack2::~grow_shrink_stack2()");
            const int_type bids_size = bids.size();
            const int_type last_pref = STXXL_MAX(int_type(bids_size) - int_type(pref_aggr), (int_type)0);
            int_type i;
            for (i = bids_size - 1; i >= last_pref; --i)
            {
                // clean the prefetch buffer
                pool->invalidate(bids[i]);
            }
            typename std::vector<bid_type>::iterator cur = bids.begin();
            typename std::vector<bid_type>::const_iterator end = bids.end();
            for ( ; cur != end; ++cur)
            {
                // FIXME: read_write_pool needs something like cancel_write(bid)
                block_type* b = NULL;   // w_pool.steal(*cur);
                if (b)
                {
                    pool->add(cache);   // return buffer
                    cache = b;
                }
            }
            delete cache;
        }
        catch (const io_error&)
        { }
        block_manager::get_instance()->delete_blocks(bids.begin(), bids.end());
        delete owned_pool;
    }

    //! \}

    //! \name Capacity
    //! \{

    //! Returns the number of elements contained in the stack
    size_type size() const
    {
        return m_size;
    }
    //! Returns true if the stack is empty.
    bool empty() const
    {
        return (!m_size);
    }

    //! \}

    //! \name Accessor Functions
    //! \{

    //! Inserts an element at the top of the stack. Postconditions: size() is
    //! incremented by 1, and top() is the inserted element.
    void push(const value_type& val)
    {
        STXXL_VERBOSE3("grow_shrink_stack2::push(" << val << ")");
        assert(cache_offset <= block_type::size);

        if (UNLIKELY(cache_offset == block_type::size))
        {
            STXXL_VERBOSE2("grow_shrink_stack2::push(" << val << ") growing, size: " << m_size);

            bids.resize(bids.size() + 1);
            typename std::vector<bid_type>::iterator cur_bid = bids.end() - 1;
            block_manager::get_instance()->new_blocks(alloc_strategy, cur_bid, bids.end(), cur_bid - bids.begin());
            pool->write(cache, bids.back());
            cache = pool->steal();
            const int_type bids_size = bids.size();
            const int_type last_pref = STXXL_MAX(int_type(bids_size) - int_type(pref_aggr) - 1, (int_type)0);
            for (int_type i = bids_size - 2; i >= last_pref; --i)
            {
                // clean prefetch buffers
                pool->invalidate(bids[i]);
            }
            cache_offset = 0;
        }
        (*cache)[cache_offset] = val;
        ++m_size;
        ++cache_offset;

        assert(cache_offset > 0);
        assert(cache_offset <= block_type::size);
    }

    //! Return mutable reference to the element at the top of the
    //! stack. Precondition: stack is not empty().
    value_type & top()
    {
        assert(m_size > 0);
        assert(cache_offset > 0);
        assert(cache_offset <= block_type::size);
        return (*cache)[cache_offset - 1];
    }

    //! Return constant reference to the element at the top of the
    //! stack. Precondition: stack is not empty().
    const value_type & top() const
    {
        assert(m_size > 0);
        assert(cache_offset > 0);
        assert(cache_offset <= block_type::size);
        return (*cache)[cache_offset - 1];
    }

    //! Removes the element at the top of the stack. Precondition: stack is not
    //! empty(). Postcondition: size() is decremented.
    void pop()
    {
        STXXL_VERBOSE3("grow_shrink_stack2::pop()");
        assert(m_size > 0);
        assert(cache_offset > 0);
        assert(cache_offset <= block_type::size);
        if (UNLIKELY(cache_offset == 1 && (!bids.empty())))
        {
            STXXL_VERBOSE2("grow_shrink_stack2::pop() shrinking, size = " << m_size);

            bid_type last_block = bids.back();
            bids.pop_back();
            pool->read(cache, last_block)->wait();
            block_manager::get_instance()->delete_block(last_block);
            rehint();
            cache_offset = block_type::size + 1;
        }

        --cache_offset;
        --m_size;
    }

    //! \}

    //! \name Miscellaneous
    //! \{

    //! Sets level of prefetch aggressiveness (number of blocks from the
    //! prefetch pool used for prefetching).
    //! \param new_p new value for the prefetch aggressiveness
    void set_prefetch_aggr(unsigned_type new_p)
    {
        if (pref_aggr > new_p)
        {
            const int_type bids_size = bids.size();
            const int_type last_pref = STXXL_MAX(int_type(bids_size) - int_type(pref_aggr), (int_type)0);
            for (int_type i = bids_size - new_p - 1; i >= last_pref; --i)
            {
                // clean prefetch buffers
                pool->invalidate(bids[i]);
            }
        }
        pref_aggr = new_p;
        rehint();
    }

    //! Returns number of blocks used for prefetching.
    unsigned_type get_prefetch_aggr() const
    {
        return pref_aggr;
    }

    //! \}

private:
    //! hint the last pref_aggr external blocks.
    void rehint()
    {
        const int_type bids_size = bids.size();
        const int_type last_pref = STXXL_MAX(int_type(bids_size) - int_type(pref_aggr), (int_type)0);
        for (int_type i = bids_size - 1; i >= last_pref; --i)
        {
            pool->hint(bids[i]);  // prefetch
        }
    }
};

//! A stack that migrates from internal memory to external when its size exceeds a certain threshold.
//!
//! For semantics of the methods see documentation of the STL \c std::stack.
template <unsigned_type CritSize, class ExternalStack, class InternalStack>
class migrating_stack : private noncopyable
{
public:
    typedef typename ExternalStack::cfg cfg;
    //! type of the elements stored in the stack
    typedef typename cfg::value_type value_type;
    //! type for sizes (64-bit)
    typedef typename cfg::size_type size_type;
    enum {
        blocks_per_page = cfg::blocks_per_page,
        block_size = cfg::block_size
    };

    typedef InternalStack int_stack_type;
    typedef ExternalStack ext_stack_type;

private:
    enum { critical_size = CritSize };

    int_stack_type* int_impl;
    ext_stack_type* ext_impl;

    //! Copy-construction from a another stack of any type.
    //! \warning not implemented yet!
    template <class StackType>
    migrating_stack(const StackType& stack_);

public:
    //! \name Constructors/Destructors
    //! \{

    //! Default constructor: creates empty stack.
    migrating_stack()
        : int_impl(new int_stack_type()), ext_impl(NULL)
    { }

    virtual ~migrating_stack()
    {
        delete int_impl;
        delete ext_impl;
    }

    //! \}

    //! \name Accessor Functions
    //! \{

    void swap(migrating_stack& obj)
    {
        std::swap(int_impl, obj.int_impl);
        std::swap(ext_impl, obj.ext_impl);
    }

    //! \}

    //! \name Miscellaneous
    //! \{

    //! Returns true if current implementation is internal, otherwise false.
    bool internal() const
    {
        assert((int_impl && !ext_impl) || (!int_impl && ext_impl));
        return (int_impl != NULL);
    }
    //! Returns true if current implementation is external, otherwise false.
    bool external() const
    {
        assert((int_impl && !ext_impl) || (!int_impl && ext_impl));
        return (ext_impl != NULL);
    }

    //! \}

    //! \name Capacity
    //! \{

    //! Returns true if the stack is empty.
    bool empty() const
    {
        assert((int_impl && !ext_impl) || (!int_impl && ext_impl));
        return (int_impl) ? int_impl->empty() : ext_impl->empty();
    }
    //! Returns the number of elements contained in the stack
    size_type size() const
    {
        assert((int_impl && !ext_impl) || (!int_impl && ext_impl));
        return (int_impl) ? size_type(int_impl->size()) : ext_impl->size();
    }

    //! \}

    //! \name Accessor Functions
    //! \{

    //! Return mutable reference to the element at the top of the
    //! stack. Precondition: stack is not empty().
    value_type & top()
    {
        assert((int_impl && !ext_impl) || (!int_impl && ext_impl));
        return (int_impl) ? int_impl->top() : ext_impl->top();
    }

    //! Return constant reference to the element at the top of the
    //! stack. Precondition: stack is not empty().
    const value_type & top() const
    {
        assert((int_impl && !ext_impl) || (!int_impl && ext_impl));
        return (int_impl) ? int_impl->top() : ext_impl->top();
    }

    //! Inserts an element at the top of the stack. Postconditions: size() is
    //! incremented by 1, and top() is the inserted element.
    void push(const value_type& val)
    {
        assert((int_impl && !ext_impl) || (!int_impl && ext_impl));

        if (int_impl)
        {
            int_impl->push(val);
            if (UNLIKELY(int_impl->size() == critical_size))
            {
                // migrate to external stack
                ext_impl = new ext_stack_type(*int_impl);
                delete int_impl;
                int_impl = NULL;
            }
        }
        else
            ext_impl->push(val);
    }

    //! Removes the element at the top of the stack. Precondition: stack is not
    //! empty(). Postcondition: size() is decremented.
    void pop()
    {
        assert((int_impl && !ext_impl) || (!int_impl && ext_impl));

        if (int_impl)
            int_impl->pop();
        else
            ext_impl->pop();
    }

    //! \}
};

enum stack_externality { external, migrating, internal };
enum stack_behaviour { normal, grow_shrink, grow_shrink2 };

//! Stack type generator \n
//! <b> Introduction </b> to stack container: see \ref tutorial_stack tutorial. \n
//! <b> Design and Internals </b> of stack container: see \ref design_stack.
//!
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//!
//! \tparam Externality selects stack implementation, default: \b external. One of
//!   - \c external, external container, implementation is chosen according to \c Behaviour parameter.
//!   - \c migrating, migrates from internal implementation given by \c IntStackType parameter
//!     to external implementation given by \c Behaviour parameter when size exceeds \c MigrCritSize
//!   - \c internal, choses \c IntStackType implementation
//!
//! \tparam Behaviour chooses \b external implementation, default: \b stxxl::normal_stack. One of:
//!   - \c normal, conservative version, implemented in \c stxxl::normal_stack
//!   - \c grow_shrink, efficient version, implemented in \c stxxl::grow_shrink_stack
//!   - \c grow_shrink2, efficient version, implemented in \c stxxl::grow_shrink_stack2
//!
//! \tparam BlocksPerPage defines how many blocks has one page of internal cache of an
//!      \b external implementation, default is \b 4. All \b external implementations have
//!      \b two pages.
//!
//! \tparam BlockSize external block size in bytes, default is <b>2 MiB</b>.
//!
//! \tparam IntStackType type of internal stack used for some implementations, default: \b std::stack.
//!
//! \tparam MigrCritSize threshold value for number of elements when
//!   stxxl::migrating_stack migrates to the external memory, default: <b>2 x BlocksPerPage x BlockSize</b>.
//!
//! \tparam AllocStr one of allocation strategies: striping, RC, SR, or FR. Default is \b RC.
//!
//! \tparam SizeType size type, default is \b stxxl::uint64.
//!
//! The configured stack type is available as STACK_GENERATOR<>::result.
//!
template <
    class ValueType,
    stack_externality Externality = external,
    stack_behaviour Behaviour = normal,
    unsigned BlocksPerPage = 4,
    unsigned BlockSize = STXXL_DEFAULT_BLOCK_SIZE(ValueType),

    class IntStackType = std::stack<ValueType>,
    unsigned_type MigrCritSize = (2* BlocksPerPage* BlockSize),

    class AllocStr = STXXL_DEFAULT_ALLOC_STRATEGY,
    class SizeType = stxxl::uint64
    >
class STACK_GENERATOR
{
    typedef stack_config_generator<ValueType, BlocksPerPage, BlockSize, AllocStr, SizeType> cfg;

    typedef typename IF<Behaviour == grow_shrink,
                        grow_shrink_stack<cfg>,
                        grow_shrink_stack2<cfg> >::result GrShrTp;
    typedef typename IF<Behaviour == normal, normal_stack<cfg>, GrShrTp>::result ExtStackType;
    typedef typename IF<Externality == migrating,
                        migrating_stack<MigrCritSize, ExtStackType, IntStackType>, ExtStackType>::result MigrOrNotStackType;

public:
    typedef typename IF<Externality == internal, IntStackType, MigrOrNotStackType>::result result;
};

//! \}

STXXL_END_NAMESPACE

namespace std {

template <class StackConfig>
void swap(stxxl::normal_stack<StackConfig>& a,
          stxxl::normal_stack<StackConfig>& b)
{
    a.swap(b);
}

template <class StackConfig>
void swap(stxxl::grow_shrink_stack<StackConfig>& a,
          stxxl::grow_shrink_stack<StackConfig>& b)
{
    a.swap(b);
}

template <class StackConfig>
void swap(stxxl::grow_shrink_stack2<StackConfig>& a,
          stxxl::grow_shrink_stack2<StackConfig>& b)
{
    a.swap(b);
}

template <stxxl::unsigned_type CritSize, class ExternalStack, class InternalStack>
void swap(stxxl::migrating_stack<CritSize, ExternalStack, InternalStack>& a,
          stxxl::migrating_stack<CritSize, ExternalStack, InternalStack>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_CONTAINERS_STACK_HEADER
// vim: et:ts=4:sw=4
