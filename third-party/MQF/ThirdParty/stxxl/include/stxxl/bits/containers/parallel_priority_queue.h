/***************************************************************************
 *  include/stxxl/bits/containers/parallel_priority_queue.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014-2015 Thomas Keh <thomas.keh@student.kit.edu>
 *  Copyright (C) 2014-2015 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_PARALLEL_PRIORITY_QUEUE_HEADER
#define STXXL_CONTAINERS_PARALLEL_PRIORITY_QUEUE_HEADER

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <list>
#include <utility>
#include <numeric>
#include <vector>

#if STXXL_PARALLEL
    #include <omp.h>
#endif

#if __cplusplus >= 201103L
#define STXXL_MOVE(T) std::move(T)
#else
#define STXXL_MOVE(T) T
#endif

#include <stxxl/bits/common/winner_tree.h>
#include <stxxl/bits/common/custom_stats.h>
#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/common/timer.h>
#include <stxxl/bits/common/is_heap.h>
#include <stxxl/bits/common/swap_vector.h>
#include <stxxl/bits/common/rand.h>
#include <stxxl/bits/config.h>
#include <stxxl/bits/io/request_operations.h>
#include <stxxl/bits/mng/block_alloc.h>
#include <stxxl/bits/mng/buf_ostream.h>
#include <stxxl/bits/mng/prefetch_pool.h>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/read_write_pool.h>
#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/parallel.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/types>

STXXL_BEGIN_NAMESPACE

namespace ppq_local {

/*!
 * A random-access iterator class for block oriented data.  The iterator is
 * intended to be provided by the internal_array and external_array classes
 * and to be used by the multiway_merge algorithm as input iterators.
 *
 * \tparam ValueType the value type
 */
template <class ValueType>
class ppq_iterator
{
public:
    typedef ValueType value_type;
    typedef value_type& reference;
    typedef value_type* pointer;
    typedef ptrdiff_t difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    typedef std::vector<std::pair<pointer, pointer> > block_pointers_type;

protected:
    typedef ppq_iterator self_type;

    //! pointer to a vector of begin/end pointer pairs
    //! They allow access to the data blocks.
    const block_pointers_type* m_block_pointers;

    //! pointer to the current element
    pointer m_current;

    //! index of the current element
    size_t m_index;

    //! index of the current element's block
    size_t m_block_index;

    //! size of each data block
    size_t m_block_items;

public:
    //! default constructor (should not be used directly)
    ppq_iterator()
        : m_block_pointers(NULL)
    { }

    //! constructor
    //!
    //! \param block_pointers   A reference to the properly initialized vector of begin and end pointers.
    //!                         One pair for each block. The pointers should be valid for all blocks that
    //!                         are expected to be accessed with this iterator.
    //! \param block_items      The size of a single block. If there is only one block (e.g. if the iterator
    //!                         belongs to an internal_array), use the total size here.
    //! \param index            The index of the current element (global - index 0 belongs to the first element
    //!                         in the first block, no matter if the values are still valid)
    ppq_iterator(const block_pointers_type* block_pointers, size_t block_items,
                 size_t index)
        : m_block_pointers(block_pointers),
          m_index(index),
          m_block_items(block_items)
    {
        update();
    }

    //! returns the value's index in the internal or external array
    size_t get_index() const
    {
        return m_index;
    }

    reference operator * () const
    {
        assert(m_current);
        return *m_current;
    }

    pointer operator -> () const
    {
        return &(operator * ());
    }

    reference operator [] (difference_type relative_index) const
    {
        const difference_type index = m_index + relative_index;

        const size_t block_index = index / m_block_items;
        const size_t local_index = index % m_block_items;

        assert(block_index < m_block_pointers->size());
        assert((*m_block_pointers)[block_index].first + local_index
               < (*m_block_pointers)[block_index].second);

        return *((*m_block_pointers)[block_index].first + local_index);
    }

    //! prefix-increment operator
    self_type& operator ++ ()
    {
        ++m_index;
        ++m_current;

        if (UNLIKELY(m_current == (*m_block_pointers)[m_block_index].second)) {
            if (m_block_index + 1 < m_block_pointers->size()) {
                m_current = (*m_block_pointers)[++m_block_index].first;
            }
            else {
                // global end
                assert(m_block_index + 1 == m_block_pointers->size());
                m_current = (*m_block_pointers)[m_block_index++].second;
            }
        }

        return *this;
    }
    //! prefix-decrement operator
    self_type& operator -- ()
    {
        assert(m_index > 0);
        --m_index;

        if (m_block_index >= m_block_pointers->size()
            || m_current == (*m_block_pointers)[m_block_index].first) {
            // begin of current block or global end
            assert(m_block_index > 0);
            assert(m_block_index <= m_block_pointers->size());
            m_current = (*m_block_pointers)[--m_block_index].second - 1;
        }
        else {
            --m_current;
        }

        return *this;
    }

    self_type operator + (difference_type addend) const
    {
        return self_type(m_block_pointers, m_block_items, m_index + addend);
    }
    self_type& operator += (difference_type addend)
    {
        m_index += addend;
        update();
        return *this;
    }
    self_type operator - (difference_type subtrahend) const
    {
        return self_type(m_block_pointers, m_block_items, m_index - subtrahend);
    }
    difference_type operator - (const self_type& o) const
    {
        return (m_index - o.m_index);
    }
    self_type& operator -= (difference_type subtrahend)
    {
        m_index -= subtrahend;
        update();
        return *this;
    }
    bool operator == (const self_type& o) const
    {
        return m_index == o.m_index;
    }
    bool operator != (const self_type& o) const
    {
        return m_index != o.m_index;
    }
    bool operator < (const self_type& o) const
    {
        return m_index < o.m_index;
    }
    bool operator <= (const self_type& o) const
    {
        return m_index <= o.m_index;
    }
    bool operator > (const self_type& o) const
    {
        return m_index > o.m_index;
    }
    bool operator >= (const self_type& o) const
    {
        return m_index >= o.m_index;
    }

    friend std::ostream& operator << (std::ostream& os, const ppq_iterator& i)
    {
        return os << "[" << i.m_index << "]";
    }

private:
    //! updates m_block_index and m_current based on m_index
    inline void update()
    {
        m_block_index = m_index / m_block_items;
        const size_t local_index = m_index % m_block_items;

        if (m_block_index < m_block_pointers->size()) {
            m_current = (*m_block_pointers)[m_block_index].first + local_index;
            assert(m_current <= (*m_block_pointers)[m_block_index].second);
        }
        else {
            // global end if end is beyond the last real block
            assert(m_block_index == m_block_pointers->size());
            assert(local_index == 0);
            //-tb old: m_current = (*m_block_pointers)[m_block_index - 1].second;
            m_current = NULL;
        }
    }
};

/*!
 * Internal arrays store a sorted sequence of values in RAM, which will be
 * merged together into the deletion buffer when it needs to be
 * refilled. Internal arrays are constructed from the insertions heaps when
 * they overflow.
 */
template <class ValueType>
class internal_array : private noncopyable
{
public:
    typedef ValueType value_type;
    typedef ppq_iterator<value_type> iterator;

protected:
    typedef typename iterator::block_pointers_type block_pointers_type;

    //! Contains the items of the sorted sequence.
    std::vector<value_type> m_values;

    //! Index of the current head
    unsigned_type m_min_index;

    //! Level of internal array (Sander's PQ: group number)
    unsigned_type m_level;

    //! Begin and end pointers of the array
    //! This is used by the iterator
    block_pointers_type m_block_pointers;

public:
    //! Default constructor. Don't use this directy. Needed for regrowing in
    //! surrounding vector.
    internal_array() : m_min_index(0) { }

    //! Constructor which takes a value vector. The value vector is empty
    //! afterwards.
    internal_array(std::vector<value_type>& values,
                   unsigned_type min_index = 0,
                   unsigned_type level = 0)
        : m_values(), m_min_index(min_index), m_level(level),
          m_block_pointers(1)
    {
        std::swap(m_values, values);
        STXXL_ASSERT(values.size() > 0);
        m_block_pointers[0] = std::make_pair(&(*m_values.begin()), &(*m_values.begin()) + m_values.size());
    }

    //! Swap internal_array with another one.
    void swap(internal_array& o)
    {
        using std::swap;

        swap(m_values, o.m_values);
        swap(m_min_index, o.m_min_index);
        swap(m_level, o.m_level);
        swap(m_block_pointers, o.m_block_pointers);
    }

    //! Swap internal_array with another one.
    friend void swap(internal_array& a, internal_array& b)
    {
        a.swap(b);
    }

    //! Random access operator
    inline value_type& operator [] (size_t i)
    {
        return m_values[i];
    }

    //! Use inc_min(diff) if multiple values have been extracted.
    inline void inc_min(size_t diff = 1)
    {
        m_min_index += diff;
    }

    //! The currently smallest element in the array.
    inline const value_type & get_min() const
    {
        return m_values[m_min_index];
    }

    //! The index of the currently smallest element in the array.
    inline size_t get_min_index() const
    {
        return m_min_index;
    }

    //! The index of the largest element in the array.
    inline size_t get_max_index() const
    {
        return (m_values.size() - 1);
    }

    //! Returns if the array has run empty.
    inline bool empty() const
    {
        return (m_min_index >= m_values.size());
    }

    //! Make this array empty.
    inline void make_empty()
    {
        m_min_index = m_values.size();
    }

    //! Returns the current size of the array.
    inline size_t size() const
    {
        return (m_values.size() - m_min_index);
    }

    //! Returns the initial size of the array.
    inline size_t capacity() const
    {
        return m_values.size();
    }

    //! Returns the level (group number) of the array.
    inline unsigned_type level() const
    {
        return m_level;
    }

    //! Return the amount of internal memory used by an array with the capacity
    //! in number of items.
    static size_t int_memory(size_t capacity)
    {
        return sizeof(internal_array) + capacity * sizeof(value_type);
    }

    //! Return the amount of internal memory used by the array
    inline size_t int_memory() const
    {
        return int_memory(m_values.capacity());
    }

    //! Begin iterator
    inline iterator begin() const
    {
        // not const, unfortunately.
        return iterator(&m_block_pointers, capacity(), m_min_index);
    }

    //! End iterator
    inline iterator end() const
    {
        // not const, unfortunately.
        return iterator(&m_block_pointers, capacity(), capacity());
    }
};

template <class ExternalArrayType>
class external_array_writer;

/*!
 * External array stores a sorted sequence of values on the hard disk and
 * allows access to the first block (containing the smallest values).  The
 * class uses buffering and prefetching in order to improve the performance.
 *
 * \tparam ValueType Type of the contained objects (POD with no references to
 * internal memory).
 *
 * \tparam BlockSize External block size. Default =
 * STXXL_DEFAULT_BLOCK_SIZE(ValueType).
 *
 * \tparam AllocStrategy Allocation strategy for the external memory. Default =
 * STXXL_DEFAULT_ALLOC_STRATEGY.
 */
template <
    class ValueType,
    unsigned_type BlockSize = STXXL_DEFAULT_BLOCK_SIZE(ValueType),
    class AllocStrategy = STXXL_DEFAULT_ALLOC_STRATEGY
    >
class external_array : private noncopyable
{
public:
    typedef ValueType value_type;
    typedef ppq_iterator<value_type> iterator;

    typedef external_array<value_type, BlockSize, AllocStrategy> self_type;
    typedef typed_block<BlockSize, value_type> block_type;
    typedef read_write_pool<block_type> pool_type;
    typedef std::vector<BID<BlockSize> > bid_vector;
    typedef typename bid_vector::iterator bid_iterator;
    typedef std::vector<block_type*> block_vector;
    typedef std::vector<request_ptr> request_vector;
    typedef std::vector<value_type> minima_vector;
    typedef typename iterator::block_pointers_type block_pointers_type;
    typedef external_array_writer<self_type> writer_type;

    //! The number of elements fitting into one block
    enum {
        block_size = BlockSize,
        block_items = BlockSize / sizeof(value_type)
    };

    static const bool debug = false;

protected:
    //! The total size of the external array in items. Cannot be changed
    //! after construction.
    external_size_type m_capacity;

    //! Number of blocks, again: calculated at construction time.
    unsigned_type m_num_blocks;

    //! Level of external array (Sander's PQ: group number)
    unsigned_type m_level;

    //! Common prefetch and write buffer pool
    pool_type* m_pool;

    //! The IDs of each block in external memory.
    bid_vector m_bids;

    //! A vector of size m_num_blocks with block_type pointers, some of them
    //! will be filled while writing, but most are NULL.
    block_vector m_blocks;

    //! Begin and end pointers for each block, used for merging with
    //! ppq_iterator.
    block_pointers_type m_block_pointers;

    //! The read request pointers are used to wait until the block has been
    //! completely fetched.
    request_vector m_requests;

    //! stores the minimum value of each block
    minima_vector m_minima;

    //! Is array in write phase? True = write phase, false = read phase.
    bool m_write_phase;

    //! The total number of elements minus the number of extracted values
    external_size_type m_size;

    //! The read position in the array.
    external_size_type m_index;

    //! The index behind the last element that is located in RAM (or is at
    //! least requested to be so)
    external_size_type m_end_index;

    //! The first unhinted block index.
    unsigned_type m_unhinted_block;

    //! The first unhinted block index as it was before the
    //! prepare_rebuilding_hints() call. Used for removal of hints which aren't
    //! needed anymore.
    unsigned_type m_old_unhinted_block;

    //! allow writer to access to all variables
    friend class external_array_writer<self_type>;

public:
    /*!
     * Constructs an external array
     *
     * \param size The total number of elements. Cannot be changed after
     * construction.
     *
     * \param num_prefetch_blocks Number of blocks to prefetch from hard disk
     *
     * \param num_write_buffer_blocks Size of the write buffer in number of
     * blocks
     */
    external_array(external_size_type size, pool_type* pool, unsigned_type level = 0)
        :   // constants
          m_capacity(size),
          m_num_blocks((size_t)div_ceil(m_capacity, block_items)),
          m_level(level),
          m_pool(pool),

          // vectors
          m_bids(m_num_blocks),
          m_blocks(m_num_blocks, reinterpret_cast<block_type*>(1)),
          m_block_pointers(m_num_blocks),
          m_requests(m_num_blocks, NULL),
          m_minima(m_num_blocks),

          // state
          m_write_phase(true),

          // indices
          m_size(0),
          m_index(0),
          m_end_index(0),
          m_unhinted_block(0),
          m_old_unhinted_block(0)
    {
        assert(m_capacity > 0);
        // allocate blocks in EM.
        block_manager* bm = block_manager::get_instance();
        bm->new_blocks(AllocStrategy(), m_bids.begin(), m_bids.end());
    }

    //! Default constructor. Don't use this directy. Needed for regrowing in
    //! surrounding vector.
    external_array()
        :   // constants
          m_capacity(0),
          m_num_blocks(0),
          m_level(0),
          m_pool(NULL),

          // vectors
          m_bids(0),
          m_blocks(0),
          m_block_pointers(0),
          m_requests(0),
          m_minima(0),

          // state
          m_write_phase(false),

          // indices
          m_size(0),
          m_index(0),
          m_end_index(0),
          m_unhinted_block(0),
          m_old_unhinted_block(0)
    { }

    //! Swap external_array with another one.
    void swap(external_array& o)
    {
        using std::swap;

        // constants
        swap(m_capacity, o.m_capacity);
        swap(m_num_blocks, o.m_num_blocks);
        swap(m_level, o.m_level);
        swap(m_pool, o.m_pool);

        // vectors
        swap(m_bids, o.m_bids);
        swap(m_requests, o.m_requests);
        swap(m_blocks, o.m_blocks);
        swap(m_block_pointers, o.m_block_pointers);
        swap(m_minima, o.m_minima);

        // state
        swap(m_write_phase, o.m_write_phase);

        // indices
        swap(m_size, o.m_size);
        swap(m_index, o.m_index);
        swap(m_end_index, o.m_end_index);
        swap(m_unhinted_block, o.m_unhinted_block);
        swap(m_old_unhinted_block, o.m_old_unhinted_block);
    }

    //! Swap external_array with another one.
    friend void swap(external_array& a, external_array& b)
    {
        a.swap(b);
    }

    //! Destructor
    ~external_array()
    {
        if (m_size == 0) return;

        // not all data has been read! this only happen when the PPQ is
        // destroyed while containing data.

        const unsigned_type block_index = m_index / block_items;
        const unsigned_type end_block_index = get_end_block_index();

        // released blocks currently held in RAM
        for (size_t i = block_index; i < end_block_index; ++i) {
            m_pool->add_prefetch(m_blocks[i]);
            // cannot report the number of freed blocks to PPQ.
        }

        // cancel currently hinted blocks
        for (size_t i = end_block_index; i < m_unhinted_block; ++i) {
            STXXL_DEBUG("ea[" << this << "]: discarding prefetch hint on"
                        " block " << i);

            m_requests[i]->cancel();
            m_requests[i]->wait();
            // put block back into pool
            m_pool->add_prefetch(m_blocks[i]);
            // invalidate block entry
            m_blocks[i] = NULL;
            m_requests[i] = request_ptr();
        }

        // figure out first block that is still allocated in EM.
        bid_iterator i_begin = m_bids.begin() + block_index;
        block_manager::get_instance()->delete_blocks(i_begin, m_bids.end());

        // check that all is empty
        for (size_t i = block_index; i < end_block_index; ++i)
            assert(m_blocks[i] == NULL);
    }

    //! Returns the capacity in items.
    size_t capacity() const
    {
        return m_capacity;
    }

    //! Returns the current size in items.
    size_t size() const
    {
        return m_size;
    }

    //! Returns true if the array is empty.
    bool empty() const
    {
        return (m_size == 0);
    }

    //! Returns the level (group number) of the array.
    inline unsigned_type level() const
    {
        return m_level;
    }

    //! Return the number of blocks.
    size_t num_blocks() const
    {
        return m_num_blocks;
    }

    //! Returns memory usage of EA with given capacity, excluding blocks loaded
    //! in RAM. Blocks belong to prefetch pool.
    static size_t int_memory(size_t capacity)
    {
        size_t num_blocks = div_ceil(capacity, block_items);

        return sizeof(external_array)
               + num_blocks * sizeof(typename bid_vector::value_type)
               + num_blocks * sizeof(typename block_vector::value_type)
               + num_blocks * sizeof(typename block_pointers_type::value_type)
               + num_blocks * sizeof(typename request_vector::value_type)
               + num_blocks * sizeof(typename minima_vector::value_type);
    }

    //! Return the amount of internal memory used by the EA.
    inline size_t int_memory() const
    {
        return int_memory(m_capacity);
    }

    //! Returns the number elements available in internal memory
    size_t buffer_size() const
    {
        return (m_end_index - m_index);
    }

    //! Returns the block beyond the block in which *(m_end_index-1) is located.
    unsigned_type get_end_block_index() const
    {
        unsigned_type end_block_index = m_end_index / block_items;

        // increase block index if inside the block
        if (m_end_index % block_items != 0) ++end_block_index;
        assert(end_block_index <= m_num_blocks);

        return end_block_index;
    }

    //! Returns the block in which m_index is located.
    inline unsigned_type get_current_block_index() const
    {
        return (m_index / block_items);
    }

    //! Returns a random-access iterator to the begin of the data
    //! in internal memory.
    iterator begin() const
    {
        //-TODO?: assert(block_valid(m_index / block_items) || m_index == m_capacity);
        return iterator(&m_block_pointers, block_items, m_index);
    }

    //! Returns a random-access iterator 1 behind the end of the data
    //! in internal memory.
    iterator end() const
    {
        //-TODO? assert(!block_valid(m_end_index / block_items) || m_end_index == m_capacity);
        return iterator(&m_block_pointers, block_items, m_end_index);
    }

    //! Returns the smallest element in the array
    const value_type & get_min()
    {
        return *begin();
    }

    //! Returns if there is data in EM, that's not randomly accessible.
    bool has_em_data() const
    {
        return (get_end_block_index() < m_num_blocks);
    }

    //! Returns the smallest element of the first block NOT in internal memory
    //! (or at least requested to be in internal memory)
    const value_type & get_next_block_min() const
    {
        assert(get_end_block_index() < m_num_blocks);
        return m_minima[get_end_block_index()];
    }

    //! Returns if the data requested to be in internal memory is
    //! completely fetched. True if wait() has been called before.
    bool valid() const
    {
        bool result = true;
        const unsigned_type block_index = m_index / block_items;
        const unsigned_type end_block_index = get_end_block_index();
        for (unsigned_type i = block_index; i < end_block_index; ++i) {
            result = result && block_valid(i);
        }
        return result;
    }

    //! Random access operator for data in internal memory
    //! You should call wait() once after fetching data from EM.
    value_type& operator [] (size_t i) const
    {
        assert(i < m_capacity);
        const size_t block_index = i / block_items;
        const size_t local_index = i % block_items;
        assert(i < m_capacity);
        assert(block_valid(block_index));
        return m_blocks[block_index]->elem[local_index];
    }

public:
    //! prepare the pool for writing external arrays with given number of
    //! threads
    static void prepare_write_pool(pool_type& pool, unsigned_type num_threads)
    {
        unsigned_type write_blocks = num_threads;
        // need at least one
        if (write_blocks == 0) write_blocks = 1;
        // for holding boundary blocks
        write_blocks *= 2;
        // more disks than threads?
        if (write_blocks < config::get_instance()->disks_number())
            write_blocks = config::get_instance()->disks_number();
#if STXXL_DEBUG_ASSERTIONS
        // required for re-reading the external array
        write_blocks = 2 * write_blocks;
#endif
        if (pool.size_write() < write_blocks) {
            STXXL_ERRMSG("WARNING: enlarging PPQ write pool to " <<
                         write_blocks << " blocks = " <<
                         write_blocks * block_size / 1024 / 1024 << " MiB");
            pool.resize_write(write_blocks);
        }
    }

protected:
    //! prepare the external_array for writing using multiway_merge() with
    //! num_threads. this method is called by the external_array_writer's
    //! constructor.
    void prepare_write(unsigned_type num_threads)
    {
        prepare_write_pool(*m_pool, num_threads);
    }

    //! finish the writing phase after multiway_merge() filled the vector. this
    //! method is called by the external_array_writer's destructor..
    void finish_write()
    {
        // check that all blocks where written
        for (unsigned_type i = 0; i < m_num_blocks; ++i)
            assert(m_blocks[i] == NULL);

        // compatibility to the block write interface
        m_size = m_capacity;
        m_index = 0;
        m_end_index = 0;
        m_unhinted_block = 0;

        m_write_phase = false;
    }

    //! Called by the external_array_writer to read a block from disk into
    //! m_blocks[]. If the block is marked as uninitialized, then no read is
    //! performed. This is the usual case, and in theory, no block ever has be
    //! re-read from disk, since all can be written fully. However, we do
    //! support re-reading blocks for debugging purposes inside
    //! multiway_merge(), in a full performance build re-reading never occurs.
    void read_block(size_t block_index)
    {
        assert(block_index < m_num_blocks);
        assert(m_blocks[block_index] == NULL ||
               m_blocks[block_index] == reinterpret_cast<block_type*>(1));

        if (m_blocks[block_index] == reinterpret_cast<block_type*>(1))
        {
            // special marker: this block is uninitialized -> no need to read
            // from disk.
            m_blocks[block_index] = m_pool->steal();
        }
        else
        {
            // block was already written, have to read from EM.
            STXXL_DEBUG("ea[" << this << "]: "
                        "read_block needs to re-read block index=" << block_index);

            static bool s_warned = false;
            if (!s_warned)
            {
                s_warned = true;
                STXXL_ERRMSG("ppq::external_array[" << this << "] "
                             "writer requested to re-read block from EM.");
                STXXL_ERRMSG("This should never occur in full-performance mode, "
                             "verify that you run in debug mode.");
            }

            // this re-reading is not necessary for full performance builds, so
            // we immediately wait for the I/O to be completed.
            m_blocks[block_index] = m_pool->steal();
            request_ptr req = m_pool->read(m_blocks[block_index], m_bids[block_index]);
            req->wait();
            assert(req->poll());
            assert(m_blocks[block_index]);
        }
    }

    //! Called by the external_array_writer to write a block from m_blocks[] to
    //! disk. Prior to writing and releasing the memory, extra information is
    //! preserved.
    void write_block(size_t block_index)
    {
        assert(block_index < m_num_blocks);
        assert(m_blocks[block_index] != NULL &&
               m_blocks[block_index] != reinterpret_cast<block_type*>(1));

        // calculate minimum and maximum values
        const internal_size_type this_block_items =
            std::min<internal_size_type>(block_items, m_capacity - block_index * (external_size_type)block_items);

        STXXL_DEBUG("ea[" << this << "]: write_block index=" << block_index <<
                    " this_block_items=" << this_block_items);

        assert(this_block_items > 0);
        block_type& this_block = *m_blocks[block_index];

        m_minima[block_index] = this_block[0];

        // write out block (in background)
        m_pool->write(m_blocks[block_index], m_bids[block_index]);

        m_blocks[block_index] = NULL;
    }

public:
    //! \name Prefetching Hints
    //! \{

    //! Prefetch the next unhinted block, requires one free read block from the
    //! global pool.
    void hint_next_block()
    {
        assert(m_unhinted_block < m_num_blocks);

        // will read (prefetch) block i
        size_t i = m_unhinted_block++;

        STXXL_DEBUG("ea[" << this << "]: prefetching block_index=" << i);

        assert(m_pool->size_write() > 0);
        assert(m_blocks[i] == NULL);

        // steal block from pool, but also perform read via pool, since this
        // checks the associated write_pool.
        m_blocks[i] = m_pool->steal_prefetch();
        m_requests[i] = m_pool->read(m_blocks[i], m_bids[i]);
    }

    //! Returns if there is data in EM, that's not already hinted
    //! to the prefetcher.
    bool has_unhinted_em_data() const
    {
        return (m_unhinted_block < m_num_blocks);
    }

    //! Returns the smallest element of the next hint candidate (the block
    //! after the last hinted one).
    const value_type & get_next_hintable_min() const
    {
        assert(m_unhinted_block < m_num_blocks);
        return m_minima[m_unhinted_block];
    }

    //! Returns the number of hinted blocks.
    size_t num_hinted_blocks() const
    {
        assert(get_end_block_index() <= m_unhinted_block);
        return m_unhinted_block - get_end_block_index();
    }

    //! This method prepares rebuilding the hints (this is done after creating
    //! a new EA in order to always have globally the n blocks hinted which
    //! will be fetched first).  Resets m_unhinted_block to the first block not
    //! in RAM. Thereafter prehint_next_block() is used to advance this index.
    //! finish_rebuilding_hints() should be called after placing all hints in
    //! order to clean up the prefetch pool.
    void rebuild_hints_prepare()
    {
        m_old_unhinted_block = m_unhinted_block;
        m_unhinted_block = get_end_block_index();
        assert(get_end_block_index() <= m_old_unhinted_block);
    }

    //! Advance m_unhinted_block index without actually prefetching.
    void rebuild_hints_prehint_next_block()
    {
        assert(m_unhinted_block < m_num_blocks);

        // will read (prefetch) block after cancellations.

        STXXL_DEBUG("ea[" << this << "]: pre-hint of" <<
                    " block_index=" << m_unhinted_block);

        ++m_unhinted_block;
    }

    //! Cancel hints which aren't needed anymore from the prefetcher and fixes
    //! it's size. prepare_rebuilding_hints() must be called before!
    void rebuild_hints_cancel()
    {
        for (size_t i = m_unhinted_block; i < m_old_unhinted_block; ++i) {
            STXXL_DEBUG("ea[" << this << "]: discarding prefetch hint on"
                        " block " << i);
            m_requests[i]->cancel();
            m_requests[i]->wait();
            // put block back into pool
            m_pool->add_prefetch(m_blocks[i]);
            // invalidate block entry
            m_blocks[i] = NULL;
            m_requests[i] = request_ptr();
        }
    }

    //! Perform real-hinting of pre-hinted blocks, since now canceled blocks
    //! are available.
    void rebuild_hints_finish()
    {
        for (size_t i = m_old_unhinted_block; i < m_unhinted_block; ++i)
        {
            STXXL_DEBUG("ea[" << this << "]: perform real-hinting of"
                        " block " << i);

            assert(m_pool->size_write() > 0);
            assert(m_blocks[i] == NULL);
            m_blocks[i] = m_pool->steal_prefetch();
            m_requests[i] = m_pool->read(m_blocks[i], m_bids[i]);
        }
    }

    //! \}

public:
    //! \name Waiting and Removal
    //! \{

    //! Waits until the next prefetched block is read into RAM, then polls for
    //! any further blocks that are done as well. Returns how many blocks were
    //! successfully read.
    unsigned_type wait_next_blocks()
    {
        size_t begin = get_end_block_index(), i = begin;

        STXXL_DEBUG("ea[" << this << "]: waiting for" <<
                    " block index=" << i <<
                    " end_index=" << m_end_index);

        assert(has_em_data());

        assert(i < m_unhinted_block);
        assert(m_bids[i].valid());
        assert(m_requests[i].valid());

        // wait for prefetched request to finish.
        m_requests[i]->wait();
        assert(m_requests[i]->poll());
        assert(m_blocks[i]);

        update_block_pointers(i);
        ++i;

        // poll further hinted blocks if already done
        while (i < m_unhinted_block && m_requests[i]->poll())
        {
            STXXL_DEBUG("ea[" << this << "]: poll-ok for" <<
                        " block index=" << i <<
                        " end_index=" << m_end_index);
            m_requests[i]->wait();
            assert(m_requests[i]->poll());
            assert(m_blocks[i]);

            update_block_pointers(i);
            ++i;
        }

        m_end_index = std::min(m_capacity, i * (external_size_type)block_items);

        return i - begin;
    }

    //! Waits until all hinted blocks are read into RAM. Returns how many
    //! blocks were successfully read.
    unsigned_type wait_all_hinted_blocks()
    {
        size_t begin = get_end_block_index(), i = begin;
        while (i < m_unhinted_block)
        {
            STXXL_DEBUG("wait_all_hinted_blocks(): ea[" << this << "]: waiting for" <<
                        " block index=" << i <<
                        " end_index=" << m_end_index);
            m_requests[i]->wait();
            assert(m_requests[i]->poll());
            assert(m_blocks[i]);
            update_block_pointers(i);
            ++i;
        }
        m_end_index = std::min(m_capacity, i * (external_size_type)block_items);
        return i - begin;
    }

    //! Returns the number of blocks loaded in RAM.
    size_t num_used_blocks() const
    {
        return get_end_block_index() - (m_index / block_items);
    }

    //! Removes the first n elements from the array. Returns the number of
    //! blocks released into the block pool.
    unsigned_type remove_items(size_t n)
    {
        assert(m_index + n <= m_capacity);
        assert(m_index + n <= m_end_index);
        assert(m_size >= n);

        STXXL_DEBUG("ea[" << this << "]: remove " << n << " items");

        if (n == 0)
            return 0;

        const size_t block_index = m_index / block_items;

        const size_t index_after = m_index + n;
        size_t block_index_after = index_after / block_items;
        size_t local_index_after = index_after % block_items;

        if (m_size == n && local_index_after != 0) // end of EA
            ++block_index_after;

        assert(block_index_after <= m_num_blocks);

        bid_iterator i_begin = m_bids.begin() + block_index;
        bid_iterator i_end = m_bids.begin() + block_index_after;
        assert(i_begin <= i_end);
        block_manager::get_instance()->delete_blocks(i_begin, i_end);

        for (size_t i = block_index; i < block_index_after; ++i) {
            assert(block_valid(i));
            // return block to pool
            m_pool->add_prefetch(m_blocks[i]);
        }

        m_index = index_after;
        m_size -= n;

        unsigned_type blocks_freed = block_index_after - block_index;

        STXXL_DEBUG("ea[" << this << "]: after remove:" <<
                    " index_after=" << index_after <<
                    " block_index_after=" << block_index_after <<
                    " local_index_after=" << local_index_after <<
                    " blocks_freed=" << blocks_freed <<
                    " num_blocks=" << m_num_blocks <<
                    " capacity=" << m_capacity);

        assert(block_index_after <= m_num_blocks);
        // at most one block outside of the currently loaded range
        assert(block_index_after <= get_end_block_index());

        return blocks_freed;
    }

    //! \}

protected:
    //! Returns if the block with the given index is completely fetched.
    bool block_valid(size_t block_index) const
    {
        if (!m_write_phase) {
            if (block_index >= m_num_blocks) return false;
            return (m_requests[block_index] && m_requests[block_index]->poll());
        }
        else {
            return (bool)m_blocks[block_index];
        }
    }

    //! Updates the m_block_pointers vector.
    //! Should be called after any steal() or read() operation.
    //! This is necessary for the iterators to work properly.
    inline void update_block_pointers(size_t block_index)
    {
        STXXL_DEBUG("ea[" << this << "]: updating block pointers for " << block_index);

        m_block_pointers[block_index].first = m_blocks[block_index]->begin();
        if (block_index + 1 != m_num_blocks)
            m_block_pointers[block_index].second = m_blocks[block_index]->end();
        else
            m_block_pointers[block_index].second =
                m_block_pointers[block_index].first
                + (m_capacity - block_index * block_items);

        assert(m_block_pointers[block_index].first != NULL);
        assert(m_block_pointers[block_index].second != NULL);
    }

    inline size_t last_block_items()
    {
        size_t mod = m_capacity % block_items;
        return (mod > 0) ? mod : (size_t)block_items;
    }
};

/**
 * An external_array can only be written using an external_array_writer
 * object. The writer objects provides iterators which are designed to be used
 * by stxxl::parallel::multiway_merge() to write the external memory blocks in
 * parallel. Thus in the writer we coordinate thread-safe access to the blocks
 * using reference counting.
 *
 * An external_array_writer::iterator has two states: normal and "live". In
 * normal mode, the iterator only has a valid index into the external array's
 * items. In normal mode, only index calculations are possible. Once
 * operator*() is called, the iterators goes into "live" mode by requesting
 * access to the corresponding block. Using reference counting the blocks is
 * written once all iterators are finished with the corresponding block. Since
 * with operator*() we cannot know if the value is going to be written or read,
 * when going to live mode, the block must be read from EM. This read overhead,
 * however, is optimized by marking blocks as uninitialized in external_array,
 * and skipping reads for then. In a full performance build, no block needs to
 * be read from disk. Reads only occur in debug mode, when the results are
 * verify.
 *
 * The iterator's normal/live mode only stays active for the individual
 * iterator object. When an iterator is copied/assigned/calculated with the
 * mode is NOT inherited! The exception is prefix operator ++, which is used by
 * multiway_merge() to fill an array. Thus the implementation of the iterator
 * heavily depends on the behavior of multiway_merge() and is optimized for it.
 */
template <class ExternalArrayType>
class external_array_writer : public noncopyable
{
public:
    typedef ExternalArrayType ea_type;

    typedef external_array_writer self_type;

    typedef typename ea_type::value_type value_type;
    typedef typename ea_type::block_type block_type;

    //! prototype declaration of nested class.
    class iterator;

    //! scope based debug variable
    static const bool debug = false;

protected:
    //! reference to the external array to be written
    ea_type& m_ea;

#ifndef NDEBUG
    //! total number of iterators referencing this writer
    unsigned int m_ref_total;
#endif

    //! reference counters for the number of live iterators on the
    //! corresponding block in external_array.
    std::vector<unsigned int> m_ref_count;

    //! mutex for reference counting array (this is actually nicer than
    //! openmp's critical)
    mutex m_mutex;

    //! optimization: hold live iterators for the expected boundary blocks of
    //! multiway_merge().
    std::vector<iterator> m_live_boundary;

protected:
    //! read block into memory and increase reference count (called when an
    //! iterator goes live on the block).
    block_type * get_block_ref(size_t block_index)
    {
        scoped_mutex_lock lock(m_mutex);

        assert(block_index < m_ea.num_blocks());
        unsigned int ref = m_ref_count[block_index]++;
#ifndef NDEBUG
        ++m_ref_total;
#endif

        if (ref == 0) {
            STXXL_DEBUG("get_block_ref block_index=" << block_index <<
                        " ref=" << ref << " reading.");
            m_ea.read_block(block_index);
        }
        else {
            STXXL_DEBUG("get_block_ref block_index=" << block_index <<
                        " ref=" << ref);
        }

        return m_ea.m_blocks[block_index];
    }

    //! decrease reference count on the block, and possibly write it to disk
    //! (called when an iterator releases live mode).
    void free_block_ref(size_t block_index)
    {
        scoped_mutex_lock lock(m_mutex);

        assert(block_index < m_ea.num_blocks());
#ifndef NDEBUG
        assert(m_ref_total > 0);
        --m_ref_total;
#endif
        unsigned int ref = --m_ref_count[block_index];

        if (ref == 0) {
            STXXL_DEBUG("free_block_ref block_index=" << block_index <<
                        " ref=" << ref << " written.");
            m_ea.write_block(block_index);
        }
        else {
            STXXL_DEBUG("free_block_ref block_index=" << block_index <<
                        " ref=" << ref);
        }
    }

    //! allow access to the block_ref functions
    friend class iterator;

public:
    /**
     * An iterator which can be used to write (and read) an external_array via
     * an external_array_writer. See the documentation of external_array_writer.
     */
    class iterator
    {
    public:
        typedef external_array_writer writer_type;
        typedef ExternalArrayType ea_type;

        typedef typename ea_type::value_type value_type;
        typedef value_type& reference;
        typedef value_type* pointer;
        typedef ptrdiff_t difference_type;
        typedef std::random_access_iterator_tag iterator_category;

        typedef iterator self_type;

        static const size_t block_items = ea_type::block_items;

        //! scope based debug variable
        static const bool debug = false;

    protected:
        //! pointer to the external array containing the elements
        writer_type* m_writer;

        //! when operator* or operator-> are called, then the iterator goes
        //! live and allocates a reference to the block's data (possibly
        //! reading it from EM).
        bool m_live;

        //! index of the current element, absolute in the external array
        external_size_type m_index;

        //! index of the current element's block in the external array's block
        //! list. undefined while m_live is false.
        internal_size_type m_block_index;

        //! pointer to the referenced block. undefined while m_live is false.
        block_type* m_block;

        //! pointer to the current element inside the referenced block.
        //! undefined while m_live is false.
        internal_size_type m_current;

    public:
        //! default constructor (should not be used directly)
        iterator()
            : m_writer(NULL), m_live(false), m_index(0)
        { }

        //! construct a new iterator
        iterator(writer_type* writer, external_size_type index)
            : m_writer(writer),
              m_live(false),
              m_index(index)
        {
            STXXL_DEBUG("Construct iterator for index " << m_index);
        }

        //! copy an iterator, the new iterator is _not_ automatically live!
        iterator(const iterator& other)
            : m_writer(other.m_writer),
              m_live(false),
              m_index(other.m_index)
        {
            STXXL_DEBUG("Copy-Construct iterator for index " << m_index);
        }

        //! assign an iterator, the assigned iterator is not automatically live!
        iterator& operator = (const iterator& other)
        {
            if (&other != this)
            {
                STXXL_DEBUG("Assign iterator to index " << other.m_index);

                if (m_live)
                    m_writer->free_block_ref(m_block_index);

                m_writer = other.m_writer;
                m_live = false;
                m_index = other.m_index;
            }

            return *this;
        }

        ~iterator()
        {
            if (!m_live) return; // no need for cleanup

            m_writer->free_block_ref(m_block_index);

            STXXL_DEBUG("Destruction of iterator for index " << m_index <<
                        " in block " << m_index / block_items);
        }

        //! return the current absolute index inside the external array.
        external_size_type get_index() const
        {
            return m_index;
        }

        //! allocates a reference to the block's data (possibly reading it from
        //! EM).
        void make_live()
        {
            assert(!m_live);

            // calculate block and index inside
            m_block_index = m_index / block_items;
            m_current = m_index % block_items;

            STXXL_DEBUG("operator*() live request for index=" << m_index <<
                        " block_index=" << m_block_index <<
                        " m_current=" << m_current);

            // get block reference
            m_block = m_writer->get_block_ref(m_block_index);
            m_live = true;
        }

        //! access the current item
        reference operator * ()
        {
            if (UNLIKELY(!m_live))
                make_live();

            return (*m_block)[m_current];
        }

        //! access the current item
        pointer operator -> ()
        {
            return &(operator * ());
        }

        //! prefix-increment operator
        self_type& operator ++ ()
        {
            ++m_index;
            if (UNLIKELY(!m_live)) return *this;

            // if index stays in the same block, everything is fine
            ++m_current;
            if (LIKELY(m_current != block_items)) return *this;

            // release current block
            m_writer->free_block_ref(m_block_index);
            m_live = false;

            return *this;
        }

        self_type operator + (difference_type addend) const
        {
            return self_type(m_writer, m_index + addend);
        }
        self_type operator - (difference_type subtrahend) const
        {
            return self_type(m_writer, m_index - subtrahend);
        }
        difference_type operator - (const self_type& o) const
        {
            return (m_index - o.m_index);
        }

        bool operator == (const self_type& o) const
        {
            return m_index == o.m_index;
        }
        bool operator != (const self_type& o) const
        {
            return m_index != o.m_index;
        }
        bool operator < (const self_type& o) const
        {
            return m_index < o.m_index;
        }
        bool operator <= (const self_type& o) const
        {
            return m_index <= o.m_index;
        }
        bool operator > (const self_type& o) const
        {
            return m_index > o.m_index;
        }
        bool operator >= (const self_type& o) const
        {
            return m_index >= o.m_index;
        }
    };

public:
    external_array_writer(ea_type& ea, unsigned int num_threads = 0)
        : m_ea(ea),
          m_ref_count(ea.num_blocks(), 0)
    {
#ifndef NDEBUG
        m_ref_total = 0;
#endif

#if STXXL_PARALLEL
        if (num_threads == 0)
            num_threads = omp_get_max_threads();
#else
        if (num_threads == 0)
            num_threads = 1;
#endif

        m_ea.prepare_write(num_threads);

        // optimization: hold live iterators for the boundary blocks which two
        // threads write to. this prohibits the blocks to be written to disk
        // and read again.

        double step = (double)m_ea.capacity() / (double)num_threads;
        m_live_boundary.resize(num_threads - 1);

        for (unsigned int i = 0; i < num_threads - 1; ++i)
        {
            external_size_type index = (external_size_type)((i + 1) * step);
            STXXL_DEBUG("hold index " << index <<
                        " in block " << index / ea_type::block_items);
            m_live_boundary[i] = iterator(this, index);
            m_live_boundary[i].make_live();
        }
    }

    ~external_array_writer()
    {
        m_live_boundary.clear(); // release block boundaries
#ifndef NDEBUG
        STXXL_ASSERT(m_ref_total == 0);
#endif
        m_ea.finish_write();
    }

    iterator begin()
    {
        return iterator(this, 0);
    }

    iterator end()
    {
        return iterator(this, m_ea.capacity());
    }
};

/*!
 * The minima_tree contains minima from all sources inside the PPQ. It contains
 * four substructures: winner trees for insertion heaps, internal and external
 * arrays, each containing the minima from all currently allocated
 * structures. These three sources, plus the deletion buffer are combined using
 * a "head" inner tree containing only up to four item.
 */
template <class ParentType>
class minima_tree
{
public:
    typedef ParentType parent_type;
    typedef minima_tree<ParentType> self_type;

    typedef typename parent_type::inv_compare_type compare_type;
    typedef typename parent_type::value_type value_type;
    typedef typename parent_type::proc_vector_type proc_vector_type;
    typedef typename parent_type::internal_arrays_type ias_type;
    typedef typename parent_type::external_arrays_type eas_type;

    static const unsigned initial_ia_size = 2;
    static const unsigned initial_ea_size = 2;

protected:
    //! WinnerTree-Comparator for the head winner tree. It accesses all
    //! relevant data structures from the priority queue.
    struct head_comp
    {
        self_type& m_parent;
        proc_vector_type& m_proc;
        ias_type& m_ias;
        const compare_type& m_compare;

        head_comp(self_type& parent, proc_vector_type& proc,
                  ias_type& ias, const compare_type& compare)
            : m_parent(parent),
              m_proc(proc),
              m_ias(ias),
              m_compare(compare)
        { }

        const value_type & get_value(int input) const
        {
            switch (input) {
            case HEAP:
                return m_proc[m_parent.m_heaps.top()]->insertion_heap[0];
            case IA:
                return m_ias[m_parent.m_ia.top()].get_min();
            case EB:
                return m_parent.m_parent.m_extract_buffer[
                    m_parent.m_parent.m_extract_buffer_index
                ];
            default:
                abort();
            }
        }

        bool operator () (const int a, const int b) const
        {
            return m_compare(get_value(a), get_value(b));
        }
    };

    //! Comparator for the insertion heaps winner tree.
    struct heaps_comp
    {
        proc_vector_type& m_proc;
        const compare_type& m_compare;

        heaps_comp(proc_vector_type& proc, const compare_type& compare)
            : m_proc(proc), m_compare(compare)
        { }

        const value_type & get_value(int index) const
        {
            return m_proc[index]->insertion_heap[0];
        }

        bool operator () (const int a, const int b) const
        {
            return m_compare(get_value(a), get_value(b));
        }
    };

    //! Comparator for the internal arrays winner tree.
    struct ia_comp
    {
        ias_type& m_ias;
        const compare_type& m_compare;

        ia_comp(ias_type& ias, const compare_type& compare)
            : m_ias(ias), m_compare(compare)
        { }

        bool operator () (const int a, const int b) const
        {
            return m_compare(m_ias[a].get_min(), m_ias[b].get_min());
        }
    };

protected:
    //! The priority queue
    parent_type& m_parent;

    //! value_type comparator
    const compare_type& m_compare;

    //! Comperator instances
    head_comp m_head_comp;
    heaps_comp m_heaps_comp;
    ia_comp m_ia_comp;

    //! The winner trees
    winner_tree<head_comp> m_head;
    winner_tree<heaps_comp> m_heaps;
    winner_tree<ia_comp> m_ia;

public:
    //! Entries in the head winner tree.
    enum Types {
        HEAP = 0,
        IA = 1,
        EB = 2,
        TYPE_ERROR = 3
    };

    //! Construct the tree of minima sources.
    minima_tree(parent_type& parent)
        : m_parent(parent),
          m_compare(parent.m_inv_compare),
          // construct comparators
          m_head_comp(*this, parent.m_proc,
                      parent.m_internal_arrays, m_compare),
          m_heaps_comp(parent.m_proc, m_compare),
          m_ia_comp(parent.m_internal_arrays, m_compare),
          // construct header winner tree
          m_head(3, m_head_comp),
          m_heaps(m_parent.m_num_insertion_heaps, m_heaps_comp),
          m_ia(initial_ia_size, m_ia_comp)
    { }

    //! Return smallest items of head winner tree.
    std::pair<unsigned, unsigned> top()
    {
        unsigned type = m_head.top();
        switch (type)
        {
        case HEAP:
            return std::make_pair(HEAP, m_heaps.top());
        case IA:
            return std::make_pair(IA, m_ia.top());
        case EB:
            return std::make_pair(EB, 0);
        default:
            return std::make_pair(TYPE_ERROR, 0);
        }
    }

    //! Update minima tree after an item from the heap index was removed.
    void update_heap(int_type index)
    {
        m_heaps.notify_change(index);
        m_head.notify_change(HEAP);
    }

    //! Update minima tree after an item of the extract buffer was removed.
    void update_extract_buffer()
    {
        m_head.notify_change(EB);
    }

    //! Update minima tree after an item from an internal array was removed.
    void update_internal_array(unsigned index)
    {
        m_ia.notify_change(index);
        m_head.notify_change(IA);
    }

    //! Add a newly created internal array to the minima tree.
    void add_internal_array(unsigned index)
    {
        m_ia.activate_player(index);
        m_head.notify_change(IA);
    }

    //! Remove an insertion heap from the minima tree.
    void deactivate_heap(unsigned index)
    {
        m_heaps.deactivate_player(index);
        if (!m_heaps.empty())
            m_head.notify_change(HEAP);
        else
            m_head.deactivate_player(HEAP);
    }

    //! Remove the extract buffer from the minima tree.
    void deactivate_extract_buffer()
    {
        m_head.deactivate_player(EB);
    }

    //! Remove an internal array from the minima tree.
    void deactivate_internal_array(unsigned index)
    {
        m_ia.deactivate_player(index);
        if (!m_ia.empty())
            m_head.notify_change(IA);
        else
            m_head.deactivate_player(IA);
    }

    //! Remove all insertion heaps from the minima tree.
    void clear_heaps()
    {
        m_heaps.clear();
        m_head.deactivate_player(HEAP);
    }

    //! Remove all internal arrays from the minima tree.
    void clear_internal_arrays()
    {
        m_ia.resize_and_clear(initial_ia_size);
        m_head.deactivate_player(IA);
    }

    void rebuild_internal_arrays()
    {
        if (!m_parent.m_internal_arrays.empty())
        {
            m_ia.resize_and_rebuild(m_parent.m_internal_arrays.size());
            m_head.notify_change(IA);
        }
        else
        {
            m_head.deactivate_player(IA);
        }
    }

    //! Return size of internal arrays minima tree
    size_t ia_slots() const
    {
        return m_ia.num_slots();
    }

    //! Returns a readable representation of the winner tree as string.
    std::string to_string() const
    {
        std::ostringstream ss;
        ss << "Head:" << std::endl << m_head.to_string() << std::endl;
        ss << "Heaps:" << std::endl << m_heaps.to_string() << std::endl;
        ss << "IA:" << std::endl << m_ia.to_string() << std::endl;
        return ss.str();
    }

    //! Prints statistical data.
    void print_stats() const
    {
        STXXL_MSG("Head winner tree stats:");
        m_head.print_stats();
        STXXL_MSG("Heaps winner tree stats:");
        m_heaps.print_stats();
        STXXL_MSG("IA winner tree stats:");
        m_ia.print_stats();
    }
};

} // namespace ppq_local

/*!
 * Parallelized External Memory Priority Queue.
 *
 * \tparam ValueType Type of the contained objects (POD with no references to
 * internal memory).
 *
 * \tparam CompareType The comparator type used to determine whether one
 * element is smaller than another element.
 *
 * \tparam DefaultMemSize Maximum memory consumption by the queue. Can be
 * overwritten by the constructor. Default = 1 GiB.
 *
 * \tparam MaxItems Maximum number of elements the queue contains at one
 * time. Default = 0 = unlimited. This is no hard limit and only used for
 * optimization. Can be overwritten by the constructor.
 *
 * \tparam BlockSize External block size. Default =
 * STXXL_DEFAULT_BLOCK_SIZE(ValueType).
 *
 * \tparam AllocStrategy Allocation strategy for the external memory. Default =
 * STXXL_DEFAULT_ALLOC_STRATEGY.
 */
template <
    class ValueType,
    class CompareType = std::less<ValueType>,
    class AllocStrategy = STXXL_DEFAULT_ALLOC_STRATEGY,
    uint64 BlockSize = STXXL_DEFAULT_BLOCK_SIZE(ValueType),
    uint64 DefaultMemSize = 1* 1024L* 1024L* 1024L,
    uint64 MaxItems = 0
    >
class parallel_priority_queue : private noncopyable
{
    //! \name Types
    //! \{

public:
    typedef ValueType value_type;
    typedef CompareType compare_type;
    typedef AllocStrategy alloc_strategy;
    static const uint64 block_size = BlockSize;
    typedef uint64 size_type;

    typedef typed_block<block_size, value_type> block_type;
    typedef std::vector<BID<block_size> > bid_vector;
    typedef bid_vector bids_container_type;
    typedef read_write_pool<block_type> pool_type;
    typedef ppq_local::internal_array<value_type> internal_array_type;
    typedef ppq_local::external_array<value_type, block_size, AllocStrategy> external_array_type;
    typedef typename external_array_type::writer_type external_array_writer_type;
    typedef typename std::vector<value_type>::iterator value_iterator;
    typedef typename internal_array_type::iterator iterator;
    typedef std::pair<iterator, iterator> iterator_pair_type;

    static const bool debug = false;

    //! currently global public tuning parameter:
    unsigned_type c_max_internal_level_size;

    //! currently global public tuning parameter:
    unsigned_type c_max_external_level_size;

protected:
    //! type of insertion heap itself
    typedef std::vector<value_type> heap_type;

    //! type of internal arrays vector
    typedef typename stxxl::swap_vector<internal_array_type> internal_arrays_type;
    //! type of external arrays vector
    typedef typename stxxl::swap_vector<external_array_type> external_arrays_type;
    //! type of minima tree combining the structures
    typedef ppq_local::minima_tree<
            parallel_priority_queue<value_type, compare_type, alloc_strategy,
                                    block_size, DefaultMemSize, MaxItems> > minima_type;
    //! allow minima tree access to internal data structures
    friend class ppq_local::minima_tree<
            parallel_priority_queue<value_type, compare_type, alloc_strategy,
                                    block_size, DefaultMemSize, MaxItems> >;

    //! Inverse comparison functor
    struct inv_compare_type
    {
        const compare_type& compare;

        inv_compare_type(const compare_type& c)
            : compare(c)
        { }

        bool operator () (const value_type& x, const value_type& y) const
        {
            return compare(y, x);
        }
    };

    //! <-Comparator for value_type
    compare_type m_compare;

    //! >-Comparator for value_type
    inv_compare_type m_inv_compare;

    //! Defines if statistics are gathered: dummy_custom_stats_counter or
    //! custom_stats_counter
    typedef dummy_custom_stats_counter<uint64> stats_counter;

    //! Defines if statistics are gathered: fake_timer or timer
    typedef fake_timer stats_timer;

    //! \}

    //! \name Compile-Time Parameters
    //! \{

    //! Merge sorted heaps when flushing into an internal array.
    //! Pro: Reduces the risk of a large winner tree
    //! Con: Flush insertion heaps becomes slower.
    static const bool c_merge_sorted_heaps = true;

    //! Default number of write buffer block for a new external array being
    //! filled.
    static const unsigned c_num_write_buffer_blocks = 14;

    //! Defines for how much external arrays memory should be reserved in the
    //! constructor.
    static const unsigned c_num_reserved_external_arrays = 10;

    //! Size of a single insertion heap in Byte, if not defined otherwise in
    //! the constructor. Default: 1 MiB
    static const size_type c_default_single_heap_ram = 1L * 1024L * 1024L;

    //! Default limit of the extract buffer ram consumption as share of total
    //! ram
    // C++11: constexpr static double c_default_extract_buffer_ram_part = 0.05;
    // C++98 does not allow static const double initialization here.
    // It's located in global scope instead.
    static const double c_default_extract_buffer_ram_part;

    /*!
     * Limit the size of the extract buffer to an absolute value.
     *
     * The actual size can be set using the extract_buffer_ram parameter of the
     * constructor. If this parameter is not set, the value is calculated by
     * (total_ram*c_default_extract_buffer_ram_part)
     *
     * If c_limit_extract_buffer==false, the memory consumption of the extract
     * buffer is only limited by the number of external and internal
     * arrays. This is considered in memory management using the
     * ram_per_external_array and ram_per_internal_array values. Attention:
     * Each internal array reserves space for the extract buffer in the size of
     * all heaps together.
     */
    static const bool c_limit_extract_buffer = true;

    //! For bulks of size up to c_single_insert_limit sequential single insert
    //! is faster than bulk_push.
    static const unsigned c_single_insert_limit = 100;

    //! \}

    //! \name Parameters and Sizes for Memory Allocation Policy

    //! Number of insertion heaps. Usually equal to the number of CPUs.
    const long m_num_insertion_heaps;

    //! Capacity of one inserion heap
    const unsigned m_insertion_heap_capacity;

    //! Return size of insertion heap reservation in bytes
    size_type insertion_heap_int_memory() const
    {
        return m_insertion_heap_capacity * sizeof(value_type);
    }

    //! Total amount of internal memory
    const size_type m_mem_total;

    //! Maximum size of extract buffer in number of elements
    //! Only relevant if c_limit_extract_buffer==true
    size_type m_extract_buffer_limit;

    //! Size of all insertion heaps together in bytes
    const size_type m_mem_for_heaps;

    //! Number of read/prefetch blocks per external array.
    const float m_num_read_blocks_per_ea;

    //! Total number of read/prefetch buffer blocks
    unsigned_type m_num_read_blocks;
    //! number of currently hinted prefetch blocks
    unsigned_type m_num_hinted_blocks;
    //! number of currently loaded blocks
    unsigned_type m_num_used_read_blocks;

    //! Free memory in bytes
    size_type m_mem_left;

    //! \}

    //! Flag if inside a bulk_push sequence.
    bool m_in_bulk_push;

    //! If the bulk currently being inserted is very large, this boolean is set
    //! and bulk_push just accumulate the elements for eventual sorting.
    bool m_is_very_large_bulk;

    //! First index in m_external_arrays that was not re-hinted during a
    //! bulk_push sequence.
    unsigned_type m_bulk_first_delayed_external_array;

    //! Index of the currently smallest element in the extract buffer
    size_type m_extract_buffer_index;

    //! \name Number of elements currently in the data structures
    //! \{

    //! Number of elements int the insertion heaps
    size_type m_heaps_size;

    //! Number of elements in the extract buffer
    size_type m_extract_buffer_size;

    //! Number of elements in the internal arrays
    size_type m_internal_size;

    //! Number of elements in the external arrays
    size_type m_external_size;

    //! \}

    //! \name Data Holding Structures
    //! \{

    //! A struct containing the local insertion heap and other information
    //! _local_ to a processor.
    struct ProcessorData
    {
        //! The heaps where new elements are usually inserted into
        heap_type insertion_heap;

        //! The number of items inserted into the insheap during bulk parallel
        //! access.
        size_type heap_add_size;
    };

    typedef std::vector<ProcessorData*> proc_vector_type;

    //! Array of processor local data structures, including the insertion heaps.
    proc_vector_type m_proc;

    //! Prefetch and write buffer pool for external arrays (has to be in front
    //! of m_external_arrays)
    pool_type m_pool;

    //! The extract buffer where external (and internal) arrays are merged into
    //! for extracting
    std::vector<value_type> m_extract_buffer;

    //! The sorted arrays in internal memory
    internal_arrays_type m_internal_arrays;

    //! The sorted arrays in external memory
    external_arrays_type m_external_arrays;

    //! The aggregated pushes. They cannot be extracted yet.
    std::vector<value_type> m_aggregated_pushes;

    //! The maximum number of internal array levels.
    static const unsigned_type c_max_internal_levels = 8;

    //! The number of internal arrays on each level, we use plain array.
    unsigned_type m_internal_levels[c_max_internal_levels];

    //! The maximum number of external array levels.
    static const unsigned_type c_max_external_levels = 8;

    //! The number of external arrays on each level, we use plain array.
    unsigned_type m_external_levels[c_max_external_levels];

    //! The winner tree containing the smallest values of all sources
    //! where the globally smallest element could come from.
    minima_type m_minima;

    //! Compares the largest accessible value of two external arrays.
    struct external_min_comparator {
        const external_arrays_type& m_eas;
        const inv_compare_type& m_compare;

        external_min_comparator(const external_arrays_type& eas,
                                const inv_compare_type& compare)
            : m_eas(eas), m_compare(compare) { }

        bool operator () (const size_t& a, const size_t& b) const
        {
            return m_compare(m_eas[a].get_next_block_min(),
                             m_eas[b].get_next_block_min());
        }
    } m_external_min_comparator;

    //! Tracks the largest accessible values of the external arrays if there
    //! is unaccessible data in EM. The winning array is the first one that
    //! needs to fetch further data from EM. Used in calculate_merge_sequences.
    winner_tree<external_min_comparator> m_external_min_tree;

    //! Compares the largest value of the block hinted the latest of two
    //! external arrays.
    struct hint_comparator {
        const external_arrays_type& m_eas;
        const inv_compare_type& m_compare;

        hint_comparator(const external_arrays_type& eas,
                        const inv_compare_type& compare)
            : m_eas(eas), m_compare(compare) { }

        bool operator () (const size_t& a, const size_t& b) const
        {
            return m_compare(m_eas[a].get_next_hintable_min(),
                             m_eas[b].get_next_hintable_min());
        }
    } m_hint_comparator;

    //! Tracks the largest values of the block hinted the latest of the
    //! external arrays if there is unaccessible data in EM. The winning
    //! array is the first one that needs to fetch further data from EM.
    //! Used for prefetch hints.
    winner_tree<hint_comparator> m_hint_tree;

    //! Random number generator for randomly selecting a heap in sequential
    //! push()
    random_number32_r m_rng;

    //! \}

    /*
     * Helper function to remove empty internal/external arrays.
     */

    //! Unary operator which returns true if the external array has run empty.
    struct empty_external_array_eraser {
        bool operator () (external_array_type& a) const
        { return a.empty(); }
    };

    //! Unary operator which returns true if the internal array has run empty.
    struct empty_internal_array_eraser {
        bool operator () (internal_array_type& a) const
        { return a.empty(); }
    };

    //! Clean up empty internal arrays, free their memory and capacity
    void cleanup_internal_arrays()
    {
        typename internal_arrays_type::iterator swap_end =
            stxxl::swap_remove_if(m_internal_arrays.begin(),
                                  m_internal_arrays.end(),
                                  empty_internal_array_eraser());

        for (typename internal_arrays_type::iterator ia = swap_end;
             ia != m_internal_arrays.end(); ++ia)
        {
            m_mem_left += ia->int_memory();
            --m_internal_levels[ia->level()];
        }

        if (swap_end != m_internal_arrays.end())
            STXXL_DEBUG0("cleanup_internal_arrays" <<
                         " cleaned=" << m_internal_arrays.end() - swap_end);

        m_internal_arrays.erase(swap_end, m_internal_arrays.end());
        m_minima.rebuild_internal_arrays();
    }

    //! Clean up empty external arrays, free their memory and capacity
    void cleanup_external_arrays()
    {
        typedef typename external_arrays_type::iterator ea_iterator;
        empty_external_array_eraser pred;

        // The following is a modified implementation of swap_remove_if().
        // Updates m_external_min_tree accordingly.

        ea_iterator first = m_external_arrays.begin();
        ea_iterator last = m_external_arrays.end();
        ea_iterator swap_end = first;
        size_t size = m_external_arrays.end() - m_external_arrays.begin();
        size_t first_removed = size;
        while (first != last)
        {
            if (!pred(*first))
            {
                using std::swap;
                swap(*first, *swap_end);
                ++swap_end;
            }
            else if (first_removed >= size)
            {
                first_removed = first - m_external_arrays.begin();
            }
            ++first;
        }

        // subtract memory of EAs, which will be freed
        for (ea_iterator ea = swap_end; ea != last; ++ea) {
            m_mem_left += ea->int_memory();
            --m_external_levels[ea->level()];
        }

        size_t swap_end_index = swap_end - m_external_arrays.begin();

        // Deactivating all affected players first.
        // Otherwise there might be outdated comparisons.
        for (size_t i = size; i != first_removed; ) {
            --i;
            m_external_min_tree.deactivate_player_step(i);
            // TODO delay if (m_in_bulk_push)?
            m_hint_tree.deactivate_player_step(i);
        }

        // Replay moved arrays.
        for (size_t i = first_removed; i < swap_end_index; ++i) {
            update_external_min_tree(i);
            // TODO delay if (m_in_bulk_push)?
            update_hint_tree(i);
        }

        STXXL_DEBUG("Removed " << m_external_arrays.end() - swap_end <<
                    " empty external arrays.");

        m_external_arrays.erase(swap_end, m_external_arrays.end());

        resize_read_pool(); // shrinks read/prefetch pool
    }

    /*!
     * SiftUp a new element from the last position in the heap, reestablishing
     * the heap invariant. This is identical to std::push_heap, except that it
     * returns the last element modified by siftUp. Thus we can identify if the
     * minimum may have changed.
     */
    template <typename RandomAccessIterator, typename HeapCompareType>
    static inline unsigned_type
    push_heap(RandomAccessIterator first, RandomAccessIterator last,
              HeapCompareType comp)
    {
        typedef typename std::iterator_traits<RandomAccessIterator>::value_type
            value_type;

        value_type value = STXXL_MOVE(*(last - 1));

        unsigned_type index = (last - first) - 1;
        unsigned_type parent = (index - 1) / 2;

        while (index > 0 && comp(*(first + parent), value))
        {
            *(first + index) = STXXL_MOVE(*(first + parent));
            index = parent;
            parent = (index - 1) / 2;
        }
        *(first + index) = STXXL_MOVE(value);

        return index;
    }

public:
    //! \name Initialization
    //! \{

    /*!
     * Constructor.
     *
     * \param compare Comparator for priority queue, which is a Max-PQ.
     *
     * \param total_ram Maximum RAM usage. 0 = Default = Use the template
     * value DefaultMemSize.
     *
     * \param num_read_blocks_per_ea Number of read blocks per external
     * array. Default = 1.5f
     *
     * \param num_write_buffer_blocks Number of write buffer blocks for a new
     * external array being filled. 0 = Default = c_num_write_buffer_blocks
     *
     * \param num_insertion_heaps Number of insertion heaps. 0 = Default =
     * Determine by omp_get_max_threads().
     *
     * \param single_heap_ram Memory usage for a single insertion heap.
     * Default = c_single_heap_ram.
     *
     * \param extract_buffer_ram Memory usage for the extract buffer. Only
     * relevant if c_limit_extract_buffer==true. 0 = Default = total_ram *
     * c_default_extract_buffer_ram_part.
     */
    parallel_priority_queue(
        const compare_type& compare = compare_type(),
        size_type total_ram = DefaultMemSize,
        float num_read_blocks_per_ea = 1.5f,
        unsigned_type num_write_buffer_blocks = c_num_write_buffer_blocks,
        unsigned_type num_insertion_heaps = 0,
        size_type single_heap_ram = c_default_single_heap_ram,
        size_type extract_buffer_ram = 0)
        : c_max_internal_level_size(64),
          c_max_external_level_size(64),
          m_compare(compare),
          m_inv_compare(m_compare),
          // Parameters and Sizes for Memory Allocation Policy
#if STXXL_PARALLEL
          m_num_insertion_heaps(num_insertion_heaps > 0 ? num_insertion_heaps : omp_get_max_threads()),
#else
          m_num_insertion_heaps(num_insertion_heaps > 0 ? num_insertion_heaps : 1),
#endif
          m_insertion_heap_capacity(single_heap_ram / sizeof(value_type)),
          m_mem_total(total_ram),
          m_mem_for_heaps(m_num_insertion_heaps * single_heap_ram),
          m_num_read_blocks_per_ea(num_read_blocks_per_ea),
          m_num_read_blocks(0),
          m_num_hinted_blocks(0),
          m_num_used_read_blocks(0),
          // (unnamed)
          m_in_bulk_push(false),
          m_is_very_large_bulk(false),
          m_extract_buffer_index(0),
          // Number of elements currently in the data structures
          m_heaps_size(0),
          m_extract_buffer_size(0),
          m_internal_size(0),
          m_external_size(0),
          // Data Holding Structures
          m_proc(m_num_insertion_heaps),
          m_pool(0, num_write_buffer_blocks),
          m_external_arrays(),
          m_minima(*this),
          m_external_min_comparator(m_external_arrays, m_inv_compare),
          m_external_min_tree(4, m_external_min_comparator),
          m_hint_comparator(m_external_arrays, m_inv_compare),
          m_hint_tree(4, m_hint_comparator),
          // flags
          m_limit_extract(false)
    {
#if STXXL_PARALLEL
        if (!omp_get_nested()) {
            omp_set_nested(1);
            if (!omp_get_nested()) {
                STXXL_ERRMSG("Could not enable OpenMP's nested parallelism, "
                             "however, the PPQ requires this OpenMP feature.");
                abort();
            }
        }
#else
        STXXL_ERRMSG("You are using stxxl::parallel_priority_queue without "
                     "support for OpenMP parallelism.");
        STXXL_ERRMSG("This is probably not what you want, so check the "
                     "compilation settings.");
#endif

        if (c_limit_extract_buffer) {
            m_extract_buffer_limit = (extract_buffer_ram > 0)
                                     ? extract_buffer_ram / sizeof(value_type)
                                     : static_cast<size_type>(((double)(m_mem_total) * c_default_extract_buffer_ram_part / sizeof(value_type)));
        }

        for (unsigned_type i = 0; i < c_max_internal_levels; ++i)
            m_internal_levels[i] = 0;

        for (unsigned_type i = 0; i < c_max_external_levels; ++i)
            m_external_levels[i] = 0;

        // TODO: Do we still need this line? Insertion heap memory is
        // registered below. And merge buffer is equal to the new IA...
        // total_ram - ram for the heaps - ram for the heap merger
        m_mem_left = m_mem_total - 2 * m_mem_for_heaps;

        // reverse insertion heap memory on processor-local memory
#if STXXL_PARALLEL
#pragma omp parallel for
#endif
        for (long p = 0; p < m_num_insertion_heaps; ++p)
        {
            m_proc[p] = new ProcessorData;
            m_proc[p]->insertion_heap.reserve(m_insertion_heap_capacity);
            assert(m_proc[p]->insertion_heap.capacity() * sizeof(value_type)
                   == insertion_heap_int_memory());
        }

        m_mem_left -= m_num_insertion_heaps * insertion_heap_int_memory();

        // prepare prefetch buffer pool (already done in initializer),
        // initially zero.

        // prepare write buffer pool: calculate size and subtract from mem_left
        external_array_type::prepare_write_pool(m_pool, m_num_insertion_heaps);
        m_mem_left -= m_pool.size_write() * block_size;

        // prepare internal arrays
        if (c_merge_sorted_heaps) {
            m_internal_arrays.reserve(m_mem_total / m_mem_for_heaps);
        }
        else {
            m_internal_arrays.reserve(m_mem_total * m_num_insertion_heaps / m_mem_for_heaps);
        }

        // prepare external arrays
        m_external_arrays.reserve(c_num_reserved_external_arrays);

        if (m_mem_total < m_mem_left) // checks if unsigned type wrapped.
        {
            STXXL_ERRMSG("Minimum memory requirement insufficient, "
                         "increase PPQ's memory limit or decrease buffers.");
            abort();
        }

        check_invariants();
    }

    //! Destructor.
    ~parallel_priority_queue()
    {
        // clean up data structures

        for (size_t p = 0; p < m_num_insertion_heaps; ++p)
        {
            delete m_proc[p];
        }
    }

protected:
    //! Assert many invariants of the data structures.
    void check_invariants() const
    {
#ifdef NDEBUG
        // disable in Release builds
        return;
#endif

        size_type mem_used = 0;

        mem_used += 2 * m_mem_for_heaps
                    + m_pool.size_write() * block_size
                    + m_pool.free_size_prefetch() * block_size
                    + m_num_hinted_blocks * block_size
                    + m_num_used_read_blocks * block_size;

        // count number of blocks hinted in prefetcher

        size_t num_hinted = 0, num_used_read = 0;
        for (size_t i = 0; i < m_external_arrays.size(); ++i) {
            num_hinted += m_external_arrays[i].num_hinted_blocks();
            num_used_read += m_external_arrays[i].num_used_blocks();
        }

        STXXL_CHECK(num_hinted == m_num_hinted_blocks);
        STXXL_CHECK(num_used_read == m_num_used_read_blocks);

        STXXL_CHECK_EQUAL(m_num_used_read_blocks,
                          m_num_read_blocks
                          - m_pool.free_size_prefetch()
                          - m_num_hinted_blocks);

        // test the processor local data structures

        size_type heaps_size = 0;

        for (int_type p = 0; p < m_num_insertion_heaps; ++p)
        {
            // check that each insertion heap is a heap

            // TODO: remove soon, because this is very expensive
            STXXL_CHECK(1 || stxxl::is_heap(m_proc[p]->insertion_heap.begin(),
                                            m_proc[p]->insertion_heap.end(),
                                            m_compare));

            STXXL_CHECK(m_proc[p]->insertion_heap.capacity() <= m_insertion_heap_capacity);

            heaps_size += m_proc[p]->insertion_heap.size();
            mem_used += m_proc[p]->insertion_heap.capacity() * sizeof(value_type);
        }

        if (!m_in_bulk_push)
            STXXL_CHECK_EQUAL(m_heaps_size, heaps_size);

        // count number of items and memory size of internal arrays

        size_type ia_size = 0;
        size_type ia_memory = 0;
        std::vector<unsigned_type> ia_levels(c_max_internal_levels, 0);

        for (typename internal_arrays_type::const_iterator ia =
                 m_internal_arrays.begin(); ia != m_internal_arrays.end(); ++ia)
        {
            ia_size += ia->size();
            ia_memory += ia->int_memory();
            ++ia_levels[ia->level()];
        }

        STXXL_CHECK_EQUAL(m_internal_size, ia_size);
        mem_used += ia_memory;

        for (unsigned_type i = 0; i < c_max_internal_levels; ++i)
            STXXL_CHECK_EQUAL(m_internal_levels[i], ia_levels[i]);

        // count number of items in external arrays

        size_type ea_size = 0;
        size_type ea_memory = 0;
        std::vector<unsigned_type> ea_levels(c_max_external_levels, 0);

        for (typename external_arrays_type::const_iterator ea =
                 m_external_arrays.begin(); ea != m_external_arrays.end(); ++ea)
        {
            ea_size += ea->size();
            ea_memory += ea->int_memory();
            ++ea_levels[ea->level()];
        }

        STXXL_CHECK_EQUAL(m_external_size, ea_size);
        mem_used += ea_memory;

        for (unsigned_type i = 0; i < c_max_external_levels; ++i)
            STXXL_CHECK_EQUAL(m_external_levels[i], ea_levels[i]);

        // calculate mem_used so that == mem_total - mem_left

        STXXL_CHECK_EQUAL(memory_consumption(), mem_used);
    }

    //! \}

    //! \name Properties
    //! \{

public:
    //! The number of elements in the queue.
    inline size_type size() const
    {
        return m_heaps_size + m_internal_size + m_external_size + m_extract_buffer_size;
    }

    //! Returns if the queue is empty.
    inline bool empty() const
    {
        return (size() == 0);
    }

    //! The memory consumption in Bytes.
    inline size_type memory_consumption() const
    {
        assert(m_mem_total >= m_mem_left);
        return (m_mem_total - m_mem_left);
    }

protected:
    //! Returns if the extract buffer is empty.
    inline bool extract_buffer_empty() const
    {
        return (m_extract_buffer_size == 0);
    }

    //! \}

public:
    //! \name Bulk Operations
    //! \{

    /*!
     * Start a sequence of push operations.
     * \param bulk_size Exact number of elements to push before the next pop.
     */
    void bulk_push_begin(size_type bulk_size)
    {
        assert(!m_in_bulk_push);
        m_in_bulk_push = true;
        m_bulk_first_delayed_external_array = m_external_arrays.size();

        size_type heap_capacity = m_num_insertion_heaps * m_insertion_heap_capacity;

        // if bulk_size is large: use simple aggregation instead of keeping the
        // heap property and sort everything afterwards.
        if (bulk_size > heap_capacity && 0) {
            m_is_very_large_bulk = true;
        }
        else {
            m_is_very_large_bulk = false;

            if (bulk_size + m_heaps_size > heap_capacity) {
                if (m_heaps_size > 0) {
                    //flush_insertion_heaps();
                }
            }
        }

        // zero bulk insertion counters
        for (int_type p = 0; p < m_num_insertion_heaps; ++p)
            m_proc[p]->heap_add_size = 0;
    }

    /*!
     * Push an element inside a sequence of pushes.
     * Run bulk_push_begin() before using this method.
     *
     * \param element The element to push.
     * \param p The id of the insertion heap to use (usually the thread id).
     */
    void bulk_push(const value_type& element, const unsigned_type p)
    {
        assert(m_in_bulk_push);

        heap_type& insheap = m_proc[p]->insertion_heap;

        if (!m_is_very_large_bulk && 0)
        {
            // if small bulk: if heap is full -> sort locally and put into
            // internal array list. insert items and keep heap invariant.
            if (UNLIKELY(insheap.size() >= m_insertion_heap_capacity)) {
#if STXXL_PARALLEL
#pragma omp atomic
#endif
                m_heaps_size += m_proc[p]->heap_add_size;

                m_proc[p]->heap_add_size = 0;
                flush_insertion_heap(p);
            }

            assert(insheap.size() < insheap.capacity());

            // put item onto heap and siftUp
            insheap.push_back(element);
            std::push_heap(insheap.begin(), insheap.end(), m_compare);
        }
        else if (!m_is_very_large_bulk && 1)
        {
            // if small bulk: if heap is full -> sort locally and put into
            // internal array list. insert items but DO NOT keep heap
            // invariant.
            if (UNLIKELY(insheap.size() >= m_insertion_heap_capacity)) {
#if STXXL_PARALLEL
#pragma omp atomic
#endif
                m_heaps_size += m_proc[p]->heap_add_size;

                m_proc[p]->heap_add_size = 0;
                flush_insertion_heap(p);
            }

            assert(insheap.size() < insheap.capacity());

            // put item onto heap and DO NOT siftUp
            insheap.push_back(element);
        }
        else // m_is_very_large_bulk
        {
            if (UNLIKELY(insheap.size() >= 2 * 1024 * 1024)) {
#if STXXL_PARALLEL
#pragma omp atomic
#endif
                m_heaps_size += m_proc[p]->heap_add_size;

                m_proc[p]->heap_add_size = 0;
                flush_insertion_heap(p);
            }

            assert(insheap.size() < insheap.capacity());

            // put onto insertion heap but do not keep heap property
            insheap.push_back(element);
        }

        m_proc[p]->heap_add_size++;
    }

    /*!
     * Push an element inside a bulk sequence of pushes. Run bulk_push_begin()
     * before using this method. This function uses the insertion heap id =
     * omp_get_thread_num().
     *
     * \param element The element to push.
     */
    void bulk_push(const value_type& element)
    {
#if STXXL_PARALLEL
        return bulk_push(element, (unsigned_type)omp_get_thread_num());
#else
        unsigned_type id = m_rng() % m_num_insertion_heaps;
        return bulk_push(element, id);
#endif
    }

    /*!
     * Ends a sequence of push operations. Run bulk_push_begin() and some
     * bulk_push() before this.
     */
    void bulk_push_end()
    {
        assert(m_in_bulk_push);
        m_in_bulk_push = false;

        if (!m_is_very_large_bulk && 0)
        {
            for (int_type p = 0; p < m_num_insertion_heaps; ++p)
            {
                m_heaps_size += m_proc[p]->heap_add_size;

                if (!m_proc[p]->insertion_heap.empty())
                    m_minima.update_heap(p);
            }
        }
        else if (!m_is_very_large_bulk && 1)
        {
#if STXXL_PARALLEL
#pragma omp parallel for
#endif
            for (int_type p = 0; p < m_num_insertion_heaps; ++p)
            {
                // reestablish heap property: siftUp only those items pushed
                for (unsigned_type index = m_proc[p]->heap_add_size; index != 0; ) {
                    std::push_heap(m_proc[p]->insertion_heap.begin(),
                                   m_proc[p]->insertion_heap.end() - (--index),
                                   m_compare);
                }

#if STXXL_PARALLEL
#pragma omp atomic
#endif
                m_heaps_size += m_proc[p]->heap_add_size;
            }

            for (int_type p = 0; p < m_num_insertion_heaps; ++p)
            {
                if (!m_proc[p]->insertion_heap.empty())
                    m_minima.update_heap(p);
            }
        }
        else // m_is_very_large_bulk
        {
#if STXXL_PARALLEL
#pragma omp parallel for
#endif
            for (int_type p = 0; p < m_num_insertion_heaps; ++p)
            {
                if (m_proc[p]->insertion_heap.size() >= m_insertion_heap_capacity) {
                    // flush out overfull insertion heap arrays
#if STXXL_PARALLEL
#pragma omp atomic
#endif
                    m_heaps_size += m_proc[p]->heap_add_size;

                    m_proc[p]->heap_add_size = 0;
                    flush_insertion_heap(p);
                }
                else {
                    // reestablish heap property: siftUp only those items pushed
                    for (unsigned_type index = m_proc[p]->heap_add_size; index != 0; ) {
                        std::push_heap(m_proc[p]->insertion_heap.begin(),
                                       m_proc[p]->insertion_heap.end() - (--index),
                                       m_compare);
                    }

#if STXXL_PARALLEL
#pragma omp atomic
#endif
                    m_heaps_size += m_proc[p]->heap_add_size;
                    m_proc[p]->heap_add_size = 0;
                }
            }

            for (int_type p = 0; p < m_num_insertion_heaps; ++p)
            {
                if (!m_proc[p]->insertion_heap.empty())
                    m_minima.update_heap(p);
            }
        }

        if (m_bulk_first_delayed_external_array != m_external_arrays.size()) {
            STXXL_DEBUG("bulk_push_end: run delayed re-hinting of EAs");
            rebuild_hint_tree();
        }

        check_invariants();
    }

    //! Extract up to max_size values at once.
    void bulk_pop(std::vector<value_type>& out, size_t max_size)
    {
        STXXL_DEBUG("bulk_pop_size with max_size=" << max_size);

        const size_t n_elements = std::min<size_t>(max_size, size());
        assert(n_elements < m_extract_buffer_limit);

        if (m_heaps_size > 0)
            flush_insertion_heaps();

        convert_eb_into_ia();

        refill_extract_buffer(n_elements, n_elements);

        out.resize(0);
        using std::swap;
        swap(m_extract_buffer, out);
        m_extract_buffer_index = 0;
        m_extract_buffer_size = 0;
        m_minima.deactivate_extract_buffer();

        check_invariants();
    }

    //! Extracts all elements which are greater or equal to a given limit.
    //! \param out result vector
    //! \param limit limit value
    //! \param max_size maximum number of items to extract
    //! \return true if the buffer contains all items < limit, false it was too
    //! small.
    bool bulk_pop_limit(std::vector<value_type>& out, const value_type& limit,
                        size_t max_size = std::numeric_limits<size_t>::max())
    {
        STXXL_DEBUG("bulk_pop_limit with limit=" << limit);

        convert_eb_into_ia();

        if (m_heaps_size > 0) {
            if (0)
                flush_insertion_heaps();
            else if (1)
                flush_insertion_heaps_with_limit(limit);
        }

        size_type ias = m_internal_arrays.size();
        size_type eas = m_external_arrays.size();
        std::vector<size_type> sizes(eas + ias);
        std::vector<iterator_pair_type> sequences(eas + ias);
        size_type output_size = 0;

        int limiting_ea_index = m_external_min_tree.top();

        // pop limit may have to change due to memory limit
        value_type this_limit = limit;
        bool has_full_range = true;

        // get all relevant blocks
        while (limiting_ea_index > -1)
        {
            const value_type& ea_limit =
                m_external_arrays[limiting_ea_index].get_next_block_min();

            if (m_compare(ea_limit, this_limit)) {
                // No more EM data smaller or equal to limit
                break;
            }

            if (m_external_arrays[limiting_ea_index].num_hinted_blocks() == 0) {
                // No more read/prefetch blocks available for EA
                this_limit = ea_limit;
                has_full_range = false;
                break;
            }

            wait_next_ea_blocks(limiting_ea_index);
            // consider next limiting EA
            limiting_ea_index = m_external_min_tree.top();
            STXXL_ASSERT(limiting_ea_index < (int)eas);
        }

        // build sequences
        for (size_type i = 0; i < eas + ias; ++i) {
            iterator begin, end;

            if (i < eas) {
                assert(!m_external_arrays[i].empty());
                assert(m_external_arrays[i].valid());
                begin = m_external_arrays[i].begin();
                end = m_external_arrays[i].end();
            }
            else {
                size_type j = i - eas;
                assert(!(m_internal_arrays[j].empty()));
                begin = m_internal_arrays[j].begin();
                end = m_internal_arrays[j].end();
            }

            end = std::lower_bound(begin, end, this_limit, m_inv_compare);

            sizes[i] = std::distance(begin, end);
            sequences[i] = std::make_pair(begin, end);
        }

        output_size = std::accumulate(sizes.begin(), sizes.end(), 0);
        if (output_size > max_size) {
            output_size = max_size;
            has_full_range = false;
        }
        out.resize(output_size);

        STXXL_DEBUG("bulk_pop_limit with" <<
                    " sequences=" << sequences.size() <<
                    " output_size=" << output_size <<
                    " has_full_range=" << has_full_range);

        potentially_parallel::multiway_merge(
            sequences.begin(), sequences.end(),
            out.begin(), output_size, m_inv_compare);

        advance_arrays(sequences, sizes, eas, ias);

        check_invariants();

        return has_full_range;
    }

#if TODO_MAYBE_FIXUP_LATER
    /*!
     * Insert a vector of elements at one time.
     * \param elements Vector containing the elements to push.
     * Attention: elements vector may be owned by the PQ afterwards.
     */
    void bulk_push_vector(std::vector<value_type>& elements)
    {
        size_type heap_capacity = m_num_insertion_heaps * m_insertion_heap_capacity;
        if (elements.size() > heap_capacity / 2) {
            flush_array(elements);
            return;
        }

        bulk_push_begin(elements.size());
#if STXXL_PARALLEL
        #pragma omp parallel
        {
            const unsigned thread_num = omp_get_thread_num();
            #pragma omp parallel for
            for (size_type i = 0; i < elements.size(); ++i) {
                bulk_push(elements[i], thread_num);
            }
        }
#else
        const unsigned thread_num = m_rng() % m_num_insertion_heaps;
        for (size_type i = 0; i < elements.size(); ++i) {
            bulk_push(elements[i], thread_num);
        }
#endif
        bulk_push_end();
    }
#endif

    //! \}

    //! \name Aggregation Operations
    //! \{

    /*!
     * Aggregate pushes. Use flush_aggregated_pushes() to finally push
     * them. extract_min is allowed is allowed in between the aggregation of
     * pushes if you can assure, that the extracted value is smaller than all
     * of the aggregated values.
     * \param element The element to push.
     */
    void aggregate_push(const value_type& element)
    {
        m_aggregated_pushes.push_back(element);
    }

#if TODO_MAYBE_FIXUP_LATER
    /*!
     * Insert the aggregated values into the queue using push(), bulk insert,
     * or sorting, depending on the number of aggregated values.
     */
    void flush_aggregated_pushes()
    {
        size_type size = m_aggregated_pushes.size();
        size_type ram_internal = 2 * size * sizeof(value_type); // ram for the sorted array + part of the ram for the merge buffer
        size_type heap_capacity = m_num_insertion_heaps * m_insertion_heap_capacity;

        if (ram_internal > m_mem_for_heaps / 2) {
            flush_array(m_aggregated_pushes);
        }
        else if ((m_aggregated_pushes.size() > c_single_insert_limit) && (m_aggregated_pushes.size() < heap_capacity)) {
            bulk_push_vector(m_aggregated_pushes);
        }
        else {
            for (value_iterator i = m_aggregated_pushes.begin(); i != m_aggregated_pushes.end(); ++i) {
                push(*i);
            }
        }

        m_aggregated_pushes.clear();
    }
#endif
    //! \}

    //! \name std::priority_queue compliant operations
    //! \{

    /*!
     * Insert new element
     * \param element the element to insert.
     * \param p number of insertion heap to insert item into
     */
    void push(const value_type& element, unsigned_type p = 0)
    {
        assert(!m_in_bulk_push && !m_limit_extract);

        heap_type& insheap = m_proc[p]->insertion_heap;

        if (insheap.size() >= m_insertion_heap_capacity) {
            flush_insertion_heap(p);
        }

        // push item to end of heap and siftUp
        insheap.push_back(element);
        unsigned_type index = push_heap(insheap.begin(), insheap.end(),
                                        m_compare);
        ++m_heaps_size;

        if (insheap.size() == 1 || index == 0)
            m_minima.update_heap(p);
    }

    //! Access the minimum element.
    const value_type & top()
    {
        assert(!m_in_bulk_push && !m_limit_extract);
        assert(!empty());

        if (extract_buffer_empty()) {
            refill_extract_buffer(std::min(m_extract_buffer_limit,
                                           m_internal_size + m_external_size));
        }

        static const bool debug = false;

        std::pair<unsigned, unsigned> type_and_index = m_minima.top();
        const unsigned& type = type_and_index.first;
        const unsigned& index = type_and_index.second;

        assert(type < 4);

        switch (type) {
        case minima_type::HEAP:
            STXXL_DEBUG("heap " << index <<
                        ": " << m_proc[index]->insertion_heap[0]);
            return m_proc[index]->insertion_heap[0];
        case minima_type::IA:
            STXXL_DEBUG("ia " << index <<
                        ": " << m_internal_arrays[index].get_min());
            return m_internal_arrays[index].get_min();
        case minima_type::EB:
            STXXL_DEBUG("eb " << m_extract_buffer_index <<
                        ": " << m_extract_buffer[m_extract_buffer_index]);
            return m_extract_buffer[m_extract_buffer_index];
        default:
            STXXL_ERRMSG("Unknown extract type: " << type);
            abort();
        }
    }

    //! Remove the minimum element.
    void pop()
    {
        assert(!m_in_bulk_push && !m_limit_extract);

        m_stats.num_extracts++;

        if (extract_buffer_empty()) {
            refill_extract_buffer(std::min(m_extract_buffer_limit,
                                           m_internal_size + m_external_size));
        }

        m_stats.extract_min_time.start();

        std::pair<unsigned, unsigned> type_and_index = m_minima.top();
        unsigned type = type_and_index.first;
        unsigned index = type_and_index.second;

        assert(type < 4);

        switch (type) {
        case minima_type::HEAP:
        {
            heap_type& insheap = m_proc[index]->insertion_heap;

            m_stats.pop_heap_time.start();
            std::pop_heap(insheap.begin(), insheap.end(), m_compare);
            insheap.pop_back();
            m_stats.pop_heap_time.stop();

            m_heaps_size--;

            if (!insheap.empty())
                m_minima.update_heap(index);
            else
                m_minima.deactivate_heap(index);

            break;
        }
        case minima_type::IA:
        {
            m_internal_arrays[index].inc_min();
            m_internal_size--;

            if (!(m_internal_arrays[index].empty()))
                m_minima.update_internal_array(index);
            else
                // internal array has run empty
                m_minima.deactivate_internal_array(index);

            break;
        }
        case minima_type::EB:
        {
            ++m_extract_buffer_index;
            assert(m_extract_buffer_size > 0);
            --m_extract_buffer_size;

            if (!extract_buffer_empty())
                m_minima.update_extract_buffer();
            else
                m_minima.deactivate_extract_buffer();

            break;
        }
        default:
            STXXL_ERRMSG("Unknown extract type: " << type);
            abort();
        }

        m_stats.extract_min_time.stop();

        check_invariants();
    }

    //! \}

    //! \name Bulk-Limit Operations
    //! \{

protected:
    //! current limit element
    value_type m_limit_element;

    //! flag if inside a bulk limit extract session
    bool m_limit_extract;

    //! flag if the extract buffer contains the full limit range
    bool m_limit_has_full_range;

public:
    //! Begin bulk-limit extraction session with limit element.
    void limit_begin(const value_type& limit, size_type bulk_size)
    {
        m_limit_extract = true;
        m_limit_element = limit;

        std::vector<value_type> new_extract_buffer;
        m_limit_has_full_range =
            bulk_pop_limit(new_extract_buffer, limit, m_extract_buffer_limit);
        std::swap(new_extract_buffer, m_extract_buffer);

        m_extract_buffer_index = 0;
        m_extract_buffer_size = m_extract_buffer.size();
        if (m_extract_buffer_size)
            m_minima.update_extract_buffer();
        else
            m_minima.deactivate_extract_buffer();

        bulk_push_begin(bulk_size);
    }

    //! Push new item >= bulk-limit element into insertion heap p.
    void limit_push(const value_type& element, const unsigned_type p = 0)
    {
        assert(m_limit_extract);
        assert(!m_compare(m_limit_element, element));

        return bulk_push(element, p);
    }

    //! Access the minimum element, which can only be in the extract buffer.
    const value_type & limit_top()
    {
        assert(m_limit_extract);

        // if buffer is empty and we extracted the full range last time, return
        // limit items as sentinel.
        if (m_extract_buffer_size == 0 && m_limit_has_full_range)
            return m_limit_element;

        if (extract_buffer_empty())
        {
            // extract more items
            std::vector<value_type> new_extract_buffer;
            m_limit_has_full_range =
                bulk_pop_limit(new_extract_buffer, m_limit_element,
                               m_extract_buffer_limit);

            std::swap(new_extract_buffer, m_extract_buffer);

            m_extract_buffer_index = 0;
            m_extract_buffer_size = m_extract_buffer.size();
            if (m_extract_buffer_size)
                m_minima.update_extract_buffer();
            else
                m_minima.deactivate_extract_buffer();
        }

        return m_extract_buffer[m_extract_buffer_index];
    }

    //! Remove the minimum element, only works correctly while elements < L.
    void limit_pop()
    {
        assert(m_limit_extract);

        ++m_extract_buffer_index;
        assert(m_extract_buffer_size > 0);
        --m_extract_buffer_size;

        if (extract_buffer_empty() && !m_limit_has_full_range)
        {
            // extract more items
            std::vector<value_type> new_extract_buffer;
            m_limit_has_full_range =
                bulk_pop_limit(new_extract_buffer, m_limit_element,
                               m_extract_buffer_limit);

            std::swap(new_extract_buffer, m_extract_buffer);

            m_extract_buffer_index = 0;
            m_extract_buffer_size = m_extract_buffer.size();
            if (m_extract_buffer_size)
                m_minima.update_extract_buffer();
            else
                m_minima.deactivate_extract_buffer();
        }
    }

    //! Finish bulk-limit extraction session.
    void limit_end()
    {
        assert(m_limit_extract);

        bulk_push_end();

        m_limit_extract = false;
    }

    //! \}

protected:
    //! Flushes all elements of the insertion heaps which are greater
    //! or equal to a given limit.
    //! \param limit limit value
    void flush_insertion_heaps_with_limit(const value_type& limit)
    {
        // perform extract for all items < L into back of insertion_heap
        std::vector<unsigned_type> back_size(m_num_insertion_heaps);

//#if STXXL_PARALLEL
//#pragma omp parallel for
//#endif
        for (size_t p = 0; p < m_num_insertion_heaps; ++p)
        {
            heap_type& insheap = m_proc[p]->insertion_heap;

            typename heap_type::iterator back = insheap.end();

            while (back != insheap.begin() &&
                   m_compare(limit, insheap[0]))
            {
                // while top < L, perform pop_heap: put top to back and
                // siftDown new items (shortens heap by one)
                std::pop_heap(insheap.begin(), back, m_compare);
                --back;
            }

            // range insheap.begin() + back to insheap.end() is < L, rest >= L.

            for (typename heap_type::const_iterator it = insheap.begin();
                 it != insheap.end(); ++it)
            {
                if (it < back)
                    assert(!m_compare(limit, *it));
                else
                    assert(m_compare(limit, *it));
            }

            back_size[p] = insheap.end() - back;
        }

        // put items from insertion heaps into an internal array
        unsigned_type back_sum = std::accumulate(
            back_size.begin(), back_size.end(), unsigned_type(0));

        STXXL_DEBUG("flush_insertion_heaps_with_limit(): back_sum = " << back_sum);

        if (back_sum)
        {
            // test that enough RAM is available for remaining items
            flush_ia_ea_until_memory_free(back_sum * sizeof(value_type));

            std::vector<value_type> values(back_sum);

            // copy items into values vector
            typename std::vector<value_type>::iterator vi = values.begin();
            for (size_t p = 0; p < m_num_insertion_heaps; ++p)
            {
                heap_type& insheap = m_proc[p]->insertion_heap;

                std::copy(insheap.end() - back_size[p], insheap.end(), vi);
                vi += back_size[p];
                insheap.resize(insheap.size() - back_size[p]);

                if (insheap.empty())
                    m_minima.deactivate_heap(p);
                else
                    m_minima.update_heap(p);
            }

            potentially_parallel::sort(values.begin(), values.end(), m_inv_compare);

            add_as_internal_array(values);
            m_heaps_size -= back_sum;
        }
    }

public:
    /*!
     * Merges all external arrays and all internal arrays into one external array.
     * Public for benchmark purposes.
     */
    void merge_external_arrays()
    {
        STXXL_ERRMSG("Merging external arrays. This should not happen."
                     << " You should adjust memory assignment and/or external array level size.");
        check_external_level(0, true);
        STXXL_DEBUG("Merging all external arrays done.");

        resize_read_pool();

        // Rebuild hint tree completely as the hint sequence may have changed.
        if (!m_in_bulk_push)
            rebuild_hint_tree();
        else
            assert(m_external_arrays.size() - 1 >= m_bulk_first_delayed_external_array);

        check_invariants();
    }

    //! Free up memory by flushing internal arrays and combining external
    //! arrays until enough bytes are free.
    void flush_ia_ea_until_memory_free(internal_size_type mem_free)
    {
        if (m_mem_left >= mem_free) return;

        if (m_internal_size > 0) {
            flush_internal_arrays();
        }
        else {
            merge_external_arrays();
        }

        assert(m_mem_left >= mem_free);
    }

    //! Automatically resize the read/prefetch buffer pool depending on number
    //! of external arrays.
    void resize_read_pool()
    {
        unsigned_type new_num_read_blocks =
            m_num_read_blocks_per_ea * m_external_arrays.size();

        STXXL_DEBUG("resize_read_pool:" <<
                    " m_num_read_blocks=" << m_num_read_blocks <<
                    " ea_size=" << m_external_arrays.size() <<
                    " m_num_read_blocks_per_ea=" << m_num_read_blocks_per_ea <<
                    " new_num_read_blocks=" << new_num_read_blocks <<
                    " free_size_prefetch=" << m_pool.free_size_prefetch() <<
                    " m_num_hinted_blocks=" << m_num_hinted_blocks <<
                    " m_num_used_read_blocks=" << m_num_used_read_blocks);

        // add new blocks
        if (new_num_read_blocks > m_num_read_blocks)
        {
            unsigned_type mem_needed =
                (new_num_read_blocks - m_num_read_blocks) * block_size;

            // -tb: this may recursively call this function!
            //flush_ia_ea_until_memory_free(mem_needed);
            STXXL_ASSERT(m_mem_left >= mem_needed);

            while (new_num_read_blocks > m_num_read_blocks) {
                block_type* new_block = new block_type();
                m_pool.add_prefetch(new_block);
                ++m_num_read_blocks;
            }

            m_mem_left -= mem_needed;
        }

        // steal extra blocks (as many as possible)
        if (new_num_read_blocks < m_num_read_blocks)
        {
            while (new_num_read_blocks < m_num_read_blocks &&
                   m_pool.free_size_prefetch() > 0)
            {
                block_type* del_block = m_pool.steal_prefetch();
                delete del_block;
                --m_num_read_blocks;
                m_mem_left += block_size;
            }

            if (new_num_read_blocks < m_num_read_blocks)
                STXXL_ERRMSG("WARNING: could not immediately reduce read/prefetch pool!");
        }
    }

    //! Rebuild hint tree completely as the hint sequence may have changed, and
    //! re-hint the correct block sequence.
    void rebuild_hint_tree()
    {
        m_stats.hint_time.start();

        // prepare rehinting sequence: reset hint begin pointer
        for (size_t i = 0; i < m_external_arrays.size(); ++i)
            m_external_arrays[i].rebuild_hints_prepare();

        // rebuild hint tree with first elements
        for (size_t i = 0; i < m_external_arrays.size(); ++i)
        {
            if (m_external_arrays[i].has_unhinted_em_data()) {
                m_hint_tree.activate_without_replay(i);
            }
            else {
                m_hint_tree.deactivate_without_replay(i);
            }
        }
        m_hint_tree.rebuild();

        // virtually release all hints
        unsigned_type free_prefetch_blocks =
            m_pool.free_size_prefetch() + m_num_hinted_blocks;
        m_num_hinted_blocks = 0;

        int gmin_index;
        while (free_prefetch_blocks > 0 &&
               (gmin_index = m_hint_tree.top()) >= 0)
        {
            assert((size_t)gmin_index < m_external_arrays.size());

            STXXL_DEBUG("Give pre-hint in EA[" << gmin_index << "] min " <<
                        m_external_arrays[gmin_index].get_next_hintable_min());

            m_external_arrays[gmin_index].rebuild_hints_prehint_next_block();
            --free_prefetch_blocks;
            ++m_num_hinted_blocks;

            if (m_external_arrays[gmin_index].has_unhinted_em_data()) {
                m_hint_tree.replay_on_change(gmin_index);
            }
            else {
                m_hint_tree.deactivate_player(gmin_index);
            }
        }

        // invalidate all hinted blocks no longer needed
        for (size_t i = 0; i < m_external_arrays.size(); ++i)
            m_external_arrays[i].rebuild_hints_cancel();

        // perform real hinting on pre-hinted blocks
        for (size_t i = 0; i < m_external_arrays.size(); ++i)
            m_external_arrays[i].rebuild_hints_finish();

        assert(free_prefetch_blocks == m_pool.free_size_prefetch());

        m_stats.hint_time.stop();
    }

    //! Updates the prefetch prediction tree afer a remove_items(), which frees
    //! up blocks.
    //! \param ea_index index of the external array in question
    inline void update_hint_tree(size_t ea_index)
    {
        m_stats.hint_time.start();
        if (m_external_arrays[ea_index].has_unhinted_em_data()) {
            m_hint_tree.replay_on_change(ea_index);
        }
        else {
            m_hint_tree.deactivate_player(ea_index);
        }
        m_stats.hint_time.stop();
    }

    //! Updates the external min tree afer a remove() or a
    //! wait_next_blocks() call.
    //! \param ea_index index of the external array in question
    inline void update_external_min_tree(size_t ea_index)
    {
        if (m_external_arrays[ea_index].has_em_data()) {
            m_external_min_tree.replay_on_change(ea_index);
        }
        else {
            m_external_min_tree.deactivate_player(ea_index);
        }
    }

    //! Hints EA blocks which will be needed soon. Hints at most
    //! m_num_prefetchers blocks globally.
    inline void hint_external_arrays()
    {
        m_stats.hint_time.start();

        STXXL_DEBUG("hint_external_arrays()"
                    " for free_size_prefetch=" << m_pool.free_size_prefetch());

        int gmin_index;
        while (m_pool.free_size_prefetch() > 0 &&
               (gmin_index = m_hint_tree.top()) >= 0)
        {
            assert((size_t)gmin_index < m_external_arrays.size());

            STXXL_DEBUG("Give hint in EA[" << gmin_index << "]");
            m_external_arrays[gmin_index].hint_next_block();
            ++m_num_hinted_blocks;

            if (m_external_arrays[gmin_index].has_unhinted_em_data()) {
                m_hint_tree.replay_on_change(gmin_index);
            }
            else {
                m_hint_tree.deactivate_player(gmin_index);
            }
        }

        m_stats.hint_time.stop();
    }

    //! Print statistics.
    void print_stats() const
    {
        STXXL_VARDUMP(c_merge_sorted_heaps);
        STXXL_VARDUMP(c_limit_extract_buffer);
        STXXL_VARDUMP(c_single_insert_limit);

        if (c_limit_extract_buffer) {
            STXXL_VARDUMP(m_extract_buffer_limit);
            STXXL_MEMDUMP(m_extract_buffer_limit * sizeof(value_type));
        }

#if STXXL_PARALLEL
        STXXL_VARDUMP(omp_get_max_threads());
#endif

        STXXL_MEMDUMP(m_mem_for_heaps);
        STXXL_MEMDUMP(m_mem_left);

        //if (num_extract_buffer_refills > 0) {
        //    STXXL_VARDUMP(total_extract_buffer_size / num_extract_buffer_refills);
        //    STXXL_MEMDUMP(total_extract_buffer_size / num_extract_buffer_refills * sizeof(value_type));
        //}

        STXXL_MSG(m_stats);
        m_minima.print_stats();
    }

protected:
    //! Calculates the sequences vector needed by the multiway merger,
    //! considering inaccessible data from external arrays.
    //! The sizes vector stores the size of each sequence.
    //! \param reuse_previous_lower_bounds Reuse upper bounds from previous runs.
    //!             sequences[i].second must be valid upper bound iterator from a previous run!
    //! \returns the index of the external array which is limiting factor
    //!             or m_external_arrays.size() if not limited.
    size_t calculate_merge_sequences(std::vector<size_type>& sizes,
                                     std::vector<iterator_pair_type>& sequences,
                                     bool reuse_previous_lower_bounds = false)
    {
        STXXL_DEBUG("calculate merge sequences");

        static const bool debug = false;

        const size_type eas = m_external_arrays.size();
        const size_type ias = m_internal_arrays.size();

        assert(sizes.size() == eas + ias);
        assert(sequences.size() == eas + ias);

        /*
         * determine minimum of each first block
         */

        int gmin_index = m_external_min_tree.top();
        bool needs_limit = (gmin_index >= 0) ? true : false;

// test correctness of external block min tree
#ifdef STXXL_DEBUG_ASSERTIONS

        bool test_needs_limit = false;
        int test_gmin_index = 0;
        value_type test_gmin_value;

        m_stats.refill_minmax_time.start();
        for (size_type i = 0; i < eas; ++i) {
            if (m_external_arrays[i].has_em_data()) {
                const value_type& min_value =
                    m_external_arrays[i].get_next_block_min();

                if (!test_needs_limit) {
                    test_needs_limit = true;
                    test_gmin_value = min_value;
                    test_gmin_index = i;
                }
                else {
                    STXXL_DEBUG("min[" << i << "]: " << min_value <<
                                " test: " << test_gmin_value <<
                                ": " << m_inv_compare(min_value, test_gmin_value));
                    if (m_inv_compare(min_value, test_gmin_value)) {
                        test_gmin_value = min_value;
                        test_gmin_index = i;
                    }
                }
            }
        }
        m_stats.refill_minmax_time.stop();

        STXXL_ASSERT(needs_limit == test_needs_limit);
        STXXL_ASSERT(!needs_limit || gmin_index == test_gmin_index);

#endif

        /*
         * calculate size and create sequences to merge
         */

#if STXXL_PARALLEL
//        #pragma omp parallel for if(eas + ias > m_num_insertion_heaps)
#endif
        for (size_type i = 0; i < eas + ias; ++i) {
            iterator begin, end;

            if (i < eas) {
                begin = m_external_arrays[i].begin();
                end = m_external_arrays[i].end();
            }
            else {
                size_type j = i - eas;
                begin = m_internal_arrays[j].begin();
                end = m_internal_arrays[j].end();
            }

            if (needs_limit) {
                const value_type& gmin_value =
                    m_external_arrays[gmin_index].get_next_block_min();

                // remove timer if parallel
                //stats.refill_lower_bound_time.start();
                if (reuse_previous_lower_bounds) {
                    // Be careful that sequences[i].second is really valid and
                    // set by a previous calculate_merge_sequences() run!
                    end = std::lower_bound(sequences[i].second, end,
                                           gmin_value, m_inv_compare);
                }
                else
                {
                    end = std::lower_bound(begin, end,
                                           gmin_value, m_inv_compare);
                }
                //stats.refill_lower_bound_time.stop();
            }

            sizes[i] = std::distance(begin, end);
            sequences[i] = std::make_pair(begin, end);

            STXXL_DEBUG("sequence[" << i << "] " << (i < eas ? "ea " : "ia ") <<
                        begin << " - " << end <<
                        " size " << sizes[i] <<
                        (needs_limit ? " with ub limit" : ""));
        }

        if (needs_limit) {
            STXXL_DEBUG("return with needs_limit: gmin_index=" << gmin_index);
            return gmin_index;
        }
        else {
            STXXL_DEBUG("return with needs_limit: eas=" << eas);
            return eas;
        }
    }

protected:
    //! Convert extract buffer into a new internal array.
    void convert_eb_into_ia(bool do_not_flush = false)
    {
        if (m_extract_buffer_size == 0) return;

        STXXL_DEBUG("convert_eb_into_ia");

        // tb: if in limit sequence and the EB gets flushed out to EM, then we
        // have to re-merge items into the EB instead of returning the
        // sentinel.
        m_limit_has_full_range = false;

        // TODO: memory is NOT allocated, but extract buffer is currently not
        // counted
        if (!do_not_flush)
            flush_ia_ea_until_memory_free(
                internal_array_type::int_memory(m_extract_buffer.size())
                );

        if (m_extract_buffer_size == 0) return;

        // first deactivate extract buffer to replay tree for new IA.
        m_minima.deactivate_extract_buffer();

        // add eb as internal array with current index
        add_as_internal_array(m_extract_buffer, m_extract_buffer_index);

        m_extract_buffer_index = 0;
        m_extract_buffer_size = 0;
    }

    //! Refills the extract buffer from the external arrays.
    //! \param minimum_size requested minimum size of the resulting extract buffer.
    //!         Prints a warning if there is not enough data to reach this size.
    //! \param maximum_size maximum size of the extract buffer. Using
    //!         m_extract_buffer_limit if set to 0.
    inline void refill_extract_buffer(size_t minimum_size = 0,
                                      size_t maximum_size = 0)
    {
        STXXL_DEBUG("refilling extract buffer" <<
                    " ia_size=" << m_internal_arrays.size() <<
                    " ea_size=" << m_external_arrays.size());

        if (maximum_size == 0)
            maximum_size = m_extract_buffer_limit;

        check_invariants();

        assert(extract_buffer_empty());
        m_extract_buffer_index = 0;

        cleanup_external_arrays();

        size_type ias, eas = m_external_arrays.size();

        m_minima.clear_internal_arrays();
        cleanup_internal_arrays();
        ias = m_internal_arrays.size();

        if (eas == 0 && ias == 0) {
            m_extract_buffer.resize(0);
            m_minima.deactivate_extract_buffer();
            return;
        }

        m_stats.num_extract_buffer_refills++;
        m_stats.refill_extract_buffer_time.start();
        m_stats.refill_time_before_merge.start();

        std::vector<size_type> sizes(eas + ias);
        std::vector<iterator_pair_type> sequences(eas + ias);
        size_type output_size = 0;

        if (minimum_size > 0) {
            size_t limiting_ea_index = eas + 1;
            bool reuse_lower_bounds = false;
            while (output_size < minimum_size)
            {
                STXXL_DEBUG("refill: request more data," <<
                            " output_size=" << output_size <<
                            " minimum_size=" << minimum_size <<
                            " limiting_ea_index=" << limiting_ea_index);

                if (limiting_ea_index < eas) {
                    if (m_external_arrays[limiting_ea_index].num_hinted_blocks() == 0)
                        break;

                    wait_next_ea_blocks(limiting_ea_index);
                    reuse_lower_bounds = true;
                }
                else if (limiting_ea_index == eas) {
                    // no more unaccessible EM data
                    STXXL_MSG("Warning: refill_extract_buffer(n): "
                              "minimum_size > # mergeable elements!");
                    break;
                }

                limiting_ea_index = calculate_merge_sequences(
                    sizes, sequences, reuse_lower_bounds);

                output_size = std::accumulate(sizes.begin(), sizes.end(), 0);
            }
        }
        else {
            calculate_merge_sequences(sizes, sequences);
            output_size = std::accumulate(sizes.begin(), sizes.end(), 0);
        }

        if (c_limit_extract_buffer) {
            output_size = std::min<size_t>(output_size, maximum_size);
        }

        m_stats.max_extract_buffer_size.set_max(output_size);
        m_stats.total_extract_buffer_size += output_size;

        assert(output_size > 0);
        m_extract_buffer.resize(output_size);
        m_extract_buffer_size = output_size;

        m_stats.refill_time_before_merge.stop();
        m_stats.refill_merge_time.start();

        potentially_parallel::multiway_merge(
            sequences.begin(), sequences.end(),
            m_extract_buffer.begin(), output_size, m_inv_compare);

        m_stats.refill_merge_time.stop();
        m_stats.refill_time_after_merge.start();

        advance_arrays(sequences, sizes, eas, ias);

        m_minima.update_extract_buffer();

        m_stats.refill_time_after_merge.stop();
        m_stats.refill_extract_buffer_time.stop();

        check_invariants();
    }

    //! Requests more EM data from a given EA and updates
    //! the winner trees and hints accordingly.
    inline void wait_next_ea_blocks(unsigned_type ea_index)
    {
        unsigned_type used_blocks =
            m_external_arrays[ea_index].wait_next_blocks();

        m_num_hinted_blocks -= used_blocks;
        m_num_used_read_blocks += used_blocks;

        update_external_min_tree(ea_index);
    }

    // Removes empty arrays and updates the winner trees accordingly
    inline void advance_arrays(std::vector<iterator_pair_type>& sequences,
                               std::vector<size_type>& sizes,
                               size_t eas, size_t ias)
    {
        unsigned_type total_freed_blocks = 0;

        for (size_type i = 0; i < eas + ias; ++i) {
            // dist represents the number of elements that haven't been merged
            size_type dist = std::distance(sequences[i].first,
                                           sequences[i].second);
            const size_t diff = sizes[i] - dist;
            if (diff == 0) continue;

            if (i < eas) {
                // remove items and free blocks in RAM.
                unsigned_type freed_blocks =
                    m_external_arrays[i].remove_items(diff);

                m_num_used_read_blocks -= freed_blocks;
                total_freed_blocks += freed_blocks;

                // correct item count.
                assert(m_external_size >= diff);
                m_external_size -= diff;
            }
            else {
                size_type j = i - eas;
                m_internal_arrays[j].inc_min(diff);
                assert(m_internal_size >= diff);
                m_internal_size -= diff;
            }
        }

        // remove empty arrays - important for the next round (may also reduce
        // number of prefetch buffers, so must be before hinting).
        cleanup_external_arrays();

        // prefetch new blocks from EAs using freed blocks
        if (total_freed_blocks)
            hint_external_arrays();

        m_stats.num_new_external_arrays = 0;
        cleanup_internal_arrays();
    }

    //! Flushes the insertions heap p into an internal array.
    inline void flush_insertion_heap(unsigned_type p)
    {
        assert(m_proc[p]->insertion_heap.size() != 0);

        heap_type& insheap = m_proc[p]->insertion_heap;
        size_t size = insheap.size();

        STXXL_DEBUG0(
            "Flushing insertion heap array p=" << p <<
            " size=" << insheap.size() <<
            " capacity=" << insheap.capacity() <<
            " int_memory=" << internal_array_type::int_memory(insheap.size()) <<
            " mem_left=" << m_mem_left);

        m_stats.num_insertion_heap_flushes++;
        stats_timer flush_time(true); // separate timer due to parallel sorting

        // sort locally, independent of others
        std::sort(insheap.begin(), insheap.end(), m_inv_compare);

#if STXXL_PARALLEL
#pragma omp critical (stxxl_flush_insertion_heap)
#endif
        {
            // test that enough RAM is available for merged internal array:
            // otherwise flush the existing internal arrays out to disk.
            flush_ia_ea_until_memory_free(
                internal_array_type::int_memory(insheap.size()));

            // invalidate player in minima tree (before adding the IA to tree)
            m_minima.deactivate_heap(p);

            // insheap is empty afterwards, as vector was swapped into new_array
            add_as_internal_array(insheap);

            // reserve new insertion heap
            insheap.reserve(m_insertion_heap_capacity);
            assert(insheap.capacity() * sizeof(value_type)
                   == insertion_heap_int_memory());

            // update item counts
#if STXXL_PARALLEL
#pragma omp atomic
#endif
            m_heaps_size -= size;
        }

        m_stats.insertion_heap_flush_time += flush_time;
    }

    //! Flushes all insertions heaps into an internal array.
    inline void flush_insertion_heaps()
    {
        size_type max_mem_needed;

        if (c_merge_sorted_heaps) {
            max_mem_needed = m_mem_for_heaps;
        }
        else {
            max_mem_needed = insertion_heap_int_memory();
        }

        // test that enough RAM is available for merged internal array:
        // otherwise flush the existing internal arrays out to disk.
        flush_ia_ea_until_memory_free(max_mem_needed);

        m_stats.num_insertion_heap_flushes++;
        m_stats.insertion_heap_flush_time.start();

        size_type size = m_heaps_size;
        size_type int_memory = 0;
        assert(size > 0);
        std::vector<std::pair<value_iterator, value_iterator> > sequences(m_num_insertion_heaps);

#if STXXL_PARALLEL
        #pragma omp parallel for
#endif
        for (long i = 0; i < m_num_insertion_heaps; ++i)
        {
            heap_type& insheap = m_proc[i]->insertion_heap;

            std::sort(insheap.begin(), insheap.end(), m_inv_compare);

            if (c_merge_sorted_heaps)
                sequences[i] = std::make_pair(insheap.begin(), insheap.end());

            int_memory += insheap.capacity();
        }

        if (c_merge_sorted_heaps)
        {
            m_stats.merge_sorted_heaps_time.start();
            std::vector<value_type> merged_array(size);

            potentially_parallel::multiway_merge(
                sequences.begin(), sequences.end(),
                merged_array.begin(), size, m_inv_compare);

            m_stats.merge_sorted_heaps_time.stop();

            add_as_internal_array(merged_array);

            for (int_type i = 0; i < m_num_insertion_heaps; ++i)
            {
                m_proc[i]->insertion_heap.clear();
                m_proc[i]->insertion_heap.reserve(m_insertion_heap_capacity);
            }
            m_minima.clear_heaps();
        }
        else
        {
            for (unsigned i = 0; i < m_num_insertion_heaps; ++i)
            {
                heap_type& insheap = m_proc[i]->insertion_heap;

                if (insheap.size() == 0) continue;

                add_as_internal_array(insheap);

                // reserve new insertion heap
                insheap.reserve(m_insertion_heap_capacity);
            }

            m_minima.clear_heaps();
        }

        m_heaps_size = 0;

        m_stats.insertion_heap_flush_time.stop();

        check_invariants();
    }

    //! Flushes the internal arrays into an external array.
    void flush_internal_arrays()
    {
        STXXL_DEBUG("Flushing internal arrays" <<
                    " num_arrays=" << m_internal_arrays.size());

        m_stats.num_internal_array_flushes++;
        m_stats.internal_array_flush_time.start();

        m_minima.clear_internal_arrays();

        // also flush extract buffer items out to disk.
        convert_eb_into_ia(true);

        // clean up internal arrays that have been deleted in extract_min!
        cleanup_internal_arrays();

        size_type num_arrays = m_internal_arrays.size();
        size_type size = m_internal_size;
        size_type int_memory = 0;
        std::vector<iterator_pair_type> sequences(num_arrays);

        for (unsigned i = 0; i < num_arrays; ++i)
        {
            sequences[i] = std::make_pair(m_internal_arrays[i].begin(),
                                          m_internal_arrays[i].end());

            int_memory += m_internal_arrays[i].int_memory();
        }

        // must release more RAM in IAs than the EA takes, otherwise: merge
        // external and internal arrays!
        if (int_memory < external_array_type::int_memory(size)
            + ceil(m_num_read_blocks_per_ea) * block_size)
        {
            return merge_external_arrays();
        }

        // construct new external array

        external_array_type ea(size, &m_pool, 0);

        m_stats.max_merge_buffer_size.set_max(size);

        {
            external_array_writer_type external_array_writer(ea);

            potentially_parallel::multiway_merge(
                sequences.begin(), sequences.end(),
                external_array_writer.begin(), size, m_inv_compare);
        }

        STXXL_DEBUG("Merge done of new ea " << &ea);

        m_external_arrays.swap_back(ea);

        m_internal_size = 0;
        m_external_size += size;

        // register EA in min tree
        // important for check_external_level()!
        m_external_min_tree.activate_without_replay(m_external_arrays.size() - 1);
        update_external_min_tree(m_external_arrays.size() - 1);

        // register EA in hint tree
        m_hint_tree.activate_without_replay(m_external_arrays.size() - 1);
        if (!m_in_bulk_push)
            update_hint_tree(m_external_arrays.size() - 1);
        // else: done in bulk_push_end() -> rebuild_hint_tree()

        m_internal_arrays.clear();
        m_stats.num_new_internal_arrays = 0;
        cleanup_internal_arrays();

        // TODO: is this necessary? See cleanup_internal_arrays().
        for (size_t i = 0; i < c_max_internal_levels; ++i)
            m_internal_levels[i] = 0;

        m_mem_left += int_memory;
        m_mem_left -= m_external_arrays.back().int_memory();

        m_stats.max_num_external_arrays.set_max(m_external_arrays.size());
        m_stats.internal_array_flush_time.stop();

        // update EA level and potentially merge
        ++m_external_levels[0];
        check_external_level(0);

        resize_read_pool();
        // Rebuild hint tree completely as the hint sequence may have changed.
        if (!m_in_bulk_push)
            rebuild_hint_tree();
        else
            assert(m_external_arrays.size() - 1 >= m_bulk_first_delayed_external_array);

        check_invariants();
    }

    // Compares the largest accessible value of two external arrays.
    struct s_min_tree_comparator {
        const external_arrays_type& m_eas;
        const std::vector<unsigned_type>& m_indices;
        const inv_compare_type& m_compare;

        s_min_tree_comparator(const external_arrays_type& eas,
                              const inv_compare_type& compare,
                              const std::vector<unsigned_type>& indices)
            : m_eas(eas), m_indices(indices), m_compare(compare) { }

        bool operator () (const size_t& a, const size_t& b) const
        {
            return m_compare(m_eas[m_indices[a]].get_next_hintable_min(),
                             m_eas[m_indices[b]].get_next_hintable_min());
        }
    };

    //! Merges external arrays if there are too many external arrays on
    //! the same level.
    void check_external_level(unsigned_type level, bool force_merge_all = false)
    {
        if (!force_merge_all)
            STXXL_DEBUG("Checking external level " << level);

        // return if EA level is not full
        if (m_external_levels[level] < c_max_external_level_size && !force_merge_all)
            return;

        unsigned_type level_size = 0;
        size_type int_memory = 0;
        std::vector<unsigned_type> ea_index;

        for (unsigned_type i = 0; i < m_external_arrays.size(); ++i)
        {
            if (m_external_arrays[i].level() != level && !force_merge_all) continue;
            if (m_external_arrays[i].empty()) continue;

            level_size += m_external_arrays[i].size();
            int_memory += m_external_arrays[i].int_memory();
            ea_index.push_back(i);
        }

        // return if there is not enough RAM for the new array.
        // TODO: force_merge_all==true is for freeing memory. Breaking here is not
        // helpful in this case. But one should maybe reserve some space in advance.
        if (m_mem_left < external_array_type::int_memory(level_size) && !force_merge_all)
            return;
        m_mem_left -= external_array_type::int_memory(level_size);

        STXXL_ASSERT(force_merge_all || c_max_external_level_size == ea_index.size());
        unsigned_type num_arrays_to_merge = ea_index.size();

        STXXL_DEBUG("merging external arrays" <<
                    " level=" << level <<
                    " level_size=" << level_size <<
                    " sequences=" << num_arrays_to_merge <<
                    " force_merge_all=" << force_merge_all);

        // if force_merge_all: create array in highest level to avoid merging
        // of such a large EA.
        unsigned_type new_level = force_merge_all ? c_max_external_levels - 1 : level + 1;

        // construct new external array
        external_array_type ea(level_size, &m_pool, new_level);
        {
            external_array_writer_type external_array_writer(ea);
            typename external_array_writer_type::iterator out_iter
                = external_array_writer.begin();

            // === build minima_tree over the level's arrays ===

            s_min_tree_comparator min_tree_comparator(m_external_arrays,
                                                      m_inv_compare, ea_index);

            winner_tree<s_min_tree_comparator> min_tree(num_arrays_to_merge,
                                                        min_tree_comparator);

            // =================================================

            int_type num_arrays_done = 0;

            while (num_arrays_to_merge != num_arrays_done)
            {
                STXXL_DEBUG("num_arrays_done = " << num_arrays_done);

                // === build hints ===

                for (int_type i = 0; i < num_arrays_to_merge; ++i) {
                    if (m_external_arrays[ea_index[i]].has_unhinted_em_data()) {
                        min_tree.activate_without_replay(i);
                    }
                    else {
                        min_tree.deactivate_without_replay(i);
                    }
                }

                min_tree.rebuild();

                // === fill available memory with read blocks ===
                while (m_mem_left >= block_size) {
                    block_type* new_block = new block_type();
                    m_pool.add_prefetch(new_block);
                    ++m_num_read_blocks;
                    m_mem_left -= block_size;
                }
                // ==============================================

                // cleanup hints (all arrays, not only the ones to merge)
                for (unsigned_type i = 0; i < m_external_arrays.size(); ++i) {
                    m_external_arrays[i].rebuild_hints_prepare();
                }

                // virtually release all hints
                unsigned_type free_prefetch_blocks =
                    m_pool.free_size_prefetch() + m_num_hinted_blocks;
                m_num_hinted_blocks = 0;

                int gmin_index_index; // index in ea_index
                while (free_prefetch_blocks > 0 &&
                       (gmin_index_index = min_tree.top()) >= 0)
                {
                    const unsigned_type gmin_index = ea_index[gmin_index_index];
                    assert(gmin_index < m_external_arrays.size());

                    STXXL_DEBUG0("check_external_level():Give pre-hint in EA[" << gmin_index << "] min " <<
                                 m_external_arrays[gmin_index].get_next_hintable_min());

                    m_external_arrays[gmin_index].rebuild_hints_prehint_next_block();
                    --free_prefetch_blocks;
                    ++m_num_hinted_blocks;

                    if (m_external_arrays[gmin_index].has_unhinted_em_data()) {
                        min_tree.replay_on_change(gmin_index_index);
                    }
                    else {
                        min_tree.deactivate_player(gmin_index_index);
                    }
                }

                // invalidate all hinted blocks no longer needed
                // (all arrays, not only the ones to merge)
                for (size_t i = 0; i < m_external_arrays.size(); ++i)
                    m_external_arrays[i].rebuild_hints_cancel();

                // perform real hinting on pre-hinted blocks
                // (all arrays, not only the ones to merge)
                for (size_t i = 0; i < m_external_arrays.size(); ++i)
                    m_external_arrays[i].rebuild_hints_finish();

                assert(free_prefetch_blocks == m_pool.free_size_prefetch());

                // ================================ end build hints ======

                // === wait for data ===
                for (size_type i = 0; i < num_arrays_to_merge; ++i) {
                    const unsigned_type index = ea_index[i];

                    unsigned_type used_blocks =
                        m_external_arrays[index].wait_all_hinted_blocks();

                    m_num_hinted_blocks -= used_blocks;
                    m_num_used_read_blocks += used_blocks;
                }
                // =====================

                // === build sequences ===
                std::vector<iterator_pair_type> sequences(num_arrays_to_merge);
                std::vector<size_type> sizes(num_arrays_to_merge);

                gmin_index_index = min_tree.top();
                bool needs_limit = (gmin_index_index >= 0) ? true : false;

                for (size_type i = 0; i < num_arrays_to_merge; ++i) {
                    const unsigned_type index = ea_index[i];
                    iterator begin = m_external_arrays[index].begin();
                    iterator end = m_external_arrays[index].end();

                    if (needs_limit) {
                        const unsigned_type gmin_index = ea_index[gmin_index_index];
                        const value_type& gmin_value =
                            m_external_arrays[gmin_index].get_next_block_min();

                        end = std::lower_bound(begin, end,
                                               gmin_value, m_inv_compare);
                    }

                    sizes[i] = std::distance(begin, end);
                    sequences[i] = std::make_pair(begin, end);

                    STXXL_DEBUG("sequence[" << i << "] ea " <<
                                begin << " - " << end <<
                                " size " << sizes[i] <<
                                (needs_limit ? " with ub limit" : ""));
                }
                // ==========================================

                // === merge ===

                size_type output_size = std::accumulate(sizes.begin(), sizes.end(), 0);

                out_iter = potentially_parallel::multiway_merge(
                    sequences.begin(), sequences.end(),
                    out_iter, output_size, m_inv_compare);

                for (unsigned_type i = 0; i < num_arrays_to_merge; ++i) {
                    const unsigned_type index = ea_index[i];

                    if (!m_external_arrays[index].empty()) {
                        // remove items and free blocks in RAM.
                        unsigned_type freed_blocks =
                            m_external_arrays[index].remove_items(sizes[i]);

                        m_num_used_read_blocks -= freed_blocks;

                        if (m_external_arrays[index].empty())
                            ++num_arrays_done;
                    }
                }

                // reset read buffer
                resize_read_pool();

                // cannot call clear_external_arrays() here, since it
                // corrupts ea_index.
            }

            if (m_in_bulk_push)
                m_bulk_first_delayed_external_array = 0; // TODO: workaround
        } // destroy external_array_writer

        // clean up now empty arrays
        cleanup_external_arrays();

        m_external_arrays.swap_back(ea);
        ++m_external_levels[new_level];

        // register EA in min tree
        m_external_min_tree.activate_without_replay(m_external_arrays.size() - 1);
        update_external_min_tree(m_external_arrays.size() - 1);

        // register EA in hint tree
        m_hint_tree.activate_without_replay(m_external_arrays.size() - 1);
        if (!m_in_bulk_push)
            update_hint_tree(m_external_arrays.size() - 1);
        // else: done in bulk_push_end() -> rebuild_hint_tree()

        STXXL_DEBUG("Merge done of new ea " << &ea);

        if (!force_merge_all)
            check_external_level(level + 1);

        check_invariants();
    }

    //! Add new internal array, which requires that values are sorted!
    //! automatically decreases m_mem_left! also merges internal arrays if
    //! there are too many internal arrays on the same level.
    void add_as_internal_array(std::vector<value_type>& values,
                               unsigned_type used = 0,
                               unsigned_type level = 0)
    {
        const size_t size = values.size();
        const size_t capacity = values.capacity();
        assert(size > used); // at least one element

        internal_array_type new_array(values, used, level);
        STXXL_ASSERT(new_array.int_memory() ==
                     internal_array_type::int_memory(capacity));
        m_internal_arrays.swap_back(new_array);

        if (!extract_buffer_empty()) {
            m_stats.num_new_internal_arrays++;
            m_stats.max_num_new_internal_arrays.set_max(
                m_stats.num_new_internal_arrays);
            m_minima.add_internal_array(
                static_cast<unsigned>(m_internal_arrays.size()) - 1
                );
        }

        m_internal_size += size - used;
        m_mem_left -= internal_array_type::int_memory(capacity);

        STXXL_CHECK(level < c_max_internal_levels &&
                    "Internal array level is larger than anything possible "
                    "in this universe. Increase the size of m_internal_levels");

        ++m_internal_levels[level];

        m_stats.max_num_internal_arrays.set_max(m_internal_arrays.size());

        // if IA level is too large ...
        if (m_internal_levels[level] < c_max_internal_level_size) return;

        unsigned_type level_size = 0;
        size_type int_memory = 0;
        std::vector<iterator_pair_type> sequences;
        std::vector<unsigned_type> ia_index;

        for (unsigned_type i = 0; i < m_internal_arrays.size(); ++i)
        {
            if (m_internal_arrays[i].level() != level) continue;
            if (m_internal_arrays[i].empty()) continue;

            level_size += m_internal_arrays[i].size();
            int_memory += m_internal_arrays[i].int_memory();
            sequences.push_back(std::make_pair(m_internal_arrays[i].begin(),
                                               m_internal_arrays[i].end()));
            ia_index.push_back(i);
        }

        // AND there is enough RAM to merge it (without flushing out to EA).
        if (m_mem_left < internal_array_type::int_memory(level_size)) return;

        // must free up more memory than the new array needs.
        STXXL_ASSERT(int_memory >= internal_array_type::int_memory(level_size));

        STXXL_DEBUG("merging internal arrays" <<
                    " level=" << level <<
                    " level_size=" << level_size <<
                    " sequences=" << sequences.size());

        std::vector<value_type> merged_array(level_size);

        potentially_parallel::multiway_merge(
            sequences.begin(), sequences.end(),
            merged_array.begin(), level_size, m_inv_compare);

        // release memory of old internal arrays immediately
        for (unsigned_type i = 0; i < ia_index.size(); ++i)
        {
            unsigned_type ia = ia_index[i];
            m_internal_arrays[ia].make_empty();
            // this is done in cleanup_internal_arrays()...
            //if (ia < m_minima.ia_slots())
            //    m_minima.deactivate_internal_array(ia);
        }

        cleanup_internal_arrays();

        // in add_as_internal_array the level_size is re-added!
        m_internal_size -= level_size;

        // add as new internal array at next level (and maybe recursively merge)
        add_as_internal_array(merged_array, 0, level + 1);
    }

    /*!
     * Sorts the values from values and writes them into an internal array.
     * Don't use the value vector afterwards!
     *
     * \param values the vector to sort and store
     */
    void flush_array_internal(std::vector<value_type>& values)
    {
        potentially_parallel::sort(values.begin(), values.end(), m_inv_compare);

        // flush until enough memory for new array
        flush_ia_ea_until_memory_free(
            internal_array_type::int_memory(values.size())
            );

        add_as_internal_array(values);
    }

    //! Struct of all statistical counters and timers.  Turn on/off statistics
    //! using the stats_counter and stats_timer typedefs.
    struct stats_type
    {
        //! Largest number of elements in the extract buffer at the same time
        stats_counter max_extract_buffer_size;

        //! Sum of the sizes of each extract buffer refill. Used for average
        //! size.
        stats_counter total_extract_buffer_size;

        //! Largest number of elements in the merge buffer when running
        //! flush_internal_arrays()
        stats_counter max_merge_buffer_size;

        //! Total number of extracts
        stats_counter num_extracts;

        //! Number of refill_extract_buffer() calls
        stats_counter num_extract_buffer_refills;

        //! Number of flush_insertion_heaps() calls
        stats_counter num_insertion_heap_flushes;

        //! Number of flush_directly_to_hd() calls
        stats_counter num_direct_flushes;

        //! Number of flush_internal_arrays() calls
        stats_counter num_internal_array_flushes;

        //! Number of merge_external_arrays() calls
        stats_counter num_external_array_merges;

        //! Largest number of internal arrays at the same time
        stats_counter max_num_internal_arrays;

        //! Largest number of external arrays at the same time
        stats_counter max_num_external_arrays;

        //! Temporary number of new external arrays at the same time (which
        //! were created while the extract buffer hadn't been empty)
        stats_counter num_new_external_arrays;

        //! Largest number of new external arrays at the same time (which were
        //! created while the extract buffer hadn't been empty)
        stats_counter max_num_new_external_arrays;

        //! Temporary number of new internal arrays at the same time (which
        //! were created while the extract buffer hadn't been empty)
        stats_counter num_new_internal_arrays;

        //! Largest number of new internal arrays at the same time (which were
        //! created while the extract buffer hadn't been empty)
        stats_counter max_num_new_internal_arrays;

        //! Total time for flush_insertion_heaps()
        stats_timer insertion_heap_flush_time;

        //! Total time for flush_directly_to_hd()
        stats_timer direct_flush_time;

        //! Total time for flush_internal_arrays()
        stats_timer internal_array_flush_time;

        //! Total time for merge_external_arrays()
        stats_timer external_array_merge_time;

        //! Total time for extract_min()
        stats_timer extract_min_time;

        //! Total time for refill_extract_buffer()
        stats_timer refill_extract_buffer_time;

        //! Total time for the merging in refill_extract_buffer()
        //! Part of refill_extract_buffer_time.
        stats_timer refill_merge_time;

        //! Total time for all things before merging in refill_extract_buffer()
        //! Part of refill_extract_buffer_time.
        stats_timer refill_time_before_merge;

        //! Total time for all things after merging in refill_extract_buffer()
        //! Part of refill_extract_buffer_time.
        stats_timer refill_time_after_merge;

        //! Total time of wait() calls in first part of
        //! refill_extract_buffer(). Part of refill_time_before_merge and
        //! refill_extract_buffer_time.
        stats_timer refill_wait_time;

        //! Total time for pop_heap() in extract_min().
        //! Part of extract_min_time.
        stats_timer pop_heap_time;

        //! Total time for merging the sorted heaps.
        //! Part of flush_insertion_heaps.
        stats_timer merge_sorted_heaps_time;

        //! Total time for std::lower_bound calls in refill_extract_buffer()
        //! Part of refill_extract_buffer_time and refill_time_before_merge.
        // stats_timer refill_lower_bound_time;

        //! Total time for std::accumulate calls in refill_extract_buffer()
        //! Part of refill_extract_buffer_time and refill_time_before_merge.
        stats_timer refill_accumulate_time;

        //! Total time for determining the smallest max value in refill_extract_buffer()
        //! Part of refill_extract_buffer_time and refill_time_before_merge.
        stats_timer refill_minmax_time;

        stats_timer hint_time;

        friend std::ostream& operator << (std::ostream& os, const stats_type& o)
        {
            return os << "max_extract_buffer_size=" << o.max_extract_buffer_size.as_memory_amount(sizeof(value_type)) << std::endl
                      << "total_extract_buffer_size=" << o.total_extract_buffer_size.as_memory_amount(sizeof(value_type)) << std::endl
                      << "max_merge_buffer_size=" << o.max_merge_buffer_size.as_memory_amount(sizeof(value_type)) << std::endl
                      << "num_extracts=" << o.num_extracts << std::endl
                      << "num_extract_buffer_refills=" << o.num_extract_buffer_refills << std::endl
                      << "num_insertion_heap_flushes=" << o.num_insertion_heap_flushes << std::endl
                      << "num_direct_flushes=" << o.num_direct_flushes << std::endl
                      << "num_internal_array_flushes=" << o.num_internal_array_flushes << std::endl
                      << "num_external_array_merges=" << o.num_external_array_merges << std::endl
                      << "max_num_internal_arrays=" << o.max_num_internal_arrays << std::endl
                      << "max_num_external_arrays=" << o.max_num_external_arrays << std::endl
                      << "num_new_external_arrays=" << o.num_new_external_arrays << std::endl
                      << "max_num_new_external_arrays=" << o.max_num_new_external_arrays << std::endl
                      << "num_new_internal_arrays=" << o.num_new_internal_arrays << std::endl
                      << "max_num_new_internal_arrays=" << o.max_num_new_internal_arrays << std::endl
                      << "insertion_heap_flush_time=" << o.insertion_heap_flush_time << std::endl
                      << "direct_flush_time=" << o.direct_flush_time << std::endl
                      << "internal_array_flush_time=" << o.internal_array_flush_time << std::endl
                      << "external_array_merge_time=" << o.external_array_merge_time << std::endl
                      << "extract_min_time=" << o.extract_min_time << std::endl
                      << "refill_extract_buffer_time=" << o.refill_extract_buffer_time << std::endl
                      << "refill_merge_time=" << o.refill_merge_time << std::endl
                      << "refill_time_before_merge=" << o.refill_time_before_merge << std::endl
                      << "refill_time_after_merge=" << o.refill_time_after_merge << std::endl
                      << "refill_wait_time=" << o.refill_wait_time << std::endl
                      << "pop_heap_time=" << o.pop_heap_time << std::endl
                      << "merge_sorted_heaps_time=" << o.merge_sorted_heaps_time << std::endl
                   // << "refill_lower_bound_time=" << o.refill_lower_bound_time << std::endl
                      << "refill_accumulate_time=" << o.refill_accumulate_time << std::endl
                      << "refill_minmax_time=" << o.refill_minmax_time << std::endl
                      << "hint_time=" << o.hint_time << std::endl;
        }
    };

    stats_type m_stats;
};

// For C++98 compatibility:
template <
    class ValueType,
    class CompareType,
    class AllocStrategy,
    uint64 BlockSize,
    uint64 DefaultMemSize,
    uint64 MaxItems
    >
const double parallel_priority_queue<ValueType, CompareType, AllocStrategy, BlockSize,
                                     DefaultMemSize, MaxItems>::c_default_extract_buffer_ram_part = 0.05;

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_PARALLEL_PRIORITY_QUEUE_HEADER
