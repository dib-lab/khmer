/***************************************************************************
 *  include/stxxl/bits/containers/priority_queue.h
 *
 *  Implements a data structure from "Peter Sanders. Fast Priority Queues for
 *  Cached Memory. ALENEX'99" for external memory.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 1999 Peter Sanders <sanders@mpi-sb.mpg.de>
 *  Copyright (C) 2003, 2004, 2007 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007-2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2007-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_PRIORITY_QUEUE_HEADER
#define STXXL_CONTAINERS_PRIORITY_QUEUE_HEADER

#include <stxxl/bits/containers/pq_helpers.h>
#include <stxxl/bits/containers/pq_mergers.h>
#include <stxxl/bits/containers/pq_int_merger.h>
#include <stxxl/bits/containers/pq_ext_merger.h>

STXXL_BEGIN_NAMESPACE

/*
   KNBufferSize1 = 32;
   KNN = 512; // length of group 1 sequences
   KNKMAX = 64;  // maximal arity
   LogKNKMAX = 6;  // ceil(log KNKMAX)
   KNLevels = 4; // overall capacity >= KNN*KNKMAX^KNLevels
 */

// internal memory consumption >= N_*(KMAX_^IntLevels_) + ext

template <
    class ValueType,
    class CompareType,
    unsigned BufferSize1_ = 32,                    // equalize procedure call overheads etc.
    unsigned N_ = 512,                             // length of group 1 sequences
    unsigned IntKMAX_ = 64,                        // maximal arity for internal mergers
    unsigned IntLevels_ = 4,                       // number of internal groups
    unsigned BlockSize_ = (2* 1024* 1024),         // external block size
    unsigned ExtKMAX_ = 64,                        // maximal arity for external mergers
    unsigned ExtLevels_ = 2,                       // number of external groups
    class AllocStr_ = STXXL_DEFAULT_ALLOC_STRATEGY
    >
struct priority_queue_config
{
    typedef ValueType value_type;
    typedef CompareType comparator_type;
    typedef AllocStr_ alloc_strategy_type;
    enum
    {
        delete_buffer_size = BufferSize1_,
        N = N_,
        IntKMAX = IntKMAX_,
        num_int_groups = IntLevels_,
        num_ext_groups = ExtLevels_,
        BlockSize = BlockSize_,
        ExtKMAX = ExtKMAX_,
        element_size = sizeof(ValueType)
    };
};

STXXL_END_NAMESPACE

namespace std {

template <class BlockType,
          class CompareType,
          unsigned Arity,
          class AllocStr>
void swap(stxxl::priority_queue_local::ext_merger<BlockType, CompareType, Arity, AllocStr>& a,
          stxxl::priority_queue_local::ext_merger<BlockType, CompareType, Arity, AllocStr>& b)
{
    a.swap(b);
}
template <class ValueType, class CompareType, unsigned KNKMAX>
void swap(stxxl::priority_queue_local::int_merger<ValueType, CompareType, KNKMAX>& a,
          stxxl::priority_queue_local::int_merger<ValueType, CompareType, KNKMAX>& b)
{
    a.swap(b);
}

} // namespace std

STXXL_BEGIN_NAMESPACE

//! External priority queue data structure \n
//! <b> Introduction </b> to priority queue container: see \ref tutorial_pqueue tutorial. \n
//! <b> Design and Internals </b> of priority queue container: see \ref design_pqueue.
template <class ConfigType>
class priority_queue : private noncopyable
{
public:
    typedef ConfigType Config;
    enum
    {
        delete_buffer_size = Config::delete_buffer_size,
        N = Config::N,
        IntKMAX = Config::IntKMAX,
        num_int_groups = Config::num_int_groups,
        num_ext_groups = Config::num_ext_groups,
        total_num_groups = Config::num_int_groups + Config::num_ext_groups,
        BlockSize = Config::BlockSize,
        ExtKMAX = Config::ExtKMAX
    };

    //! The type of object stored in the priority_queue.
    typedef typename Config::value_type value_type;
    //! Comparison object.
    typedef typename Config::comparator_type comparator_type;
    typedef typename Config::alloc_strategy_type alloc_strategy_type;
    //! An unsigned integral type (64 bit).
    typedef stxxl::uint64 size_type;
    //! Type of the block used in disk-memory transfers
    typedef typed_block<BlockSize, value_type> block_type;
    typedef read_write_pool<block_type> pool_type;

protected:
    typedef priority_queue_local::internal_priority_queue<value_type, std::vector<value_type>, comparator_type>
        insert_heap_type;

    typedef priority_queue_local::int_merger<
            value_type,
            comparator_type,
            IntKMAX> int_merger_type;

    typedef priority_queue_local::ext_merger<
            block_type,
            comparator_type,
            ExtKMAX,
            alloc_strategy_type> ext_merger_type;

    int_merger_type int_mergers[num_int_groups];
    pool_type* pool;
    bool pool_owned;
    ext_merger_type** ext_mergers;

    // one delete buffer for each tree => group buffer
    value_type group_buffers[total_num_groups][N + 1];          // tree->group_buffers->delete_buffer (extra space for sentinel)
    value_type* group_buffer_current_mins[total_num_groups];    // group_buffer_current_mins[i] is current start of group_buffers[i], end is group_buffers[i] + N

    // overall delete buffer
    value_type delete_buffer[delete_buffer_size + 1];
    value_type* delete_buffer_current_min;                      // current start of delete_buffer
    value_type* delete_buffer_end;                              // end of delete_buffer

    comparator_type cmp;

    // insert buffer
    insert_heap_type insert_heap;

    // how many groups are active
    unsigned_type num_active_groups;

    // total size not counting insert_heap and delete_buffer
    size_type size_;

private:
    void init();

    void refill_delete_buffer();
    size_type refill_group_buffer(unsigned_type k);

    unsigned_type make_space_available(unsigned_type level);
    void empty_insert_heap();

    value_type get_supremum() const { return cmp.min_value(); } //{ return group_buffers[0][KNN].key; }
    unsigned_type current_delete_buffer_size() const { return delete_buffer_end - delete_buffer_current_min; }
    unsigned_type current_group_buffer_size(unsigned_type i) const { return &(group_buffers[i][N]) - group_buffer_current_mins[i]; }

public:
    //! \name Constructors/Destructors
    //! \{

    //! Constructs external priority queue object.
    //! \param pool_ pool of blocks that will be used
    //! for data writing and prefetching for the disk<->memory transfers
    //! happening in the priority queue. Larger pool size
    //! helps to speed up operations.
    priority_queue(pool_type& pool_);

    //! Constructs external priority queue object.
    //! \param p_pool_ pool of blocks that will be used
    //! for data prefetching for the disk<->memory transfers
    //! happening in the priority queue. Larger pool size
    //! helps to speed up operations.
    //! \param w_pool_ pool of blocks that will be used
    //! for writing data for the memory<->disk transfers
    //! happening in the priority queue. Larger pool size
    //! helps to speed up operations.
    STXXL_DEPRECATED(
        priority_queue(prefetch_pool<block_type>& p_pool_, write_pool<block_type>& w_pool_)
        );

    //! Constructs external priority queue object.
    //! \param p_pool_mem memory (in bytes) for prefetch pool that will be used
    //! for data prefetching for the disk<->memory transfers
    //! happening in the priority queue. Larger pool size
    //! helps to speed up operations.
    //! \param w_pool_mem memory (in bytes) for buffered write pool that will be used
    //! for writing data for the memory<->disk transfers
    //! happening in the priority queue. Larger pool size
    //! helps to speed up operations.
    priority_queue(unsigned_type p_pool_mem, unsigned_type w_pool_mem);

    virtual ~priority_queue();

    //! \}

#if 0

    //! swap this priority queue with another one.
    //! Implementation correctness is questionable.
    void swap(priority_queue& obj)
    {
        //swap_1D_arrays(int_mergers,obj.int_mergers,num_int_groups); // does not work in g++ 3.4.3 :( bug?
        for (unsigned_type i = 0; i < num_int_groups; ++i)
            std::swap(int_mergers[i], obj.int_mergers[i]);

        //std::swap(pool,obj.pool);
        //std::swap(pool_owned, obj.pool_owned);
        std::swap(ext_mergers, obj.ext_mergers);
        for (unsigned_type i1 = 0; i1 < total_num_groups; ++i1)
            for (unsigned_type i2 = 0; i2 < (N + 1); ++i2)
                std::swap(group_buffers[i1][i2], obj.group_buffers[i1][i2]);

        STXXL_STATIC_ASSERT(false);
        // Shoot yourself in the foot: group_buffer_current_mins contains pointers into group_buffers ...
        // either recompute them or add/subtract (&this->group_buffers[0][0] - &obj->group_buffers[0][0])
        swap_1D_arrays(group_buffer_current_mins, obj.group_buffer_current_mins, total_num_groups);
        swap_1D_arrays(delete_buffer, obj.delete_buffer, delete_buffer_size + 1);
        std::swap(delete_buffer_current_min, obj.delete_buffer_current_min);
        std::swap(delete_buffer_end, obj.delete_buffer_end);
        std::swap(cmp, obj.cmp);
        std::swap(insert_heap, obj.insert_heap);
        std::swap(num_active_groups, obj.num_active_groups);
        std::swap(size_, obj.size_);
    }
#endif

    //! \name Capacity
    //! \{

    //! Returns number of elements contained.
    //! \return number of elements contained
    size_type size() const;

    //! Returns true if queue has no elements.
    //! \return \b true if queue has no elements, \b false otherwise
    bool empty() const { return (size() == 0); }

    //! \}

    //! \name Operators
    //! \{

    //! Returns "largest" element.
    //!
    //! Returns a const reference to the element at the top of the
    //! priority_queue. The element at the top is guaranteed to be the largest
    //! element in the \b priority queue, as determined by the comparison
    //! function \b comparator_type (the same as the second parameter of
    //! PRIORITY_QUEUE_GENERATOR utility class). That is, for every other
    //! element \b x in the priority_queue, \b comparator_type(Q.top(), x) is
    //! false. Precondition: \c empty() is false.
    const value_type & top() const;

    //! \}

    //! \name Modifiers
    //! \{

    //! Removes the element at the top.
    //!
    //! Removes the element at the top of the priority_queue, that is, the
    //! largest element in the \b priority_queue. Precondition: \c empty() is
    //! \b false. Postcondition: \c size() will be decremented by 1.
    void pop();

    //! Inserts x into the priority_queue.
    //!
    //! Inserts x into the priority_queue. Postcondition: \c size() will be
    //! incremented by 1.
    void push(const value_type& obj);

    //! \}

    //! \name Miscellaneous
    //! \{

    //! Number of bytes consumed by the \b priority_queue from the internal
    //! memory not including pools (see the constructor)
    unsigned_type mem_cons() const
    {
        unsigned_type dynam_alloc_mem = 0;
        //dynam_alloc_mem += w_pool.mem_cons();
        //dynam_alloc_mem += p_pool.mem_cons();
        for (int i = 0; i < num_int_groups; ++i)
            dynam_alloc_mem += int_mergers[i].mem_cons();

        for (int i = 0; i < num_ext_groups; ++i)
            dynam_alloc_mem += ext_mergers[i]->mem_cons();

        return (sizeof(*this) +
                sizeof(ext_merger_type) * num_ext_groups +
                dynam_alloc_mem);
    }

    void dump_sizes() const;
    void dump_params() const;

    //! \}
};

template <class ConfigType>
inline typename priority_queue<ConfigType>::size_type
priority_queue<ConfigType>::size() const
{
    return size_ +
           insert_heap.size() - 1 +
           (delete_buffer_end - delete_buffer_current_min);
}

template <class ConfigType>
inline const typename priority_queue<ConfigType>::value_type &
priority_queue<ConfigType>::top() const
{
    assert(!insert_heap.empty());

    const typename priority_queue<ConfigType>::value_type& t = insert_heap.top();
    if (/*(!insert_heap.empty()) && */ cmp(*delete_buffer_current_min, t))
        return t;
    else
        return *delete_buffer_current_min;
}

template <class ConfigType>
inline void priority_queue<ConfigType>::pop()
{
    //STXXL_VERBOSE1("priority_queue::pop()");
    assert(!insert_heap.empty());

    if (/*(!insert_heap.empty()) && */ cmp(*delete_buffer_current_min, insert_heap.top()))
        insert_heap.pop();
    else
    {
        assert(delete_buffer_current_min < delete_buffer_end);
        ++delete_buffer_current_min;
        if (delete_buffer_current_min == delete_buffer_end)
            refill_delete_buffer();
    }
}

template <class ConfigType>
inline void priority_queue<ConfigType>::push(const value_type& obj)
{
    //STXXL_VERBOSE3("priority_queue::push("<< obj <<")");
    assert(!int_mergers->is_sentinel(obj));
    if (insert_heap.size() == N + 1)
        empty_insert_heap();

    assert(!insert_heap.empty());

    insert_heap.push(obj);
}

////////////////////////////////////////////////////////////////

template <class ConfigType>
priority_queue<ConfigType>::priority_queue(pool_type& pool_)
    : pool(&pool_),
      pool_owned(false),
      delete_buffer_end(delete_buffer + delete_buffer_size),
      insert_heap(N + 2),
      num_active_groups(0), size_(0)
{
    STXXL_VERBOSE_PQ("priority_queue(pool)");
    init();
}

// DEPRECATED
template <class ConfigType>
priority_queue<ConfigType>::priority_queue(prefetch_pool<block_type>& p_pool_, write_pool<block_type>& w_pool_)
    : pool(new pool_type(p_pool_, w_pool_)),
      pool_owned(true),
      delete_buffer_end(delete_buffer + delete_buffer_size),
      insert_heap(N + 2),
      num_active_groups(0), size_(0)
{
    STXXL_VERBOSE_PQ("priority_queue(p_pool, w_pool)");
    init();
}

template <class ConfigType>
priority_queue<ConfigType>::priority_queue(unsigned_type p_pool_mem, unsigned_type w_pool_mem)
    : pool(new pool_type(p_pool_mem / BlockSize, w_pool_mem / BlockSize)),
      pool_owned(true),
      delete_buffer_end(delete_buffer + delete_buffer_size),
      insert_heap(N + 2),
      num_active_groups(0), size_(0)
{
    STXXL_VERBOSE_PQ("priority_queue(pool sizes)");
    init();
}

template <class ConfigType>
void priority_queue<ConfigType>::init()
{
    assert(!cmp(cmp.min_value(), cmp.min_value())); // verify strict weak ordering

    ext_mergers = new ext_merger_type*[num_ext_groups];
    for (unsigned_type j = 0; j < num_ext_groups; ++j) {
        ext_mergers[j] = new ext_merger_type;
        ext_mergers[j]->set_pool(pool);
    }

    value_type sentinel = cmp.min_value();
    insert_heap.push(sentinel);                                // always keep the sentinel
    delete_buffer[delete_buffer_size] = sentinel;              // sentinel
    delete_buffer_current_min = delete_buffer_end;             // empty
    for (unsigned_type i = 0; i < total_num_groups; i++)
    {
        group_buffers[i][N] = sentinel;                        // sentinel
        group_buffer_current_mins[i] = &(group_buffers[i][N]); // empty
    }
}

template <class ConfigType>
priority_queue<ConfigType>::~priority_queue()
{
    STXXL_VERBOSE_PQ("~priority_queue()");
    if (pool_owned)
        delete pool;

    for (unsigned_type j = 0; j < num_ext_groups; ++j)
        delete ext_mergers[j];
    delete[] ext_mergers;
}

//--------------------- Buffer refilling -------------------------------

// refill group_buffers[j] and return number of elements found
template <class ConfigType>
typename priority_queue<ConfigType>::size_type
priority_queue<ConfigType>::refill_group_buffer(unsigned_type group)
{
    STXXL_VERBOSE_PQ("refill_group_buffer(" << group << ")");

    value_type* target;
    size_type length;
    size_type group_size = (group < num_int_groups) ?
                           int_mergers[group].size() :
                           ext_mergers[group - num_int_groups]->size();                        // elements left in segments
    unsigned_type left_elements = group_buffers[group] + N - group_buffer_current_mins[group]; //elements left in target buffer
    if (group_size + left_elements >= size_type(N))
    {                                                                                          // buffer will be filled completely
        target = group_buffers[group];
        length = N - left_elements;
    }
    else
    {
        target = group_buffers[group] + N - group_size - left_elements;
        length = group_size;
    }

    if (length > 0)
    {
        // shift remaininig elements to front
        memmove(target, group_buffer_current_mins[group], left_elements * sizeof(value_type));
        group_buffer_current_mins[group] = target;

        // fill remaining space from group
        if (group < num_int_groups)
            int_mergers[group].multi_merge(target + left_elements,
                                           target + left_elements + length);
        else
            ext_mergers[group - num_int_groups]->multi_merge(
                target + left_elements,
                target + left_elements + length);
    }

    //STXXL_MSG(length + left_elements);
    //std::copy(target,target + length + left_elements,std::ostream_iterator<value_type>(std::cout, "\n"));
#if STXXL_CHECK_ORDER_IN_SORTS
    priority_queue_local::invert_order<typename Config::comparator_type, value_type, value_type> inv_cmp(cmp);
    if (!stxxl::is_sorted(group_buffer_current_mins[group], group_buffers[group] + N, inv_cmp))
    {
        STXXL_VERBOSE_PQ("refill_grp... length: " << length << " left_elements: " << left_elements);
        for (value_type* v = group_buffer_current_mins[group] + 1; v < group_buffer_current_mins[group] + left_elements; ++v)
        {
            if (inv_cmp(*v, *(v - 1)))
            {
                STXXL_MSG("Error in buffer " << group << " at position " << (v - group_buffer_current_mins[group] - 1) << "/" << (v - group_buffer_current_mins[group]) << "   " << *(v - 2) << " " << *(v - 1) << " " << *v << " " << *(v + 1));
            }
        }
        assert(false);
    }
#endif

    return length + left_elements;
}

template <class ConfigType>
void priority_queue<ConfigType>::refill_delete_buffer()
{
    STXXL_VERBOSE_PQ("refill_delete_buffer()");

    size_type total_group_size = 0;
    //num_active_groups is <= 4
    for (unsigned_type i = num_active_groups; i > 0; )
    {
        --i;
        if ((group_buffers[i] + N) - group_buffer_current_mins[i] < delete_buffer_size)
        {
            size_type length = refill_group_buffer(i);
            // max active level dry now?
            if (length == 0 && unsigned(i) == num_active_groups - 1)
                --num_active_groups;

            total_group_size += length;
        }
        else
            total_group_size += delete_buffer_size;  // actually only a sufficient lower bound
    }

    size_type length;
    if (total_group_size >= delete_buffer_size)      // buffer can be filled completely
    {
        length = delete_buffer_size;                 // amount to be copied
        size_ -= size_type(delete_buffer_size);      // amount left in group_buffers
    }
    else
    {
        length = total_group_size;
        assert(size_ == length); // trees and group_buffers get empty
        size_ = 0;
    }

    priority_queue_local::invert_order<typename Config::comparator_type, value_type, value_type> inv_cmp(cmp);

    // now call simplified refill routines
    // which can make the assumption that
    // they find all they are asked in the buffers
    delete_buffer_current_min = delete_buffer_end - length;
    STXXL_VERBOSE_PQ("refill_del... Active groups = " << num_active_groups);
    switch (num_active_groups)
    {
    case 0:
        break;
    case 1:
        std::copy(group_buffer_current_mins[0], group_buffer_current_mins[0] + length, delete_buffer_current_min);
        group_buffer_current_mins[0] += length;
        break;
    case 2:
#if STXXL_PARALLEL && STXXL_PARALLEL_PQ_MULTIWAY_MERGE_DELETE_BUFFER
        {
            std::pair<value_type*, value_type*> seqs[2] =
            {
                std::make_pair(group_buffer_current_mins[0], group_buffers[0] + N),
                std::make_pair(group_buffer_current_mins[1], group_buffers[1] + N)
            };

            parallel::multiway_merge_sentinels(
                seqs, seqs + 2, delete_buffer_current_min, length, inv_cmp);
            // sequence iterators are progressed appropriately

            group_buffer_current_mins[0] = seqs[0].first;
            group_buffer_current_mins[1] = seqs[1].first;
        }
#else
        priority_queue_local::merge2_iterator(
            group_buffer_current_mins[0], group_buffer_current_mins[1],
            delete_buffer_current_min, delete_buffer_current_min + length, cmp);
#endif
        break;
    case 3:
#if STXXL_PARALLEL && STXXL_PARALLEL_PQ_MULTIWAY_MERGE_DELETE_BUFFER
        {
            std::pair<value_type*, value_type*> seqs[3] =
            {
                std::make_pair(group_buffer_current_mins[0], group_buffers[0] + N),
                std::make_pair(group_buffer_current_mins[1], group_buffers[1] + N),
                std::make_pair(group_buffer_current_mins[2], group_buffers[2] + N)
            };

            parallel::multiway_merge_sentinels(
                seqs, seqs + 3, delete_buffer_current_min, length, inv_cmp);
            // sequence iterators are progressed appropriately

            group_buffer_current_mins[0] = seqs[0].first;
            group_buffer_current_mins[1] = seqs[1].first;
            group_buffer_current_mins[2] = seqs[2].first;
        }
#else
        priority_queue_local::merge3_iterator(
            group_buffer_current_mins[0],
            group_buffer_current_mins[1],
            group_buffer_current_mins[2],
            delete_buffer_current_min, delete_buffer_current_min + length, cmp);
#endif
        break;
    case 4:
#if STXXL_PARALLEL && STXXL_PARALLEL_PQ_MULTIWAY_MERGE_DELETE_BUFFER
        {
            std::pair<value_type*, value_type*> seqs[4] =
            {
                std::make_pair(group_buffer_current_mins[0], group_buffers[0] + N),
                std::make_pair(group_buffer_current_mins[1], group_buffers[1] + N),
                std::make_pair(group_buffer_current_mins[2], group_buffers[2] + N),
                std::make_pair(group_buffer_current_mins[3], group_buffers[3] + N)
            };

            parallel::multiway_merge_sentinels(
                seqs, seqs + 4, delete_buffer_current_min, length, inv_cmp);
            // sequence iterators are progressed appropriately

            group_buffer_current_mins[0] = seqs[0].first;
            group_buffer_current_mins[1] = seqs[1].first;
            group_buffer_current_mins[2] = seqs[2].first;
            group_buffer_current_mins[3] = seqs[3].first;
        }
#else
        priority_queue_local::merge4_iterator(
            group_buffer_current_mins[0],
            group_buffer_current_mins[1],
            group_buffer_current_mins[2],
            group_buffer_current_mins[3],
            delete_buffer_current_min, delete_buffer_current_min + length, cmp);
        // side effect free
#endif
        break;
    default:
        STXXL_THROW2(std::runtime_error, "priority_queue<...>::refill_delete_buffer()",
                     "Overflow! The number of buffers on 2nd level in stxxl::priority_queue is currently limited to 4");
    }

#if STXXL_CHECK_ORDER_IN_SORTS
    if (!stxxl::is_sorted(delete_buffer_current_min, delete_buffer_end, inv_cmp))
    {
        for (value_type* v = delete_buffer_current_min + 1; v < delete_buffer_end; ++v)
        {
            if (inv_cmp(*v, *(v - 1)))
            {
                STXXL_MSG("Error at position " << (v - delete_buffer_current_min - 1) << "/" << (v - delete_buffer_current_min) << "   " << *(v - 1) << " " << *v);
            }
        }
        assert(false);
    }
#endif
    //std::copy(delete_buffer_current_min,delete_buffer_current_min + length,std::ostream_iterator<value_type>(std::cout, "\n"));
}

//--------------------------------------------------------------------

// check if space is available on level k and
// empty this level if necessary leading to a recursive call.
// return the level where space was finally available
template <class ConfigType>
unsigned_type priority_queue<ConfigType>::make_space_available(unsigned_type level)
{
    STXXL_VERBOSE_PQ("make_space_available(" << level << ")");
    unsigned_type finalLevel;
    assert(level < total_num_groups);
    assert(level <= num_active_groups);

    if (level == num_active_groups)
        ++num_active_groups;

    const bool spaceIsAvailable_ =
        (level < num_int_groups) ? int_mergers[level].is_space_available()
        : (ext_mergers[level - num_int_groups]->is_space_available());

    if (spaceIsAvailable_)
    {
        finalLevel = level;
    }
    else if (level == total_num_groups - 1)
    {
        size_type capacity = N;
        for (int i = 0; i < num_int_groups; ++i)
            capacity *= IntKMAX;
        for (int i = 0; i < num_ext_groups; ++i)
            capacity *= ExtKMAX;
        STXXL_ERRMSG("priority_queue OVERFLOW - all groups full, size=" << size() <<
                     ", capacity(last externel group (" << num_int_groups + num_ext_groups - 1 << "))=" << capacity);
        dump_sizes();

        unsigned_type extLevel = level - num_int_groups;
        const size_type segmentSize = ext_mergers[extLevel]->size();
        STXXL_VERBOSE1("Inserting segment into last level external: " << level << " " << segmentSize);
        ext_merger_type* overflow_merger = new ext_merger_type;
        overflow_merger->set_pool(pool);
        overflow_merger->append_merger(*ext_mergers[extLevel], segmentSize);
        std::swap(ext_mergers[extLevel], overflow_merger);
        delete overflow_merger;
        finalLevel = level;
    }
    else
    {
        finalLevel = make_space_available(level + 1);

        if (level < num_int_groups - 1)                                           // from internal to internal tree
        {
            unsigned_type segmentSize = int_mergers[level].size();
            value_type* newSegment = new value_type[segmentSize + 1];
            int_mergers[level].multi_merge(newSegment, newSegment + segmentSize); // empty this level

            newSegment[segmentSize] = delete_buffer[delete_buffer_size];          // sentinel
            // for queues where size << #inserts
            // it might make sense to stay in this level if
            // segmentSize < alpha * KNN * k^level for some alpha < 1
            int_mergers[level + 1].append_array(newSegment, segmentSize);
        }
        else
        {
            if (level == num_int_groups - 1) // from internal to external tree
            {
                const unsigned_type segmentSize = int_mergers[num_int_groups - 1].size();
                STXXL_VERBOSE_PQ("make_space... Inserting segment into first level external: " << level << " " << segmentSize);
                ext_mergers[0]->append_merger(int_mergers[num_int_groups - 1], segmentSize);
            }
            else // from external to external tree
            {
                const size_type segmentSize = ext_mergers[level - num_int_groups]->size();
                STXXL_VERBOSE_PQ("make_space... Inserting segment into second level external: " << level << " " << segmentSize);
                ext_mergers[level - num_int_groups + 1]->append_merger(*ext_mergers[level - num_int_groups], segmentSize);
            }
        }
    }
    return finalLevel;
}

// empty the insert heap into the main data structure
template <class ConfigType>
void priority_queue<ConfigType>::empty_insert_heap()
{
    STXXL_VERBOSE_PQ("empty_insert_heap()");
    assert(insert_heap.size() == (N + 1));

    const value_type sup = get_supremum();

    // build new segment
    value_type* newSegment = new value_type[N + 1];
    value_type* newPos = newSegment;

    // put the new data there for now
    //insert_heap.sortTo(newSegment);
    value_type* SortTo = newSegment;

    insert_heap.sort_to(SortTo);

    SortTo = newSegment + N;
    insert_heap.clear();
    insert_heap.push(*SortTo);

    assert(insert_heap.size() == 1);

    newSegment[N] = sup; // sentinel

    // copy the delete_buffer and group_buffers[0] to temporary storage
    // (the temporary can be eliminated using some dirty tricks)
    const unsigned_type tempSize = N + delete_buffer_size;
    value_type temp[tempSize + 1];
    unsigned_type sz1 = current_delete_buffer_size();
    unsigned_type sz2 = current_group_buffer_size(0);
    value_type* pos = temp + tempSize - sz1 - sz2;
    std::copy(delete_buffer_current_min, delete_buffer_current_min + sz1, pos);
    std::copy(group_buffer_current_mins[0], group_buffer_current_mins[0] + sz2, pos + sz1);
    temp[tempSize] = sup; // sentinel

    // refill delete_buffer
    // (using more complicated code it could be made somewhat fuller
    // in certain circumstances)
    priority_queue_local::merge2_iterator(
        pos, newPos,
        delete_buffer_current_min, delete_buffer_current_min + sz1, cmp);

    // refill group_buffers[0]
    // (as above we might want to take the opportunity
    // to make group_buffers[0] fuller)
    priority_queue_local::merge2_iterator(
        pos, newPos,
        group_buffer_current_mins[0], group_buffer_current_mins[0] + sz2, cmp);

    // merge the rest to the new segment
    // note that merge exactly trips into the footsteps
    // of itself
    priority_queue_local::merge2_iterator(pos, newPos,
                                          newSegment, newSegment + N, cmp);

    // and insert it
    unsigned_type freeLevel = make_space_available(0);
    assert(freeLevel == 0 || int_mergers[0].size() == 0);
    int_mergers[0].append_array(newSegment, N);

    // get rid of invalid level 2 buffers
    // by inserting them into tree 0 (which is almost empty in this case)
    if (freeLevel > 0)
    {
        for (int_type i = freeLevel; i >= 0; i--)
        {
            // reverse order not needed
            // but would allow immediate refill

            newSegment = new value_type[current_group_buffer_size(i) + 1]; // with sentinel
            std::copy(group_buffer_current_mins[i], group_buffer_current_mins[i] + current_group_buffer_size(i) + 1, newSegment);
            int_mergers[0].append_array(newSegment, current_group_buffer_size(i));
            group_buffer_current_mins[i] = group_buffers[i] + N;           // empty
        }
    }

    // update size
    size_ += size_type(N);

    // special case if the tree was empty before
    if (delete_buffer_current_min == delete_buffer_end)
        refill_delete_buffer();
}

template <class ConfigType>
void priority_queue<ConfigType>::dump_sizes() const
{
    unsigned_type capacity = N;
    STXXL_MSG("pq::size()\t= " << size());
    STXXL_MSG("  insert_heap\t= " << insert_heap.size() - 1 << "/" << capacity);
    STXXL_MSG("  delete_buffer\t= " << (delete_buffer_end - delete_buffer_current_min) << "/" << delete_buffer_size);
    for (int i = 0; i < num_int_groups; ++i) {
        capacity *= IntKMAX;
        STXXL_MSG("  grp " << i << " int" <<
                  " grpbuf=" << current_group_buffer_size(i) <<
                  " size=" << int_mergers[i].size() << "/" << capacity <<
                  " (" << (int)((double)int_mergers[i].size() * 100.0 / (double)capacity) << "%)" <<
                  " space=" << int_mergers[i].is_space_available());
    }
    for (int i = 0; i < num_ext_groups; ++i) {
        capacity *= ExtKMAX;
        STXXL_MSG("  grp " << i + num_int_groups << " ext" <<
                  " grpbuf=" << current_group_buffer_size(i + num_int_groups) <<
                  " size=" << ext_mergers[i]->size() << "/" << capacity <<
                  " (" << (int)((double)ext_mergers[i]->size() * 100.0 / (double)capacity) << "%)" <<
                  " space=" << ext_mergers[i]->is_space_available());
    }
    dump_params();
}

template <class ConfigType>
void priority_queue<ConfigType>::dump_params() const
{
    STXXL_MSG("params: delete_buffer_size=" << delete_buffer_size << " N=" << N << " IntKMAX=" << IntKMAX << " num_int_groups=" << num_int_groups << " ExtKMAX=" << ExtKMAX << " num_ext_groups=" << num_ext_groups << " BlockSize=" << BlockSize);
}

namespace priority_queue_local {

struct Parameters_for_priority_queue_not_found_Increase_IntMem
{
    enum { fits = false };
    typedef Parameters_for_priority_queue_not_found_Increase_IntMem result;
};

struct dummy
{
    enum { fits = false };
    typedef dummy result;
};

template <internal_size_type ElementSize, internal_size_type IntMem,
          external_size_type MaxItems, internal_size_type BlockSize,
          unsigned_type m_, bool stop = false>
struct find_B_m
{
    typedef find_B_m<ElementSize, IntMem,
                     MaxItems, BlockSize, m_, stop> self_type;

    //! element size
    static const internal_size_type element_size = ElementSize;
    //! internal memory size of PQ
    static const internal_size_type intmem = IntMem;
    //! block size (iterates from 8 MiB downwards)
    static const internal_size_type B = BlockSize;

    //! number of blocks that fit into internal memory (M)
    static const internal_size_type k = IntMem / BlockSize;
    //! number of blocks fitting into buffers of mergers (arity of both
    //! mergers), increased from 1 to 2048 ?-tb
    static const internal_size_type m = m_;
    //! remaining blocks, (freely moving, not necessarily unused) ?-tb
    static const int_type c = k - m_;

    // memory occupied by block must be at least 10 times larger than size of ext sequence

    //! calculated boolean whether the configuration fits into internal memory.
    static const external_size_type fits =
        // need some temporary constant-size internal blocks
        (c > 10) &&
        // satisfy items requirement
        (((k - m) * m * (m * B / (ElementSize * 4 * 1024))) >= MaxItems) &&
        // if we have two ext mergers their degree must be at least 64=m/2
        ((MaxItems < ((k - m) * m / (2 * ElementSize)) * 1024) || m >= 128);

    static const unsigned_type step = 1;

    //! if not fits, recurse into configuration with +step more internal buffers
    typedef typename find_B_m<ElementSize, IntMem, MaxItems, B,
                              m + step, fits || (m + step >= k)>::result candidate1;
    //! if not fits, recurse into configuration with block size halved.
    typedef typename find_B_m<ElementSize, IntMem, MaxItems, B / 2,
                              1, fits || candidate1::fits>::result candidate2;

    //! return a fitting configuration.
    typedef typename IF<fits, self_type, typename IF<candidate1::fits, candidate1, candidate2>::result>::result result;
};

// specialization for the case when no valid parameters are found
template <internal_size_type ElementSize, internal_size_type IntMem,
          external_size_type MaxItems, bool stop>
struct find_B_m<ElementSize, IntMem, MaxItems, 2048, 1, stop>
{
    enum { fits = false };
    typedef Parameters_for_priority_queue_not_found_Increase_IntMem result;
};

// to speedup search
template <internal_size_type ElementSize, internal_size_type IntMem,
          external_size_type MaxItems, unsigned_type BlockSize,
          unsigned_type m_>
struct find_B_m<ElementSize, IntMem, MaxItems, BlockSize, m_, true>
{
    enum { fits = false };
    typedef dummy result;
};

// start search
template <internal_size_type ElementSize, internal_size_type IntMem,
          external_size_type MaxItems>
struct find_settings
{
    // start from block size (8*1024*1024) bytes
    typedef typename find_B_m<ElementSize, IntMem,
                              MaxItems, (8* 1024* 1024), 1>::result result;
};

struct Parameters_not_found_Try_to_change_the_Tune_parameter
{
    typedef Parameters_not_found_Try_to_change_the_Tune_parameter result;
};

template <unsigned_type AI_, unsigned_type X_, unsigned_type CriticalSize>
struct compute_N
{
    typedef compute_N<AI_, X_, CriticalSize> Self;

    static const unsigned_type X = X_;
    static const unsigned_type AI = AI_;
    static const unsigned_type N = X / (AI * AI);     // two stage internal

    typedef typename IF<(N >= CriticalSize), Self, typename compute_N<AI / 2, X, CriticalSize>::result>::result result;
};

template <unsigned_type X_, unsigned_type CriticalSize_>
struct compute_N<1, X_, CriticalSize_>
{
    typedef Parameters_not_found_Try_to_change_the_Tune_parameter result;
};

} // namespace priority_queue_local

//! \}

//! \addtogroup stlcont
//! \{

//! Priority queue type generator. \n
//! <b> Introduction </b> to priority queue container: see \ref tutorial_pqueue tutorial. \n
//! <b> Design and Internals </b> of priority queue container: see \ref design_pqueue.
//!
//! \tparam ValueType type of the contained objects (POD with no references to internal memory)
//!
//! \tparam CompareType the comparator type used to determine whether one element is
//! smaller than another element.
//!
//! \tparam IntMemory upper limit for internal memory consumption in bytes.
//!
//! \tparam MaxItems upper limit for number of elements contained in the priority queue (in 1024 units). <BR>
//! Example: if you are sure that priority queue contains no more than
//! one million elements in a time, then the right parameter is (1000000 / 1024) = 976.
//!
//! \tparam Tune tuning parameter for meta-program search. <BR>
//! Try to play with it if the code does not compile (larger than default
//! values might help). Code does not compile if no suitable internal
//! parameters were found for given IntMemory and MaxItems. It might also
//! happen that given IntMemory is too small for given MaxItems, try larger
//! values.
template <class ValueType,
          class CompareType,
          internal_size_type IntMemory,
          external_size_type MaxItems,
          unsigned Tune = 6>
class PRIORITY_QUEUE_GENERATOR
{
public:
    // actual calculation of B, m, k and element_size
    typedef typename priority_queue_local::find_settings<sizeof(ValueType), IntMemory, MaxItems>::result settings;
    enum {
        B = settings::B,
        m = settings::m,
        X = B * (settings::k - m) / settings::element_size,  // interpretation of result
        Buffer1Size = 32                                     // fixed
    };
    // derivation of N, AI, AE
    typedef typename priority_queue_local::compute_N<(1 << Tune), X, 4* Buffer1Size>::result ComputeN;
    enum
    {
        N = ComputeN::N,
        AI = ComputeN::AI,
        AE = (m / 2 < 2) ? 2 : (m / 2)            // at least 2
    };

    // Estimation of maximum internal memory consumption (in bytes)
    static const unsigned_type EConsumption = X * settings::element_size + settings::B * AE + ((MaxItems / X) / AE) * settings::B * 1024;

    /*
        unsigned BufferSize1_ = 32, // equalize procedure call overheads etc.
        unsigned N_ = 512,          // bandwidth
        unsigned IntKMAX_ = 64,     // maximal arity for internal mergers
        unsigned IntLevels_ = 4,
        unsigned BlockSize = (2*1024*1024),
        unsigned ExtKMAX_ = 64,     // maximal arity for external mergers
        unsigned ExtLevels_ = 2,
     */
    typedef priority_queue<priority_queue_config<ValueType, CompareType, Buffer1Size, N, AI, 2, B, AE, 2> > result;
};

//! \}

STXXL_END_NAMESPACE

namespace std {

template <class ConfigType>
void swap(stxxl::priority_queue<ConfigType>& a,
          stxxl::priority_queue<ConfigType>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_CONTAINERS_PRIORITY_QUEUE_HEADER
// vim: et:ts=4:sw=4
