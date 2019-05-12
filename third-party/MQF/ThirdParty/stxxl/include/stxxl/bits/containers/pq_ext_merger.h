/***************************************************************************
 *  include/stxxl/bits/containers/pq_ext_merger.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 1999 Peter Sanders <sanders@mpi-sb.mpg.de>
 *  Copyright (C) 2003, 2004, 2007 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007-2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2007-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_PQ_EXT_MERGER_HEADER
#define STXXL_CONTAINERS_PQ_EXT_MERGER_HEADER

#include <stxxl/bits/containers/pq_mergers.h>
#include <deque>

STXXL_BEGIN_NAMESPACE

//! \addtogroup stlcontinternals
//!
//! \{

/*! \internal
 */
namespace priority_queue_local {

/*!
 * External merger, based on the loser tree data structure.
 * \param Arity  maximum arity of merger, does not need to be a power of 2
 */
template <class BlockType, class CompareType, unsigned Arity,
          class AllocStr = STXXL_DEFAULT_ALLOC_STRATEGY>
class ext_merger : private noncopyable
{
public:
    //! class is parameterized by the block of the external arrays
    typedef BlockType block_type;
    typedef CompareType compare_type;

    // max_arity / 2  <  arity  <=  max_arity
    enum { arity = Arity, max_arity = 1UL << (LOG2<Arity>::ceil) };

    typedef AllocStr alloc_strategy;

    typedef typename block_type::bid_type bid_type;
    typedef typename block_type::value_type value_type;

    typedef typename std::deque<bid_type> bid_container_type;

    typedef read_write_pool<block_type> pool_type;

    //! our type
    typedef ext_merger<BlockType, CompareType, Arity, AllocStr> self_type;

#if STXXL_PARALLEL && STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL
    //! type of embedded adapter to parallel multiway_merge
    typedef parallel_merger_adapter<self_type, CompareType, Arity> tree_type;
#else
    //! type of embedded loser tree
    typedef loser_tree<self_type, CompareType, Arity> tree_type;
#endif

    //! size type of total number of item in merger
    typedef external_size_type size_type;

public:
    struct sequence_state : private noncopyable
    {
        block_type* block;          //!< current block
        unsigned_type current;      //!< current index in current block
        bid_container_type bids;    //!< list of blocks forming this sequence
        compare_type cmp;
        ext_merger* merger;
        bool allocated;

        //! \returns current element
        const value_type& operator * () const
        {
            return (*block)[current];
        }

        sequence_state()
            : block(NULL), current(0),
              merger(NULL),
              allocated(false)
        { }

        ~sequence_state()
        {
            STXXL_VERBOSE2("ext_merger sequence_state::~sequence_state()");

            block_manager* bm = block_manager::get_instance();
            bm->delete_blocks(bids.begin(), bids.end());
        }

        void make_inf()
        {
            current = 0;
            (*block)[0] = cmp.min_value();
        }

        bool is_sentinel(const value_type& a) const
        {
            return !(cmp(cmp.min_value(), a));
        }

        bool not_sentinel(const value_type& a) const
        {
            return cmp(cmp.min_value(), a);
        }

        void swap(sequence_state& obj)
        {
            if (&obj != this)
            {
                std::swap(current, obj.current);
                std::swap(block, obj.block);
                std::swap(bids, obj.bids);
                assert(merger == obj.merger);
                std::swap(allocated, obj.allocated);
            }
        }

        sequence_state& operator ++ ()
        {
            assert(not_sentinel((*block)[current]));
            assert(current < block->size);

            ++current;

            if (current == block->size)
            {
                STXXL_VERBOSE2("ext_merger sequence_state operator++ crossing block border ");
                // go to the next block
                if (bids.empty()) // if there is no next block
                {
                    STXXL_VERBOSE2("ext_merger sequence_state operator++ it was the last block in the sequence ");
                    // swap memory area and delete other object.
                    bid_container_type to_delete;
                    std::swap(to_delete, bids);
                    make_inf();
                }
                else
                {
                    STXXL_VERBOSE2("ext_merger sequence_state operator++ there is another block ");
                    bid_type bid = bids.front();
                    bids.pop_front();
                    merger->pool->hint(bid);
                    if (!(bids.empty()))
                    {
                        STXXL_VERBOSE2("ext_merger sequence_state operator++ more blocks exist in a sequence, hinting the next");
                        merger->pool->hint(bids.front());
                    }
                    merger->pool->read(block, bid)->wait();
                    STXXL_VERBOSE2("first element of read block " << bid << " " << *(block->begin()) << " cached in " << block);
                    if (!(bids.empty()))
                        merger->pool->hint(bids.front());  // re-hint, reading might have made a block free
                    block_manager::get_instance()->delete_block(bid);
                    current = 0;
                }
            }
            return *this;
        }
    };

    //! type of sequences in which the values are stored: external arrays
    typedef sequence_state sequence_type;

protected:
    //! loser tree instance
    tree_type tree;

    //! sequence including current position, dereference gives current element
    sequence_state states[max_arity];

    //! read and writer block pool
    pool_type* pool;

    //! a memory block filled with sentinel values
    block_type* sentinel_block;

    //! total number of elements stored
    size_type m_size;

public:
    ext_merger(const compare_type& c = compare_type()) // TODO: pass pool as parameter
        : tree(c, *this),
          pool(NULL),
          m_size(0)
    {
        init();

        tree.initialize();
    }

    virtual ~ext_merger()
    {
        STXXL_VERBOSE1("ext_merger::~ext_merger()");
        for (unsigned_type i = 0; i < arity; ++i)
        {
            delete states[i].block;
        }
        delete sentinel_block;
    }

    void set_pool(pool_type* pool_)
    {
        pool = pool_;
    }

public:
    //! \name Interface for loser_tree
    //! \{

    //! is this segment empty ?
    bool is_array_empty(unsigned_type slot) const
    {
        return is_sentinel(*(states[slot]));
    }

    //! Is this segment allocated? Otherwise it's empty, already on the stack
    //! of free segment indices and can be reused.
    bool is_array_allocated(unsigned_type slot) const
    {
        return states[slot].allocated;
    }

    //! Return the item sequence of the given slot
    sequence_type & get_array(unsigned_type slot)
    {
        return states[slot];
    }

    //! Swap contents of arrays a and b
    void swap_arrays(unsigned_type a, unsigned_type b)
    {
        states[a].swap(states[b]);
    }

    //! Set a usually empty array to the sentinel
    void make_array_sentinel(unsigned_type a)
    {
        states[a].make_inf();
    }

    //! free an empty segment.
    void free_array(unsigned_type slot)
    {
        STXXL_VERBOSE1("ext_merger::free_array() deleting array " << slot << " allocated=" << int(is_array_allocated(slot)));
        assert(is_array_allocated(slot));
        states[slot].allocated = false;
        states[slot].make_inf();

        // free player in loser tree
        tree.free_player(slot);
    }

    //! Hint (prefetch) first non-internal (actually second) block of each
    //! sequence.
    void prefetch_arrays()
    {
        for (unsigned_type i = 0; i < tree.k; ++i)
        {
            if (!states[i].bids.empty())
                pool->hint(states[i].bids.front());
        }
    }

    //! \}

protected:
    void init()
    {
        STXXL_VERBOSE2("ext_merger::init()");

        sentinel_block = NULL;
        if (arity < max_arity)
        {
            sentinel_block = new block_type;
            for (unsigned_type i = 0; i < block_type::size; ++i)
                (*sentinel_block)[i] = tree.cmp.min_value();
            if (arity + 1 == max_arity) {
                // same memory consumption, but smaller merge width, better use arity = max_arity
                STXXL_ERRMSG("inefficient PQ parameters for ext_merger: arity + 1 == max_arity");
            }
        }

        for (unsigned_type i = 0; i < max_arity; ++i)
        {
            states[i].merger = this;
            if (i < arity)
                states[i].block = new block_type;
            else
                states[i].block = sentinel_block;

            states[i].make_inf();
        }
    }

#if 0
    void swap(ext_merger& obj)
    {
        std::swap(cmp, obj.cmp);
        std::swap(k, obj.k);
        std::swap(logK, obj.logK);
        std::swap(m_size, obj.m_size);
        swap_1D_arrays(states, obj.states, max_arity);

        // std::swap(pool,obj.pool);
    }
#endif

public:
    unsigned_type mem_cons() const // only rough estimation
    {
        return (STXXL_MIN<unsigned_type>(arity + 1, max_arity) * block_type::raw_size);
    }

    //! Whether there is still space for new array
    bool is_space_available() const
    {
        return tree.is_space_available();
    }

    //! True if a is the sentinel value
    bool is_sentinel(const value_type& a) const
    {
        return tree.is_sentinel(a);
    }

    //! Return the number of items in the arrays
    size_type size() const
    {
        return m_size;
    }

    /*!
      \param bidlist list of blocks to insert
      \param first_block the first block of the sequence, before bidlist
      \param first_size number of elements in the first block
      \param slot slot to insert into
    */
    void insert_segment(bid_container_type& bidlist, block_type* first_block,
                        unsigned_type first_size, unsigned_type slot)
    {
        STXXL_VERBOSE1("ext_merger::insert_segment(bidlist,...) " << this << " " << bidlist.size() << " " << slot);
        assert(!is_array_allocated(slot));
        assert(first_size > 0);

        sequence_state& new_sequence = states[slot];
        new_sequence.current = block_type::size - first_size;
        std::swap(new_sequence.block, first_block);
        delete first_block;
        std::swap(new_sequence.bids, bidlist);
        new_sequence.allocated = true;
        assert(is_array_allocated(slot));
    }

    //! Merge all items from another merger and insert the resulting external
    //! array into the merger. Requires: is_space_available() == 1
    template <class Merger>
    void append_merger(Merger& another_merger, size_type segment_size)
    {
        STXXL_VERBOSE1("ext_merger::append_merger(merger,...)" << this);

        if (segment_size == 0)
        {
            // deallocate memory ?
            STXXL_VERBOSE1("Merged segment with zero size.");
            return;
        }

        // allocate a new player slot
        unsigned_type index = tree.new_player();

        // construct new sorted array from merger
        assert(segment_size);
        unsigned_type nblocks = (unsigned_type)(segment_size / block_type::size);
        //assert(nblocks); // at least one block
        STXXL_VERBOSE1("ext_merger::insert_segment nblocks=" << nblocks);
        if (nblocks == 0)
        {
            STXXL_VERBOSE1("ext_merger::insert_segment(merger,...) WARNING: inserting a segment with " <<
                           nblocks << " blocks");
            STXXL_VERBOSE1("THIS IS INEFFICIENT: TRY TO CHANGE PRIORITY QUEUE PARAMETERS");
        }
        unsigned_type first_size = (unsigned_type)(segment_size % block_type::size);
        if (first_size == 0)
        {
            first_size = block_type::size;
            --nblocks;
        }

        // allocate blocks
        block_manager* bm = block_manager::get_instance();
        bid_container_type bids(nblocks);
        bm->new_blocks(alloc_strategy(), bids.begin(), bids.end());
        block_type* first_block = new block_type;

        another_merger.multi_merge(
            first_block->begin() + (block_type::size - first_size),
            first_block->end());

        STXXL_VERBOSE1("last element of first block " << *(first_block->end() - 1));
        assert(!tree.cmp(*(first_block->begin() + (block_type::size - first_size)), *(first_block->end() - 1)));

        assert(pool->size_write() > 0);

        for (typename bid_container_type::iterator curbid = bids.begin(); curbid != bids.end(); ++curbid)
        {
            block_type* b = pool->steal();
            another_merger.multi_merge(b->begin(), b->end());
            STXXL_VERBOSE1("first element of following block " << *curbid << " " << *(b->begin()));
            STXXL_VERBOSE1("last element of following block " << *curbid << " " << *(b->end() - 1));
            assert(!tree.cmp(*(b->begin()), *(b->end() - 1)));
            pool->write(b, *curbid);
            STXXL_VERBOSE1("written to block " << *curbid << " cached in " << b);
        }

        insert_segment(bids, first_block, first_size, index);

        m_size += segment_size;

        // propagate new information up the tree
        tree.update_on_insert((index + tree.k) >> 1, *(states[index]), index);
    }

    // delete the (length = end-begin) smallest elements and write them to [begin..end)
    // empty segments are deallocated
    // requires:
    // - there are at least length elements
    // - segments are ended by sentinels
    template <class OutputIterator>
    void multi_merge(OutputIterator begin, OutputIterator end)
    {
        assert((size_type)(end - begin) <= m_size);

#if STXXL_PARALLEL && STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL
        multi_merge_parallel(begin, end);
#else       // STXXL_PARALLEL && STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL
        tree.multi_merge(begin, end);
        m_size -= end - begin;
#endif
    }

#if STXXL_PARALLEL && STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL

protected:
    //! extract the (length = end - begin) smallest elements using parallel
    //! multiway_merge.

    template <class OutputIterator>
    void multi_merge_parallel(OutputIterator begin, OutputIterator end)
    {
        const unsigned_type& k = tree.k;

        if (begin == end)
            return;

        typedef stxxl::int64 diff_type;

        typedef std::pair<typename block_type::iterator, typename block_type::iterator> sequence;

        std::vector<sequence> seqs;
        std::vector<unsigned_type> orig_seq_index;

        invert_order<compare_type, value_type, value_type> inv_cmp(tree.cmp);

        for (unsigned_type i = 0; i < k; ++i) //initialize sequences
        {
            if (states[i].current == states[i].block->size || is_sentinel(*states[i]))
                continue;

            seqs.push_back(std::make_pair(states[i].block->begin() + states[i].current, states[i].block->end()));
            orig_seq_index.push_back(i);

#if STXXL_CHECK_ORDER_IN_SORTS
            if (!is_sentinel(*seqs.back().first) && !stxxl::is_sorted(seqs.back().first, seqs.back().second, inv_cmp))
            {
                STXXL_VERBOSE1("length " << i << " " << (seqs.back().second - seqs.back().first));
                for (value_type* v = seqs.back().first + 1; v < seqs.back().second; ++v)
                {
                    if (inv_cmp(*v, *(v - 1)))
                    {
                        STXXL_VERBOSE1("Error at position " << i << "/" << (v - seqs.back().first - 1) << "/" << (v - seqs.back().first) << "   " << *(v - 1) << " " << *v);
                    }
                    if (is_sentinel(*v))
                    {
                        STXXL_VERBOSE1("Wrong sentinel at position " << (v - seqs.back().first));
                    }
                }
                assert(false);
            }
#endif

            // Hint first non-internal (actually second) block of this sequence.
            if (!states[i].bids.empty())
                pool->hint(states[i].bids.front());
        }

        assert(seqs.size() > 0);

#if STXXL_CHECK_ORDER_IN_SORTS
        value_type last_elem;
#endif

        // elements still to merge for this output block
        diff_type rest = end - begin;

        while (rest > 0)
        {
            // minimum of the sequences' last elements
            value_type min_last = tree.cmp.min_value();

            diff_type total_size = 0;

            for (unsigned_type i = 0; i < seqs.size(); ++i)
            {
                diff_type seq_i_size = seqs[i].second - seqs[i].first;
                if (seq_i_size > 0)
                {
                    total_size += seq_i_size;
                    if (inv_cmp(*(seqs[i].second - 1), min_last))
                        min_last = *(seqs[i].second - 1);

                    STXXL_VERBOSE1("front block of seq " << i << ": front=" << *(seqs[i].first) << " back=" << *(seqs[i].second - 1) << " len=" << seq_i_size);
                }
                else {
                    STXXL_VERBOSE1("front block of seq " << i << ": empty");
                }
            }

            assert(total_size > 0);
            assert(!is_sentinel(min_last));

            STXXL_VERBOSE1("min_last " << min_last << " total size " << total_size << " num_seq " << seqs.size());

            diff_type less_equal_than_min_last = 0;
            //locate this element in all sequences
            for (unsigned_type i = 0; i < seqs.size(); ++i)
            {
                //assert(seqs[i].first < seqs[i].second);

                typename block_type::iterator position =
                    std::upper_bound(seqs[i].first, seqs[i].second, min_last, inv_cmp);

                //no element larger than min_last is merged

                STXXL_VERBOSE1("seq " << i << ": " << (position - seqs[i].first) << " greater equal than " << min_last);

                less_equal_than_min_last += (position - seqs[i].first);
            }

            // at most rest elements
            diff_type output_size = STXXL_MIN(less_equal_than_min_last, rest);

            STXXL_VERBOSE1("output_size=" << output_size << " = min(leq_t_ml=" << less_equal_than_min_last << ", rest=" << rest << ")");

            assert(output_size > 0);

            //main call

            // sequence iterators are progressed appropriately:
            begin = parallel::multiway_merge(
                seqs.begin(), seqs.end(), begin, output_size, inv_cmp);

            rest -= output_size;
            m_size -= output_size;

            for (unsigned_type i = 0; i < seqs.size(); ++i)
            {
                sequence_state& state = states[orig_seq_index[i]];

                state.current = seqs[i].first - state.block->begin();

                assert(seqs[i].first <= seqs[i].second);

                if (seqs[i].first == seqs[i].second)
                {
                    // has run empty?

                    assert(state.current == state.block->size);
                    if (state.bids.empty())
                    {
                        // if there is no next block
                        STXXL_VERBOSE1("seq " << i << ": ext_merger::multi_merge(...) it was the last block in the sequence ");
                        state.make_inf();
                    }
                    else
                    {
#if STXXL_CHECK_ORDER_IN_SORTS
                        last_elem = *(seqs[i].second - 1);
#endif
                        STXXL_VERBOSE1("seq " << i << ": ext_merger::multi_merge(...) there is another block ");
                        bid_type bid = state.bids.front();
                        state.bids.pop_front();
                        pool->hint(bid);
                        if (!(state.bids.empty()))
                        {
                            STXXL_VERBOSE2("seq " << i << ": ext_merger::multi_merge(...) more blocks exist, hinting the next");
                            pool->hint(state.bids.front());
                        }
                        pool->read(state.block, bid)->wait();
                        STXXL_VERBOSE1("seq " << i << ": first element of read block " << bid << " " << *(state.block->begin()) << " cached in " << state.block);
                        if (!(state.bids.empty()))
                            pool->hint(state.bids.front());  // re-hint, reading might have made a block free
                        state.current = 0;
                        seqs[i] = std::make_pair(state.block->begin() + state.current, state.block->end());
                        block_manager::get_instance()->delete_block(bid);

#if STXXL_CHECK_ORDER_IN_SORTS
                        STXXL_VERBOSE1("before " << last_elem << " after " << *seqs[i].first << " newly loaded block " << bid);
                        if (!stxxl::is_sorted(seqs[i].first, seqs[i].second, inv_cmp))
                        {
                            STXXL_VERBOSE1("length " << i << " " << (seqs[i].second - seqs[i].first));
                            for (value_type* v = seqs[i].first + 1; v < seqs[i].second; ++v)
                            {
                                if (inv_cmp(*v, *(v - 1)))
                                {
                                    STXXL_VERBOSE1("Error at position " << i << "/" << (v - seqs[i].first - 1) << "/" << (v - seqs[i].first) << "   " << *(v - 1) << " " << *v);
                                }
                                if (is_sentinel(*v))
                                {
                                    STXXL_VERBOSE1("Wrong sentinel at position " << (v - seqs[i].first));
                                }
                            }
                            assert(false);
                        }
#endif
                    }
                }
            }
        }   //while (rest > 1)

        for (unsigned_type i = 0; i < seqs.size(); ++i)
        {
            unsigned_type seg = orig_seq_index[i];
            if (is_array_empty(seg))
            {
                STXXL_VERBOSE1("deallocated " << seg);
                free_array(seg);
            }
        }

        tree.maybe_compact();
    }
#endif  // STXXL_PARALLEL && STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL
}; // class ext_merger

} // namespace priority_queue_local

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_PQ_EXT_MERGER_HEADER
// vim: et:ts=4:sw=4
