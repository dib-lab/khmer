/***************************************************************************
 *  include/stxxl/bits/algo/sort_helper.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_SORT_HELPER_HEADER
#define STXXL_ALGO_SORT_HELPER_HEADER

#include <algorithm>
#include <functional>
#include <stxxl/bits/algo/run_cursor.h>
#include <stxxl/bits/verbose.h>

STXXL_BEGIN_NAMESPACE

//! \internal
namespace sort_helper {

template <typename StrictWeakOrdering>
inline void verify_sentinel_strict_weak_ordering(StrictWeakOrdering cmp)
{
    STXXL_ASSERT(!cmp(cmp.min_value(), cmp.min_value()));
    STXXL_ASSERT(cmp(cmp.min_value(), cmp.max_value()));
    STXXL_ASSERT(!cmp(cmp.max_value(), cmp.min_value()));
    STXXL_ASSERT(!cmp(cmp.max_value(), cmp.max_value()));
}

template <typename BlockType, typename ValueType = typename BlockType::value_type>
struct trigger_entry
{
    typedef BlockType block_type;
    typedef typename block_type::bid_type bid_type;
    typedef ValueType value_type;

    bid_type bid;
    value_type value;

    operator bid_type ()
    {
        return bid;
    }
};

template <typename TriggerEntryType, typename ValueCmp>
struct trigger_entry_cmp
    : public std::binary_function<TriggerEntryType, TriggerEntryType, bool>
{
    typedef TriggerEntryType trigger_entry_type;
    ValueCmp cmp;
    trigger_entry_cmp(ValueCmp c) : cmp(c) { }
    trigger_entry_cmp(const trigger_entry_cmp& a) : cmp(a.cmp) { }
    bool operator () (const trigger_entry_type& a, const trigger_entry_type& b) const
    {
        return cmp(a.value, b.value);
    }
};

template <typename BlockType,
          typename PrefetcherType,
          typename ValueCmp>
struct run_cursor2_cmp
    : public std::binary_function<
          run_cursor2<BlockType, PrefetcherType>,
          run_cursor2<BlockType, PrefetcherType>,
          bool
          >
{
    typedef BlockType block_type;
    typedef PrefetcherType prefetcher_type;
    typedef ValueCmp value_cmp;

    typedef run_cursor2<block_type, prefetcher_type> cursor_type;
    value_cmp cmp;

    run_cursor2_cmp(value_cmp c) : cmp(c) { }
    run_cursor2_cmp(const run_cursor2_cmp& a) : cmp(a.cmp) { }
    inline bool operator () (const cursor_type& a, const cursor_type& b) const
    {
        if (UNLIKELY(b.empty()))
            return true;
        // sentinel emulation
        if (UNLIKELY(a.empty()))
            return false;
        // sentinel emulation

        return (cmp(a.current(), b.current()));
    }
};

// this function is used by parallel mergers
template <typename SequenceVector, typename ValueType, typename Comparator>
inline unsigned_type
count_elements_less_equal(const SequenceVector& seqs,
                          const ValueType& bound, Comparator cmp)
{
    typedef typename SequenceVector::size_type seqs_size_type;
    typedef typename SequenceVector::value_type::first_type iterator;
    unsigned_type count = 0;

    for (seqs_size_type i = 0; i < seqs.size(); ++i)
    {
        iterator position = std::upper_bound(seqs[i].first, seqs[i].second, bound, cmp);
        STXXL_VERBOSE1("less equal than " << position - seqs[i].first);
        count += position - seqs[i].first;
    }
    STXXL_VERBOSE1("finished loop");
    return count;
}

// this function is used by parallel mergers
template <typename SequenceVector, typename BufferPtrVector, typename Prefetcher>
inline void
refill_or_remove_empty_sequences(SequenceVector& seqs,
                                 BufferPtrVector& buffers,
                                 Prefetcher& prefetcher)
{
    typedef typename SequenceVector::size_type seqs_size_type;

    for (seqs_size_type i = 0; i < seqs.size(); ++i)
    {
        if (seqs[i].first == seqs[i].second)                    // run empty
        {
            if (prefetcher.block_consumed(buffers[i]))
            {
                seqs[i].first = buffers[i]->begin();            // reset iterator
                seqs[i].second = buffers[i]->end();
                STXXL_VERBOSE1("block ran empty " << i);
            }
            else
            {
                seqs.erase(seqs.begin() + i);                   // remove this sequence
                buffers.erase(buffers.begin() + i);
                STXXL_VERBOSE1("seq removed " << i);
                --i;                                            // don't skip the next sequence
            }
        }
    }
}

} // namespace sort_helper

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_SORT_HELPER_HEADER
// vim: et:ts=4:sw=4
