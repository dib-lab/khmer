/***************************************************************************
 *  include/stxxl/bits/stream/sorted_runs.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2005 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_STREAM_SORTED_RUNS_HEADER
#define STXXL_STREAM_SORTED_RUNS_HEADER

#include <vector>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/algo/adaptor.h>
#include <stxxl/bits/common/counting_ptr.h>

STXXL_BEGIN_NAMESPACE

namespace stream {

//! \addtogroup streampack Stream Package
//! \{

////////////////////////////////////////////////////////////////////////
//     SORTED RUNS                                                    //
////////////////////////////////////////////////////////////////////////

//! All sorted runs of a sort operation.
template <typename TriggerEntryType, typename CompareType>
struct sorted_runs : private noncopyable, public counted_object
{
    typedef TriggerEntryType trigger_entry_type;
    typedef typename trigger_entry_type::block_type block_type;
    //! may differ from trigger_entry_type::value_type
    typedef typename block_type::value_type value_type;
    typedef std::vector<trigger_entry_type> run_type;
    typedef std::vector<value_type> small_run_type;
    typedef stxxl::external_size_type size_type;
    typedef typename std::vector<run_type>::size_type run_index_type;

    typedef CompareType cmp_type;

    //! total number of elements in all runs
    size_type elements;

    //! vector of runs (containing vectors of block ids)
    std::vector<run_type> runs;

    //! vector of the number of elements in each individual run
    std::vector<size_type> runs_sizes;

    //! Small sort optimization:
    // if the input is small such that its total size is at most B
    // (block_type::size) then input is sorted internally and kept in the
    // array "small_run"
    small_run_type small_run;

public:
    sorted_runs()
        : elements(0)
    { }

    ~sorted_runs()
    {
        deallocate_blocks();
    }

    //! Clear the internal state of the object: release all runs and reset.
    void clear()
    {
        deallocate_blocks();

        elements = 0;
        runs.clear();
        runs_sizes.clear();
        small_run.clear();
    }

    //! Add a new run with given number of elements
    void add_run(const run_type& run, size_type run_size)
    {
        runs.push_back(run);
        runs_sizes.push_back(run_size);
        elements += run_size;
    }

    //! Swap contents with another object. This is used by the recursive
    //! merger to swap in a sorted_runs object with fewer runs.
    void swap(sorted_runs& b)
    {
        std::swap(elements, b.elements);
        std::swap(runs, b.runs);
        std::swap(runs_sizes, b.runs_sizes);
        std::swap(small_run, b.small_run);
    }

private:
    //! Deallocates the blocks which the runs occupy.
    //!
    //! \remark There is no need in calling this method, the blocks are
    //! deallocated by the destructor. However, if you wish to reuse the
    //! object, then this function can be used to clear its state.
    void deallocate_blocks()
    {
        block_manager* bm = block_manager::get_instance();
        for (unsigned_type i = 0; i < runs.size(); ++i)
        {
            bm->delete_blocks(make_bid_iterator(runs[i].begin()),
                              make_bid_iterator(runs[i].end()));
        }
    }
};

//! \}

} // namespace stream

STXXL_END_NAMESPACE

#endif // !STXXL_STREAM_SORTED_RUNS_HEADER
// vim: et:ts=4:sw=4
