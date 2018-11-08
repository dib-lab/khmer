/***************************************************************************
 *  include/stxxl/bits/algo/sort.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2006 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2008-2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_SORT_HEADER
#define STXXL_ALGO_SORT_HEADER

#include <functional>

#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/common/rand.h>
#include <stxxl/bits/mng/adaptor.h>
#include <stxxl/bits/common/simple_vector.h>
#include <stxxl/bits/common/settings.h>
#include <stxxl/bits/common/error_handling.h>
#include <stxxl/bits/mng/block_alloc_interleaved.h>
#include <stxxl/bits/io/request_operations.h>
#include <stxxl/bits/algo/sort_base.h>
#include <stxxl/bits/algo/sort_helper.h>
#include <stxxl/bits/algo/intksort.h>
#include <stxxl/bits/algo/adaptor.h>
#include <stxxl/bits/algo/async_schedule.h>
#include <stxxl/bits/mng/block_prefetcher.h>
#include <stxxl/bits/mng/buf_writer.h>
#include <stxxl/bits/algo/run_cursor.h>
#include <stxxl/bits/algo/losertree.h>
#include <stxxl/bits/algo/inmemsort.h>
#include <stxxl/bits/parallel.h>
#include <stxxl/bits/common/is_sorted.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup stlalgo
//! \{

/*! \internal
 */
namespace sort_local {

template <typename BlockType, typename BidType>
struct read_next_after_write_completed
{
    typedef BlockType block_type;
    block_type* block;
    BidType bid;
    request_ptr* req;
    void operator () (request* /*completed_req*/)
    {
        * req = block->read(bid);
    }
};

template <
    typename BlockType,
    typename RunType,
    typename InputBidIterator,
    typename ValueCmp>
void
create_runs(
    InputBidIterator it,
    RunType** runs,
    int_type nruns,
    int_type _m,
    ValueCmp cmp)
{
    typedef BlockType block_type;
    typedef RunType run_type;

    typedef typename block_type::bid_type bid_type;
    STXXL_VERBOSE1("stxxl::create_runs nruns=" << nruns << " m=" << _m);

    int_type m2 = _m / 2;
    block_manager* bm = block_manager::get_instance();
    block_type* Blocks1 = new block_type[m2];
    block_type* Blocks2 = new block_type[m2];
    bid_type* bids1 = new bid_type[m2];
    bid_type* bids2 = new bid_type[m2];
    request_ptr* read_reqs1 = new request_ptr[m2];
    request_ptr* read_reqs2 = new request_ptr[m2];
    request_ptr* write_reqs = new request_ptr[m2];
    read_next_after_write_completed<block_type, bid_type>* next_run_reads =
        new read_next_after_write_completed<block_type, bid_type>[m2];

    disk_queues::get_instance()->set_priority_op(request_queue::WRITE);

    int_type i;
    int_type run_size = 0;

    assert(nruns >= 2);

    run_size = runs[0]->size();
    assert(run_size == m2);

    for (i = 0; i < run_size; ++i)
    {
        STXXL_VERBOSE1("stxxl::create_runs posting read " << Blocks1[i].elem);
        bids1[i] = *(it++);
        read_reqs1[i] = Blocks1[i].read(bids1[i]);
    }

    run_size = runs[1]->size();

    for (i = 0; i < run_size; ++i)
    {
        STXXL_VERBOSE1("stxxl::create_runs posting read " << Blocks2[i].elem);
        bids2[i] = *(it++);
        read_reqs2[i] = Blocks2[i].read(bids2[i]);
    }

    for (int_type k = 0; k < nruns - 1; ++k)
    {
        run_type* run = runs[k];
        run_size = run->size();
        assert(run_size == m2);
        int_type next_run_size = runs[k + 1]->size();
        STXXL_ASSERT((next_run_size == m2) || (next_run_size <= m2 && k == nruns - 2));

        STXXL_VERBOSE1("stxxl::create_runs start waiting read_reqs1");
        wait_all(read_reqs1, run_size);
        STXXL_VERBOSE1("stxxl::create_runs finish waiting read_reqs1");
        for (i = 0; i < run_size; ++i)
            bm->delete_block(bids1[i]);

        check_sort_settings();
        potentially_parallel::
        sort(make_element_iterator(Blocks1, 0),
             make_element_iterator(Blocks1, run_size * block_type::size),
             cmp);

        STXXL_VERBOSE1("stxxl::create_runs start waiting write_reqs");
        if (k > 0)
            wait_all(write_reqs, m2);
        STXXL_VERBOSE1("stxxl::create_runs finish waiting write_reqs");

        int_type runplus2size = (k < nruns - 2) ? runs[k + 2]->size() : 0;
        for (i = 0; i < m2; ++i)
        {
            STXXL_VERBOSE1("stxxl::create_runs posting write " << Blocks1[i].elem);
            (*run)[i].value = Blocks1[i][0];
            if (i >= runplus2size) {
                write_reqs[i] = Blocks1[i].write((*run)[i].bid);
            }
            else
            {
                next_run_reads[i].block = Blocks1 + i;
                next_run_reads[i].req = read_reqs1 + i;
                bids1[i] = next_run_reads[i].bid = *(it++);
                write_reqs[i] = Blocks1[i].write((*run)[i].bid, next_run_reads[i]);
            }
        }
        std::swap(Blocks1, Blocks2);
        std::swap(bids1, bids2);
        std::swap(read_reqs1, read_reqs2);
    }

    run_type* run = runs[nruns - 1];
    run_size = run->size();
    STXXL_VERBOSE1("stxxl::create_runs start waiting read_reqs1");
    wait_all(read_reqs1, run_size);
    STXXL_VERBOSE1("stxxl::create_runs finish waiting read_reqs1");
    for (i = 0; i < run_size; ++i)
        bm->delete_block(bids1[i]);

    check_sort_settings();
    potentially_parallel::
    sort(make_element_iterator(Blocks1, 0),
         make_element_iterator(Blocks1, run_size * block_type::size),
         cmp);

    STXXL_VERBOSE1("stxxl::create_runs start waiting write_reqs");
    wait_all(write_reqs, m2);
    STXXL_VERBOSE1("stxxl::create_runs finish waiting write_reqs");

    for (i = 0; i < run_size; ++i)
    {
        STXXL_VERBOSE1("stxxl::create_runs posting write " << Blocks1[i].elem);
        (*run)[i].value = Blocks1[i][0];
        write_reqs[i] = Blocks1[i].write((*run)[i].bid);
    }

    STXXL_VERBOSE1("stxxl::create_runs start waiting write_reqs");
    wait_all(write_reqs, run_size);
    STXXL_VERBOSE1("stxxl::create_runs finish waiting write_reqs");

    delete[] Blocks1;
    delete[] Blocks2;
    delete[] bids1;
    delete[] bids2;
    delete[] read_reqs1;
    delete[] read_reqs2;
    delete[] write_reqs;
    delete[] next_run_reads;
}

template <typename BlockType, typename RunType, typename ValueCmp>
bool check_sorted_runs(RunType** runs,
                       unsigned_type nruns,
                       unsigned_type m,
                       ValueCmp cmp)
{
    typedef BlockType block_type;
    typedef typename block_type::value_type value_type;

    STXXL_MSG("check_sorted_runs  Runs: " << nruns);
    unsigned_type irun = 0;
    for (irun = 0; irun < nruns; ++irun)
    {
        const unsigned_type nblocks_per_run = runs[irun]->size();
        unsigned_type blocks_left = nblocks_per_run;
        block_type* blocks = new block_type[m];
        request_ptr* reqs = new request_ptr[m];
        value_type last = cmp.min_value();

        for (unsigned_type off = 0; off < nblocks_per_run; off += m)
        {
            const unsigned_type nblocks = STXXL_MIN(blocks_left, m);
            const unsigned_type nelements = nblocks * block_type::size;
            blocks_left -= nblocks;

            for (unsigned_type j = 0; j < nblocks; ++j)
            {
                reqs[j] = blocks[j].read((*runs[irun])[j + off].bid);
            }
            wait_all(reqs, reqs + nblocks);

            if (off && cmp(blocks[0][0], last))
            {
                STXXL_MSG("check_sorted_runs  wrong first value in the run " << irun);
                STXXL_MSG(" first value: " << blocks[0][0]);
                STXXL_MSG(" last  value: " << last);
                for (unsigned_type k = 0; k < block_type::size; ++k)
                    STXXL_MSG("Element " << k << " in the block is :" << blocks[0][k]);

                delete[] reqs;
                delete[] blocks;
                return false;
            }

            for (unsigned_type j = 0; j < nblocks; ++j)
            {
                if (!(blocks[j][0] == (*runs[irun])[j + off].value))
                {
                    STXXL_MSG("check_sorted_runs  wrong trigger in the run " << irun << " block " << (j + off));
                    STXXL_MSG("                   trigger value: " << (*runs[irun])[j + off].value);
                    STXXL_MSG("Data in the block:");
                    for (unsigned_type k = 0; k < block_type::size; ++k)
                        STXXL_MSG("Element " << k << " in the block is :" << blocks[j][k]);

                    STXXL_MSG("BIDS:");
                    for (unsigned_type k = 0; k < nblocks; ++k)
                    {
                        if (k == j)
                            STXXL_MSG("Bad one comes next.");
                        STXXL_MSG("BID " << (k + off) << " is: " << ((*runs[irun])[k + off].bid));
                    }

                    delete[] reqs;
                    delete[] blocks;
                    return false;
                }
            }
            if (!stxxl::is_sorted(make_element_iterator(blocks, 0),
                                  make_element_iterator(blocks, nelements),
                                  cmp))
            {
                STXXL_MSG("check_sorted_runs  wrong order in the run " << irun);
                STXXL_MSG("Data in blocks:");
                for (unsigned_type j = 0; j < nblocks; ++j)
                {
                    for (unsigned_type k = 0; k < block_type::size; ++k)
                        STXXL_MSG("     Element " << k << " in block " << (j + off) << " is :" << blocks[j][k]);
                }
                STXXL_MSG("BIDS:");
                for (unsigned_type k = 0; k < nblocks; ++k)
                {
                    STXXL_MSG("BID " << (k + off) << " is: " << ((*runs[irun])[k + off].bid));
                }

                delete[] reqs;
                delete[] blocks;
                return false;
            }

            last = blocks[nblocks - 1][block_type::size - 1];
        }

        assert(blocks_left == 0);
        delete[] reqs;
        delete[] blocks;
    }

    return true;
}

template <typename BlockType, typename RunType, typename ValueCmp>
void merge_runs(RunType** in_runs, int_type nruns,
                RunType* out_run, unsigned_type _m, ValueCmp cmp)
{
    typedef BlockType block_type;
    typedef RunType run_type;
    typedef ValueCmp value_cmp;
    typedef typename run_type::value_type trigger_entry_type;
    typedef block_prefetcher<block_type, typename run_type::iterator> prefetcher_type;
    typedef run_cursor2<block_type, prefetcher_type> run_cursor_type;
    typedef sort_helper::run_cursor2_cmp<block_type, prefetcher_type, value_cmp> run_cursor2_cmp_type;

    run_type consume_seq(out_run->size());

    int_type* prefetch_seq = new int_type[out_run->size()];

    typename run_type::iterator copy_start = consume_seq.begin();
    for (int_type i = 0; i < nruns; i++)
    {
        // \todo: try to avoid copy
        copy_start = std::copy(
            in_runs[i]->begin(),
            in_runs[i]->end(),
            copy_start);
    }

    std::stable_sort(consume_seq.begin(), consume_seq.end(),
                     sort_helper::trigger_entry_cmp<trigger_entry_type, value_cmp>(cmp) _STXXL_SORT_TRIGGER_FORCE_SEQUENTIAL);

    int_type disks_number = config::get_instance()->disks_number();

#ifdef PLAY_WITH_OPT_PREF
    const int_type n_write_buffers = 4 * disks_number;
#else
    const int_type n_prefetch_buffers = STXXL_MAX(2 * disks_number, (3 * (int_type(_m) - nruns) / 4));
    const int_type n_write_buffers = STXXL_MAX(2 * disks_number, int_type(_m) - nruns - n_prefetch_buffers);
 #if STXXL_SORT_OPTIMAL_PREFETCHING
    // heuristic
    const int_type n_opt_prefetch_buffers = 2 * disks_number + (3 * (n_prefetch_buffers - 2 * disks_number)) / 10;
 #endif
#endif

#if STXXL_SORT_OPTIMAL_PREFETCHING
    compute_prefetch_schedule(
        consume_seq,
        prefetch_seq,
        n_opt_prefetch_buffers,
        config::get_instance()->get_max_device_id());
#else
    for (unsigned_type i = 0; i < out_run->size(); i++)
        prefetch_seq[i] = i;

#endif

    prefetcher_type prefetcher(consume_seq.begin(),
                               consume_seq.end(),
                               prefetch_seq,
                               nruns + n_prefetch_buffers);

    buffered_writer<block_type> writer(n_write_buffers, n_write_buffers / 2);

    int_type out_run_size = out_run->size();

    block_type* out_buffer = writer.get_free_block();

//If parallelism is activated, one can still fall back to the
//native merge routine by setting stxxl::SETTINGS::native_merge= true, //otherwise, it is used anyway.

    if (do_parallel_merge())
    {
#if STXXL_PARALLEL_MULTIWAY_MERGE

// begin of STL-style merging

        typedef stxxl::int64 diff_type;
        typedef std::pair<typename block_type::iterator, typename block_type::iterator> sequence;
        std::vector<sequence> seqs(nruns);
        std::vector<block_type*> buffers(nruns);

        for (int_type i = 0; i < nruns; i++)                    // initialize sequences
        {
            buffers[i] = prefetcher.pull_block();               // get first block of each run
            seqs[i] = std::make_pair(buffers[i]->begin(), buffers[i]->end());
            // this memory location stays the same, only the data is exchanged
        }

 #if STXXL_CHECK_ORDER_IN_SORTS
        value_type last_elem = cmp.min_value();
 #endif
        diff_type num_currently_mergeable = 0;

        for (int_type j = 0; j < out_run_size; ++j)                     // for the whole output run, out_run_size is in blocks
        {
            diff_type rest = block_type::size;                          // elements still to merge for this output block

            STXXL_VERBOSE1("output block " << j);
            do {
                if (num_currently_mergeable < rest)
                {
                    if (prefetcher.empty())
                    {
                        // anything remaining is already in memory
                        num_currently_mergeable = (out_run_size - j) * block_type::size
                                                  - (block_type::size - rest);
                    }
                    else
                    {
                        num_currently_mergeable = sort_helper::count_elements_less_equal(
                            seqs, consume_seq[prefetcher.pos()].value, cmp);
                    }
                }

                diff_type output_size = STXXL_MIN(num_currently_mergeable, rest);       // at most rest elements

                STXXL_VERBOSE1("before merge " << output_size);

                parallel::multiway_merge(
                    seqs.begin(), seqs.end(),
                    out_buffer->end() - rest, output_size, cmp);
                // sequence iterators are progressed appropriately

                rest -= output_size;
                num_currently_mergeable -= output_size;

                STXXL_VERBOSE1("after merge");

                sort_helper::refill_or_remove_empty_sequences(seqs, buffers, prefetcher);
            } while (rest > 0 && seqs.size() > 0);

 #if STXXL_CHECK_ORDER_IN_SORTS
            if (!stxxl::is_sorted(out_buffer->begin(), out_buffer->end(), cmp))
            {
                for (value_type* i = out_buffer->begin() + 1; i != out_buffer->end(); i++)
                    if (cmp(*i, *(i - 1)))
                    {
                        STXXL_VERBOSE1("Error at position " << (i - out_buffer->begin()));
                    }
                assert(false);
            }

            if (j > 0)     // do not check in first iteration
                assert(cmp((*out_buffer)[0], last_elem) == false);

            last_elem = (*out_buffer)[block_type::size - 1];
 #endif

            (*out_run)[j].value = (*out_buffer)[0];                                  // save smallest value

            out_buffer = writer.write(out_buffer, (*out_run)[j].bid);
        }

// end of STL-style merging

#else
        STXXL_THROW_UNREACHABLE();
#endif
    }
    else
    {
// begin of native merging procedure

        loser_tree<run_cursor_type, run_cursor2_cmp_type>
        losers(&prefetcher, nruns, run_cursor2_cmp_type(cmp));

#if STXXL_CHECK_ORDER_IN_SORTS
        value_type last_elem = cmp.min_value();
#endif

        for (int_type i = 0; i < out_run_size; ++i)
        {
            losers.multi_merge(out_buffer->elem, out_buffer->elem + block_type::size);
            (*out_run)[i].value = *(out_buffer->elem);

#if STXXL_CHECK_ORDER_IN_SORTS
            assert(stxxl::is_sorted(
                       out_buffer->begin(),
                       out_buffer->end(),
                       cmp));

            if (i)
                assert(cmp(*(out_buffer->elem), last_elem) == false);

            last_elem = (*out_buffer).elem[block_type::size - 1];
#endif

            out_buffer = writer.write(out_buffer, (*out_run)[i].bid);
        }

// end of native merging procedure
    }

    delete[] prefetch_seq;

    block_manager* bm = block_manager::get_instance();
    for (int_type i = 0; i < nruns; ++i)
    {
        unsigned_type sz = in_runs[i]->size();
        for (unsigned_type j = 0; j < sz; ++j)
            bm->delete_block((*in_runs[i])[j].bid);

        delete in_runs[i];
    }
}

template <typename BlockType,
          typename AllocStrategy,
          typename InputBidIterator,
          typename ValueCmp>
simple_vector<sort_helper::trigger_entry<BlockType> >*
sort_blocks(InputBidIterator input_bids,
            unsigned_type _n,
            unsigned_type _m,
            ValueCmp cmp)
{
    typedef BlockType block_type;
    typedef AllocStrategy alloc_strategy;
    typedef InputBidIterator input_bid_iterator;
    typedef ValueCmp value_cmp;
    typedef typename block_type::bid_type bid_type;
    typedef sort_helper::trigger_entry<block_type> trigger_entry_type;
    typedef simple_vector<trigger_entry_type> run_type;
    typedef typename interleaved_alloc_traits<alloc_strategy>::strategy interleaved_alloc_strategy;

    unsigned_type m2 = _m / 2;
    unsigned_type full_runs = _n / m2;
    unsigned_type partial_runs = ((_n % m2) ? 1 : 0);
    unsigned_type nruns = full_runs + partial_runs;
    unsigned_type i;

    block_manager* mng = block_manager::get_instance();

    //STXXL_VERBOSE ("n=" << _n << " nruns=" << nruns << "=" << full_runs << "+" << partial_runs);

    double begin = timestamp(), after_runs_creation, end;

    run_type** runs = new run_type*[nruns];

    for (i = 0; i < full_runs; i++)
        runs[i] = new run_type(m2);

    if (partial_runs)
        runs[i] = new run_type(_n - full_runs * m2);

    for (i = 0; i < nruns; ++i)
        mng->new_blocks(alloc_strategy(), make_bid_iterator(runs[i]->begin()), make_bid_iterator(runs[i]->end()));

    sort_local::create_runs<block_type,
                            run_type,
                            input_bid_iterator,
                            value_cmp>(input_bids, runs, nruns, _m, cmp);

    after_runs_creation = timestamp();

    double io_wait_after_rf = stats::get_instance()->get_io_wait_time();

    disk_queues::get_instance()->set_priority_op(request_queue::WRITE);

    const int_type merge_factor = optimal_merge_factor(nruns, _m);
    run_type** new_runs;

    while (nruns > 1)
    {
        int_type new_nruns = div_ceil(nruns, merge_factor);
        STXXL_VERBOSE("Starting new merge phase: nruns: " << nruns <<
                      " opt_merge_factor: " << merge_factor << " m:" << _m << " new_nruns: " << new_nruns);

        new_runs = new run_type*[new_nruns];

        int_type runs_left = nruns;
        int_type cur_out_run = 0;
        int_type blocks_in_new_run = 0;

        while (runs_left > 0)
        {
            int_type runs2merge = STXXL_MIN(runs_left, merge_factor);
            blocks_in_new_run = 0;
            for (unsigned_type i = nruns - runs_left; i < (nruns - runs_left + runs2merge); i++)
                blocks_in_new_run += runs[i]->size();

            // allocate run
            new_runs[cur_out_run++] = new run_type(blocks_in_new_run);
            runs_left -= runs2merge;
        }
        // allocate blocks for the new runs
        if (cur_out_run == 1 && blocks_in_new_run == int_type(_n) && !input_bids->is_managed())
        {
            // if we sort a file we can reuse the input bids for the output
            input_bid_iterator cur = input_bids;
            for (int_type i = 0; cur != (input_bids + _n); ++cur)
            {
                (*new_runs[0])[i++].bid = *cur;
            }

            bid_type& firstBID = (*new_runs[0])[0].bid;
            if (firstBID.is_managed())
            {
                // the first block does not belong to the file
                // need to reallocate it
                mng->new_block(FR(), firstBID);
            }
            bid_type& lastBID = (*new_runs[0])[_n - 1].bid;
            if (lastBID.is_managed())
            {
                // the first block does not belong to the file
                // need to reallocate it
                mng->new_block(FR(), lastBID);
            }
        }
        else
        {
            mng->new_blocks(interleaved_alloc_strategy(new_nruns, alloc_strategy()),
                            runs2bid_array_adaptor2<block_type::raw_size, run_type>(new_runs, 0, new_nruns, blocks_in_new_run),
                            runs2bid_array_adaptor2<block_type::raw_size, run_type>(new_runs, _n, new_nruns, blocks_in_new_run));
        }
        // merge all
        runs_left = nruns;
        cur_out_run = 0;
        while (runs_left > 0)
        {
            int_type runs2merge = STXXL_MIN(runs_left, merge_factor);
#if STXXL_CHECK_ORDER_IN_SORTS
            assert((check_sorted_runs<block_type, run_type, value_cmp>(runs + nruns - runs_left, runs2merge, m2, cmp)));
#endif
            STXXL_VERBOSE("Merging " << runs2merge << " runs");
            merge_runs<block_type, run_type>(runs + nruns - runs_left,
                                             runs2merge, *(new_runs + (cur_out_run++)), _m, cmp
                                             );
            runs_left -= runs2merge;
        }

        nruns = new_nruns;
        delete[] runs;
        runs = new_runs;
    }

    run_type* result = *runs;
    delete[] runs;

    end = timestamp();

    STXXL_VERBOSE("Elapsed time        : " << end - begin << " s. Run creation time: " <<
                  after_runs_creation - begin << " s");
    STXXL_VERBOSE("Time in I/O wait(rf): " << io_wait_after_rf << " s");
    STXXL_VERBOSE(*stats::get_instance());

    return result;
}

} // namespace sort_local

/*!
 * Sort records comparison-based, see \ref design_algo_sort.
 *
 * stxxl::sort sorts the elements in [first, last) into ascending order,
 * meaning that if \c i and \c j are any two valid iterators in [first, last)
 * such that \c i precedes \c j, then \c *j is not less than \c *i. Note: as
 * std::sort, stxxl::sort is not guaranteed to be stable. That is, suppose that
 * \c *i and \c *j are equivalent: neither one is less than the other. It is
 * not guaranteed that the relative order of these two elements will be
 * preserved by stxxl::sort.
 *
 * The order is defined by the \c cmp parameter. The sorter's internal memory
 * consumption is bounded by \a M bytes.
 *
 * \param first object of model of \c ext_random_access_iterator concept
 * \param last object of model of \c ext_random_access_iterator concept
 * \param cmp comparison object of \ref StrictWeakOrdering
 * \param M amount of memory for internal use (in bytes)
 */
template <typename ExtIterator, typename StrictWeakOrdering>
void sort(ExtIterator first, ExtIterator last, StrictWeakOrdering cmp, unsigned_type M)
{
    sort_helper::verify_sentinel_strict_weak_ordering(cmp);

    typedef typename ExtIterator::vector_type::value_type value_type;
    typedef typename ExtIterator::block_type block_type;
    typedef typename ExtIterator::bid_type bid_type;
    typedef typename ExtIterator::vector_type::alloc_strategy_type alloc_strategy_type;
    typedef typename ExtIterator::bids_container_iterator bids_container_iterator;

    typedef simple_vector<sort_helper::trigger_entry<block_type> > run_type;

    unsigned_type n = 0;

    block_manager* mng = block_manager::get_instance();

    first.flush();

    if ((last - first) * sizeof(value_type) * sort_memory_usage_factor() < M)
    {
        stl_in_memory_sort(first, last, cmp);
    }
    else
    {
        if (!(2 * block_type::raw_size * (unsigned_type)sort_memory_usage_factor() <= M)) {
            throw bad_parameter("stxxl::sort(): INSUFFICIENT MEMORY provided, please increase parameter 'M'");
        }

        if (first.block_offset())
        {
            if (last.block_offset())              // first and last element are
            // not the first elements of their block
            {
                block_type* first_block = new block_type;
                block_type* last_block = new block_type;
                typename ExtIterator::bid_type first_bid, last_bid;
                request_ptr req;

                req = first_block->read(*first.bid());
                mng->new_block(FR(), first_bid);                // try to overlap
                mng->new_block(FR(), last_bid);
                req->wait();

                req = last_block->read(*last.bid());

                unsigned_type i = 0;
                for ( ; i < first.block_offset(); ++i)
                {
                    first_block->elem[i] = cmp.min_value();
                }

                req->wait();

                req = first_block->write(first_bid);
                for (i = last.block_offset(); i < block_type::size; ++i)
                {
                    last_block->elem[i] = cmp.max_value();
                }

                req->wait();

                req = last_block->write(last_bid);

                n = last.bid() - first.bid() + 1;

                std::swap(first_bid, *first.bid());
                std::swap(last_bid, *last.bid());

                req->wait();

                delete first_block;
                delete last_block;

                run_type* out =
                    sort_local::sort_blocks<
                        block_type, alloc_strategy_type, bids_container_iterator
                        >(first.bid(), n,
                          M / sort_memory_usage_factor() / block_type::raw_size, cmp);

                first_block = new block_type;
                last_block = new block_type;
                block_type* sorted_first_block = new block_type;
                block_type* sorted_last_block = new block_type;
                request_ptr* reqs = new request_ptr[2];

                reqs[0] = first_block->read(first_bid);
                reqs[1] = sorted_first_block->read((*(out->begin())).bid);

                reqs[0]->wait();
                reqs[1]->wait();

                reqs[0] = last_block->read(last_bid);
                reqs[1] = sorted_last_block->read(((*out)[out->size() - 1]).bid);

                for (i = first.block_offset(); i < block_type::size; i++)
                {
                    first_block->elem[i] = sorted_first_block->elem[i];
                }

                reqs[0]->wait();
                reqs[1]->wait();

                req = first_block->write(first_bid);

                for (i = 0; i < last.block_offset(); ++i)
                {
                    last_block->elem[i] = sorted_last_block->elem[i];
                }

                req->wait();

                req = last_block->write(last_bid);

                mng->delete_block(out->begin()->bid);
                mng->delete_block((*out)[out->size() - 1].bid);

                *first.bid() = first_bid;
                *last.bid() = last_bid;

                typename run_type::iterator it = out->begin();
                ++it;
                bids_container_iterator cur_bid = first.bid();
                ++cur_bid;

                for ( ; cur_bid != last.bid(); ++cur_bid, ++it)
                {
                    *cur_bid = (*it).bid;
                }

                delete first_block;
                delete sorted_first_block;
                delete sorted_last_block;
                delete[] reqs;
                delete out;

                req->wait();

                delete last_block;
            }
            else
            {
                // first element is
                // not the first element of its block
                block_type* first_block = new block_type;
                bid_type first_bid;
                request_ptr req;

                req = first_block->read(*first.bid());
                mng->new_block(FR(), first_bid);                // try to overlap
                req->wait();

                unsigned_type i = 0;
                for ( ; i < first.block_offset(); ++i)
                {
                    first_block->elem[i] = cmp.min_value();
                }

                req = first_block->write(first_bid);

                n = last.bid() - first.bid();

                std::swap(first_bid, *first.bid());

                req->wait();

                delete first_block;

                run_type* out =
                    sort_local::sort_blocks<
                        block_type, alloc_strategy_type, bids_container_iterator
                        >(first.bid(), n,
                          M / sort_memory_usage_factor() / block_type::raw_size, cmp);

                first_block = new block_type;

                block_type* sorted_first_block = new block_type;

                request_ptr* reqs = new request_ptr[2];

                reqs[0] = first_block->read(first_bid);
                reqs[1] = sorted_first_block->read((*(out->begin())).bid);

                reqs[0]->wait();
                reqs[1]->wait();

                for (i = first.block_offset(); i < block_type::size; ++i)
                {
                    first_block->elem[i] = sorted_first_block->elem[i];
                }

                req = first_block->write(first_bid);

                mng->delete_block(out->begin()->bid);

                *first.bid() = first_bid;

                typename run_type::iterator it = out->begin();
                ++it;
                bids_container_iterator cur_bid = first.bid();
                ++cur_bid;

                for ( ; cur_bid != last.bid(); ++cur_bid, ++it)
                {
                    *cur_bid = (*it).bid;
                }

                *cur_bid = (*it).bid;

                delete sorted_first_block;
                delete[] reqs;
                delete out;

                req->wait();

                delete first_block;
            }
        }
        else
        {
            if (last.block_offset())            // last is
            // not the first element of its block
            {
                block_type* last_block = new block_type;
                bid_type last_bid;
                request_ptr req;
                unsigned_type i;

                req = last_block->read(*last.bid());
                mng->new_block(FR(), last_bid);
                req->wait();

                for (i = last.block_offset(); i < block_type::size; ++i)
                {
                    last_block->elem[i] = cmp.max_value();
                }

                req = last_block->write(last_bid);

                n = last.bid() - first.bid() + 1;

                std::swap(last_bid, *last.bid());

                req->wait();

                delete last_block;

                run_type* out =
                    sort_local::sort_blocks<
                        block_type, alloc_strategy_type, bids_container_iterator
                        >(first.bid(), n,
                          M / sort_memory_usage_factor() / block_type::raw_size, cmp);

                last_block = new block_type;
                block_type* sorted_last_block = new block_type;
                request_ptr* reqs = new request_ptr[2];

                reqs[0] = last_block->read(last_bid);
                reqs[1] = sorted_last_block->read(((*out)[out->size() - 1]).bid);

                reqs[0]->wait();
                reqs[1]->wait();

                for (i = 0; i < last.block_offset(); ++i)
                {
                    last_block->elem[i] = sorted_last_block->elem[i];
                }

                req = last_block->write(last_bid);

                mng->delete_block((*out)[out->size() - 1].bid);

                *last.bid() = last_bid;

                typename run_type::iterator it = out->begin();
                bids_container_iterator cur_bid = first.bid();

                for ( ; cur_bid != last.bid(); ++cur_bid, ++it)
                {
                    *cur_bid = (*it).bid;
                }

                delete sorted_last_block;
                delete[] reqs;
                delete out;

                req->wait();

                delete last_block;
            }
            else
            {
                // first and last element are first elements of their of blocks
                n = last.bid() - first.bid();

                run_type* out =
                    sort_local::sort_blocks<
                        block_type, alloc_strategy_type, bids_container_iterator
                        >(first.bid(), n,
                          M / sort_memory_usage_factor() / block_type::raw_size, cmp);

                typename run_type::iterator it = out->begin();
                bids_container_iterator cur_bid = first.bid();

                for ( ; cur_bid != last.bid(); ++cur_bid, ++it)
                {
                    *cur_bid = (*it).bid;
                }

                delete out;
            }
        }
    }

#if STXXL_CHECK_ORDER_IN_SORTS
    typedef typename ExtIterator::const_iterator const_iterator;
    assert(stxxl::is_sorted(const_iterator(first), const_iterator(last), cmp));
#endif
}

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_SORT_HEADER
// vim: et:ts=4:sw=4
