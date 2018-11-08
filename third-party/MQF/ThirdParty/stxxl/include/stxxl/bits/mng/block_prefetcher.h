/***************************************************************************
 *  include/stxxl/bits/mng/block_prefetcher.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009, 2010 Johannes Singler <singler@kit.edu>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_BLOCK_PREFETCHER_HEADER
#define STXXL_MNG_BLOCK_PREFETCHER_HEADER

#include <vector>
#include <queue>

#include <stxxl/bits/common/onoff_switch.h>
#include <stxxl/bits/io/request.h>
#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/noncopyable.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup schedlayer
//! \{

class set_switch_handler
{
    onoff_switch& switch_;
    completion_handler on_compl;

public:
    set_switch_handler(onoff_switch& _switch, const completion_handler& on_compl)
        : switch_(_switch), on_compl(on_compl)
    { }

    void operator () (request* req)
    {
        // call before setting switch to on, otherwise, user has no way to wait
        // for the completion handler to be executed
        on_compl(req);
        switch_.on();
    }
};

//! Encapsulates asynchronous prefetching engine.
//!
//! \c block_prefetcher overlaps I/Os with consumption of read data.
//! Utilizes optimal asynchronous prefetch scheduling (by Peter Sanders et.al.)
template <typename BlockType, typename BidIteratorType>
class block_prefetcher : private noncopyable
{
public:
    typedef BlockType block_type;
    typedef BidIteratorType bid_iterator_type;

    typedef typename block_type::bid_type bid_type;

protected:
    bid_iterator_type consume_seq_begin;
    bid_iterator_type consume_seq_end;
    unsigned_type seq_length;

    int_type* prefetch_seq;

    unsigned_type nextread;
    unsigned_type nextconsume;

    const int_type nreadblocks;

    block_type* read_buffers;
    request_ptr* read_reqs;
    bid_type* read_bids;

    onoff_switch* completed;
    int_type* pref_buffer;

    completion_handler do_after_fetch;

    block_type * wait(int_type iblock)
    {
        STXXL_VERBOSE1("block_prefetcher: waiting block " << iblock);
        {
            stats::scoped_wait_timer wait_timer(stats::WAIT_OP_READ);

            completed[iblock].wait_for_on();
        }
        STXXL_VERBOSE1("block_prefetcher: finished waiting block " << iblock);
        int_type ibuffer = pref_buffer[iblock];
        STXXL_VERBOSE1("block_prefetcher: returning buffer " << ibuffer);
        assert(ibuffer >= 0 && ibuffer < nreadblocks);
        return (read_buffers + ibuffer);
    }

public:
    //! Constructs an object and immediately starts prefetching.
    //! \param _cons_begin \c bid_iterator pointing to the \c bid of the first block to be consumed
    //! \param _cons_end \c bid_iterator pointing to the \c bid of the ( \b last + 1 ) block of consumption sequence
    //! \param _pref_seq gives the prefetch order, is a pointer to the integer array that contains
    //!        the indices of the blocks in the consumption sequence
    //! \param _prefetch_buf_size amount of prefetch buffers to use
    //! \param do_after_fetch unknown
    block_prefetcher(
        bid_iterator_type _cons_begin,
        bid_iterator_type _cons_end,
        int_type* _pref_seq,
        int_type _prefetch_buf_size,
        completion_handler do_after_fetch = completion_handler())
        : consume_seq_begin(_cons_begin),
          consume_seq_end(_cons_end),
          seq_length(_cons_end - _cons_begin),
          prefetch_seq(_pref_seq),
          nextread(STXXL_MIN(unsigned_type(_prefetch_buf_size), seq_length)),
          nextconsume(0),
          nreadblocks(nextread),
          do_after_fetch(do_after_fetch)
    {
        STXXL_VERBOSE1("block_prefetcher: seq_length=" << seq_length);
        STXXL_VERBOSE1("block_prefetcher: _prefetch_buf_size=" << _prefetch_buf_size);
        assert(seq_length > 0);
        assert(_prefetch_buf_size > 0);
        int_type i;
        read_buffers = new block_type[nreadblocks];
        read_reqs = new request_ptr[nreadblocks];
        read_bids = new bid_type[nreadblocks];
        pref_buffer = new int_type[seq_length];

        std::fill(pref_buffer, pref_buffer + seq_length, -1);

        completed = new onoff_switch[seq_length];

        for (i = 0; i < nreadblocks; ++i)
        {
            assert(prefetch_seq[i] < int_type(seq_length));
            assert(prefetch_seq[i] >= 0);
            read_bids[i] = *(consume_seq_begin + prefetch_seq[i]);
            STXXL_VERBOSE1("block_prefetcher: reading block " << i <<
                           " prefetch_seq[" << i << "]=" << prefetch_seq[i] <<
                           " @ " << &read_buffers[i] <<
                           " @ " << read_bids[i]);
            read_reqs[i] = read_buffers[i].read(
                read_bids[i],
                set_switch_handler(*(completed + prefetch_seq[i]), do_after_fetch));
            pref_buffer[prefetch_seq[i]] = i;
        }
    }
    //! Pulls next unconsumed block from the consumption sequence.
    //! \return Pointer to the already prefetched block from the internal buffer pool
    block_type * pull_block()
    {
        STXXL_VERBOSE1("block_prefetcher: pulling a block");
        return wait(nextconsume++);
    }
    //! Exchanges buffers between prefetcher and application.
    //! \param buffer pointer to the consumed buffer. After call if return value is true \c buffer
    //!        contains valid pointer to the next unconsumed prefetched buffer.
    //! \remark parameter \c buffer must be value returned by \c pull_block() or \c block_consumed() methods
    //! \return \c false if there are no blocks to prefetch left, \c true if consumption sequence is not emptied
    bool block_consumed(block_type*& buffer)
    {
        int_type ibuffer = buffer - read_buffers;
        STXXL_VERBOSE1("block_prefetcher: buffer " << ibuffer << " consumed");
        if (read_reqs[ibuffer].valid())
            read_reqs[ibuffer]->wait();

        read_reqs[ibuffer] = NULL;

        if (nextread < seq_length)
        {
            assert(ibuffer >= 0 && ibuffer < nreadblocks);
            int_type next_2_prefetch = prefetch_seq[nextread++];
            STXXL_VERBOSE1("block_prefetcher: prefetching block " << next_2_prefetch);

            assert((next_2_prefetch < int_type(seq_length)) && (next_2_prefetch >= 0));
            assert(!completed[next_2_prefetch].is_on());

            pref_buffer[next_2_prefetch] = ibuffer;
            read_bids[ibuffer] = *(consume_seq_begin + next_2_prefetch);
            read_reqs[ibuffer] = read_buffers[ibuffer].read(
                read_bids[ibuffer],
                set_switch_handler(*(completed + next_2_prefetch), do_after_fetch)
                );
        }

        if (nextconsume >= seq_length)
            return false;

        buffer = wait(nextconsume++);

        return true;
    }

    //! No more consumable blocks available, but can't delete the prefetcher,
    //! because not all blocks may have been returned, yet.
    bool empty() const
    {
        return nextconsume >= seq_length;
    }

    //! Index of the next element in the consume sequence.
    unsigned_type pos() const
    {
        return nextconsume;
    }

    //! Frees used memory.
    ~block_prefetcher()
    {
        for (int_type i = 0; i < nreadblocks; ++i)
            if (read_reqs[i].valid())
                read_reqs[i]->wait();

        delete[] read_reqs;
        delete[] read_bids;
        delete[] completed;
        delete[] pref_buffer;
        delete[] read_buffers;
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_MNG_BLOCK_PREFETCHER_HEADER
