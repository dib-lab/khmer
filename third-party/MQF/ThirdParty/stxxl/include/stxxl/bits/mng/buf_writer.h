/***************************************************************************
 *  include/stxxl/bits/mng/buf_writer.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_BUF_WRITER_HEADER
#define STXXL_MNG_BUF_WRITER_HEADER

#include <vector>
#include <queue>

#include <stxxl/bits/io/request_operations.h>
#include <stxxl/bits/io/disk_queues.h>
#include <stxxl/bits/noncopyable.h>

STXXL_BEGIN_NAMESPACE

//! \defgroup schedlayer Block Scheduling Sublayer
//! \ingroup mnglayer
//! Group of classes which help in scheduling
//! sequences of read and write requests
//! via prefetching and buffered writing
//! \{

//! Encapsulates asynchronous buffered block writing engine.
//!
//! \c buffered_writer overlaps I/Os with filling of output buffer.
template <typename BlockType>
class buffered_writer : private noncopyable
{
    typedef BlockType block_type;
    typedef typename block_type::bid_type bid_type;

protected:
    const unsigned_type nwriteblocks;
    block_type* write_buffers;
    bid_type* write_bids;
    request_ptr* write_reqs;
    const unsigned_type writebatchsize;

    std::vector<int_type> free_write_blocks;            // contains free write blocks
    std::vector<int_type> busy_write_blocks;            // blocks that are in writing, notice that if block is not in free_
    // an not in busy then block is not yet filled

    struct batch_entry
    {
        stxxl::int64 offset;
        int_type ibuffer;
        batch_entry(stxxl::int64 o, int_type b) : offset(o), ibuffer(b) { }
    };
    struct batch_entry_cmp
    {
        bool operator () (const batch_entry& a, const batch_entry& b) const
        {
            return (a.offset > b.offset);
        }
    };

    typedef std::priority_queue<batch_entry, std::vector<batch_entry>, batch_entry_cmp> batch_type;
    batch_type batch_write_blocks;      // sorted sequence of blocks to write

public:
    //! Constructs an object.
    //! \param write_buf_size number of write buffers to use
    //! \param write_batch_size number of blocks to accumulate in
    //!        order to flush write requests (bulk buffered writing)
    buffered_writer(unsigned_type write_buf_size, unsigned_type write_batch_size)
        : nwriteblocks((write_buf_size > 2) ? write_buf_size : 2),
          writebatchsize(write_batch_size ? write_batch_size : 1)
    {
        write_buffers = new block_type[nwriteblocks];
        write_reqs = new request_ptr[nwriteblocks];

        write_bids = new bid_type[nwriteblocks];

        for (unsigned_type i = 0; i < nwriteblocks; i++)
            free_write_blocks.push_back(i);

        disk_queues::get_instance()->set_priority_op(request_queue::WRITE);
    }
    //! Returns free block from the internal buffer pool.
    //! \return pointer to the block from the internal buffer pool
    block_type * get_free_block()
    {
        int_type ibuffer;
        for (std::vector<int_type>::iterator it = busy_write_blocks.begin();
             it != busy_write_blocks.end(); ++it)
        {
            if (write_reqs[ibuffer = (*it)]->poll())
            {
                busy_write_blocks.erase(it);
                free_write_blocks.push_back(ibuffer);

                break;
            }
        }
        if (UNLIKELY(free_write_blocks.empty()))
        {
            int_type size = busy_write_blocks.size();
            request_ptr* reqs = new request_ptr[size];
            int_type i = 0;
            for ( ; i < size; ++i)
            {
                reqs[i] = write_reqs[busy_write_blocks[i]];
            }
            int_type completed = wait_any(reqs, size);
            int_type completed_global = busy_write_blocks[completed];
            delete[] reqs;
            busy_write_blocks.erase(busy_write_blocks.begin() + completed);

            return (write_buffers + completed_global);
        }
        ibuffer = free_write_blocks.back();
        free_write_blocks.pop_back();

        return (write_buffers + ibuffer);
    }
    //! Submits block for writing.
    //! \param filled_block pointer to the block
    //! \remark parameter \c filled_block must be value returned by \c get_free_block() or \c write() methods
    //! \param bid block identifier, a place to write data of the \c filled_block
    //! \return pointer to the new free block from the pool
    block_type * write(block_type* filled_block, const bid_type& bid)          // writes filled_block and returns a new block
    {
        if (batch_write_blocks.size() >= writebatchsize)
        {
            // flush batch
            while (!batch_write_blocks.empty())
            {
                int_type ibuffer = batch_write_blocks.top().ibuffer;
                batch_write_blocks.pop();

                if (write_reqs[ibuffer].valid())
                    write_reqs[ibuffer]->wait();

                write_reqs[ibuffer] = write_buffers[ibuffer].write(write_bids[ibuffer]);

                busy_write_blocks.push_back(ibuffer);
            }
        }
        //    STXXL_MSG("Adding write request to batch");

        int_type ibuffer = filled_block - write_buffers;
        write_bids[ibuffer] = bid;
        batch_write_blocks.push(batch_entry(bid.offset, ibuffer));

        return get_free_block();
    }
    //! Flushes not yet written buffers.
    void flush()
    {
        int_type ibuffer;
        while (!batch_write_blocks.empty())
        {
            ibuffer = batch_write_blocks.top().ibuffer;
            batch_write_blocks.pop();

            if (write_reqs[ibuffer].valid())
                write_reqs[ibuffer]->wait();

            write_reqs[ibuffer] = write_buffers[ibuffer].write(write_bids[ibuffer]);

            busy_write_blocks.push_back(ibuffer);
        }
        for (std::vector<int_type>::const_iterator it =
                 busy_write_blocks.begin();
             it != busy_write_blocks.end(); it++)
        {
            ibuffer = *it;
            write_reqs[ibuffer]->wait();
        }

        assert(batch_write_blocks.empty());
        free_write_blocks.clear();
        busy_write_blocks.clear();

        for (unsigned_type i = 0; i < nwriteblocks; i++)
            free_write_blocks.push_back(i);
    }

    //! Flushes not yet written buffers and frees used memory.
    ~buffered_writer()
    {
        int_type ibuffer;
        while (!batch_write_blocks.empty())
        {
            ibuffer = batch_write_blocks.top().ibuffer;
            batch_write_blocks.pop();

            if (write_reqs[ibuffer].valid())
                write_reqs[ibuffer]->wait();

            write_reqs[ibuffer] = write_buffers[ibuffer].write(write_bids[ibuffer]);

            busy_write_blocks.push_back(ibuffer);
        }
        for (std::vector<int_type>::const_iterator it =
                 busy_write_blocks.begin();
             it != busy_write_blocks.end(); it++)
        {
            ibuffer = *it;
            write_reqs[ibuffer]->wait();
        }

        delete[] write_reqs;
        delete[] write_buffers;
        delete[] write_bids;
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_MNG_BUF_WRITER_HEADER
