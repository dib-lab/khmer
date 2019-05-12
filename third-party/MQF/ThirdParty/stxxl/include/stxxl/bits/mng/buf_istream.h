/***************************************************************************
 *  include/stxxl/bits/mng/buf_istream.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_BUF_ISTREAM_HEADER
#define STXXL_MNG_BUF_ISTREAM_HEADER

#include <stxxl/bits/mng/config.h>
#include <stxxl/bits/mng/block_prefetcher.h>
#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/algo/async_schedule.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup schedlayer
//! \{

// a paranoid check
#define BUF_ISTREAM_CHECK_END

//! Buffered input stream.
//!
//! Reads data records from the stream of blocks.
//! \remark Reading performed in the background, i.e. with overlapping of I/O and computation
template <typename BlockType, typename BidIteratorType>
class buf_istream : private noncopyable
{
public:
    typedef BlockType block_type;
    typedef BidIteratorType bid_iterator_type;

private:
    buf_istream() { }

protected:
    typedef block_prefetcher<block_type, bid_iterator_type> prefetcher_type;
    prefetcher_type* prefetcher;
    int_type current_elem;
    block_type* current_blk;
    int_type* prefetch_seq;
#ifdef BUF_ISTREAM_CHECK_END
    bool not_finished;
#endif

public:
    typedef typename block_type::reference reference;
    typedef buf_istream<block_type, bid_iterator_type> self_type;

    //! Constructs input stream object.
    //! \param begin \c bid_iterator pointing to the first block of the stream
    //! \param end \c bid_iterator pointing to the ( \b last + 1 ) block of the stream
    //! \param nbuffers number of buffers for internal use
    buf_istream(bid_iterator_type begin, bid_iterator_type end, unsigned_type nbuffers)
        : current_elem(0)
#ifdef BUF_ISTREAM_CHECK_END
          , not_finished(true)
#endif
    {
        const unsigned_type ndisks = config::get_instance()->disks_number();
        const unsigned_type mdevid = config::get_instance()->get_max_device_id();
        const int_type seq_length = end - begin;
        prefetch_seq = new int_type[seq_length];

        // obvious schedule
        //for(int_type i = 0; i < seq_length; ++i)
        //      prefetch_seq[i] = i;

        // optimal schedule
        nbuffers = STXXL_MAX(2 * ndisks, unsigned_type(nbuffers - 1));
        compute_prefetch_schedule(begin, end, prefetch_seq,
                                  nbuffers, mdevid);

        prefetcher = new prefetcher_type(begin, end, prefetch_seq, nbuffers);

        current_blk = prefetcher->pull_block();
    }

    //! Input stream operator, reads in \c record.
    //! \param record reference to the block record type,
    //!        contains value of the next record in the stream after the call of the operator
    //! \return reference to itself (stream object)
    self_type& operator >> (reference record)
    {
#ifdef BUF_ISTREAM_CHECK_END
        assert(not_finished);
#endif

        record = current_blk->elem[current_elem++];

        if (UNLIKELY(current_elem >= block_type::size))
        {
            current_elem = 0;
#ifdef BUF_ISTREAM_CHECK_END
            not_finished = prefetcher->block_consumed(current_blk);
#else
            prefetcher->block_consumed(current_blk);
#endif
        }

        return (*this);
    }

    //! Returns reference to the current record in the stream.
    reference current()     /* const */
    {
        return current_blk->elem[current_elem];
    }

    //! Returns reference to the current record in the stream.
    reference operator * ()     /* const */
    {
        return current_blk->elem[current_elem];
    }

    //! Moves to the next record in the stream.
    //! \return reference to itself after the advance
    self_type& operator ++ ()
    {
#ifdef BUF_ISTREAM_CHECK_END
        assert(not_finished);
#endif

        current_elem++;

        if (UNLIKELY(current_elem >= block_type::size))
        {
            current_elem = 0;
#ifdef BUF_ISTREAM_CHECK_END
            not_finished = prefetcher->block_consumed(current_blk);
#else
            prefetcher->block_consumed(current_blk);
#endif
        }
        return *this;
    }

    //! Frees used internal objects.
    ~buf_istream()
    {
        delete prefetcher;
        delete[] prefetch_seq;
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_MNG_BUF_ISTREAM_HEADER
