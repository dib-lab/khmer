/***************************************************************************
 *  include/stxxl/bits/mng/buf_ostream.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MNG_BUF_OSTREAM_HEADER
#define STXXL_MNG_BUF_OSTREAM_HEADER

#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/mng/buf_writer.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup schedlayer
//! \{

//! Buffered output stream.
//!
//! Writes data records to the stream of blocks.
//! \remark Writing performed in the background, i.e. with overlapping of I/O and computation
template <typename BlockType, typename BidIteratorType>
class buf_ostream : private noncopyable
{
public:
    typedef BlockType block_type;
    typedef BidIteratorType bid_iterator_type;

protected:
    buffered_writer<block_type> writer;
    bid_iterator_type current_bid;
    int_type current_elem;
    block_type* current_blk;

public:
    typedef typename block_type::const_reference const_reference;
    typedef typename block_type::reference reference;
    typedef buf_ostream<block_type, bid_iterator_type> self_type;

    //! Constructs output stream object.
    //! \param first_bid \c bid_iterator pointing to the first block of the stream
    //! \param nbuffers number of buffers for internal use
    buf_ostream(bid_iterator_type first_bid, int_type nbuffers)
        : writer(nbuffers, nbuffers / 2), current_bid(first_bid),
          current_elem(0)
    {
        current_blk = writer.get_free_block();
    }

    //! Output stream operator, writes out \c record.
    //! \param record const reference to block record type, containing a value of record to write to the stream
    //! \return reference to itself (stream object)
    self_type& operator << (const_reference record)
    {
        current_blk->elem[current_elem++] = record;
        if (UNLIKELY(current_elem >= block_type::size))
        {
            current_elem = 0;
            current_blk = writer.write(current_blk, *(current_bid++));
        }
        return *this;
    }

    //! Returns reference to the current record.
    //! \return reference to the current record
    reference current()
    {
        return current_blk->elem[current_elem];
    }

    //! Returns reference to the current record.
    //! \return reference to the current record
    reference operator * ()
    {
        return current_blk->elem[current_elem];
    }

    //! Moves to the next record in the stream.
    //! \return reference to itself after the advance
    self_type& operator ++ ()
    {
        ++current_elem;
        if (UNLIKELY(current_elem >= block_type::size))
        {
            current_elem = 0;
            current_blk = writer.write(current_blk, *(current_bid++));
        }
        return *this;
    }

    //! Fill current block with padding and flush
    self_type & fill(const_reference record)
    {
        while (current_elem != 0)
        {
            operator << (record);
        }
        return *this;
    }

    //! Force flush of current block, for finishing writing within a block.
    //! \warning Use with caution as the block may contain uninitialized data
    self_type & flush()
    {
        current_elem = 0;
        current_blk = writer.write(current_blk, *(current_bid++));
        return *this;
    }

    //! Deallocates internal objects.
    ~buf_ostream()
    {
        assert(current_elem == 0);
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_MNG_BUF_OSTREAM_HEADER
