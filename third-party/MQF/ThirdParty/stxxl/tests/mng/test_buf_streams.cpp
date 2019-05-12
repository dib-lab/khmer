/***************************************************************************
 *  tests/mng/test_buf_streams.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example mng/test_buf_streams.cpp
//! This is an example of use of \c stxxl::buf_istream and \c stxxl::buf_ostream

#include <iostream>
#include <stxxl/mng>
#include <stxxl/bits/mng/buf_ostream.h>
#include <stxxl/bits/mng/buf_istream.h>
#include <stxxl/bits/mng/buf_istream_reverse.h>

#define BLOCK_SIZE (1024 * 512)

typedef stxxl::typed_block<BLOCK_SIZE, unsigned> block_type;
typedef stxxl::BIDArray<BLOCK_SIZE>::iterator bid_iterator_type;

typedef stxxl::buf_ostream<block_type, bid_iterator_type> buf_ostream_type;
typedef stxxl::buf_istream<block_type, bid_iterator_type> buf_istream_type;
typedef stxxl::buf_istream_reverse<block_type, bid_iterator_type> buf_istream_reverse_type;

// forced instantiations
template class stxxl::buf_ostream<block_type, stxxl::BIDArray<BLOCK_SIZE>::iterator>;
template class stxxl::buf_istream<block_type, stxxl::BIDArray<BLOCK_SIZE>::iterator>;
template class stxxl::buf_istream_reverse<block_type, stxxl::BIDArray<BLOCK_SIZE>::iterator>;

int main()
{
    const unsigned nblocks = 128;
    const unsigned nelements = nblocks * block_type::size;
    stxxl::BIDArray<BLOCK_SIZE> bids(nblocks);

    stxxl::block_manager* bm = stxxl::block_manager::get_instance();
    bm->new_blocks(stxxl::striping(), bids.begin(), bids.end());
    {
        buf_ostream_type out(bids.begin(), 2);
        for (unsigned i = 0; i < nelements; i++)
            out << i;
    }
    {
        buf_istream_type in(bids.begin(), bids.end(), 2);
        for (unsigned i = 0; i < nelements; i++)
        {
            unsigned prevalue = *in;
            unsigned value;
            in >> value;

            STXXL_CHECK2(value == i,
                         "Error at position " << std::hex << i << " (" << value << ") block " << (i / block_type::size));
            STXXL_CHECK(prevalue == value);
        }
    }
    {
        buf_istream_reverse_type in(bids.begin(), bids.end(), 2);
        for (unsigned i = 0; i < nelements; i++)
        {
            unsigned prevalue = *in;
            unsigned value;
            in >> value;

            STXXL_CHECK2(value == nelements - i - 1,
                         "Error at position " << std::hex << i << " (" << value << ") block " << (i / block_type::size));
            STXXL_CHECK(prevalue == value);
        }
    }
    bm->delete_blocks(bids.begin(), bids.end());

    return 0;
}
