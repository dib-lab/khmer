/***************************************************************************
 *  tests/mng/test_prefetch_pool.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example mng/test_prefetch_pool.cpp

#include <iostream>
#include <stxxl/mng>
#include <stxxl/bits/mng/prefetch_pool.h>

#define BLOCK_SIZE (1024 * 512)

struct MyType
{
    int integer;
    char chars[5];
};

typedef stxxl::typed_block<BLOCK_SIZE, MyType> block_type;

// forced instantiation
template class stxxl::prefetch_pool<block_type>;

int main()
{
    stxxl::prefetch_pool<block_type> pool(2);
    pool.resize(10);
    pool.resize(5);

    block_type* blk = new block_type;
    (*blk)[0].integer = 42;
    block_type::bid_type bids[2];
    stxxl::block_manager::get_instance()->new_blocks(stxxl::single_disk(), bids, bids + 2);
    blk->write(bids[0])->wait();
    blk->write(bids[1])->wait();

    pool.hint(bids[0]);
    pool.read(blk, bids[0])->wait();
    pool.read(blk, bids[1])->wait();

    delete blk;
}
