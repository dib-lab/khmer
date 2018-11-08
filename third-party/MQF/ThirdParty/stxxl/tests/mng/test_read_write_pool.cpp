/***************************************************************************
 *  tests/mng/test_read_write_pool.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example mng/test_read_write_pool.cpp

#include <iostream>
#include <stxxl/mng>
#include <stxxl/bits/mng/read_write_pool.h>

#define BLOCK_SIZE (1024 * 512)

struct MyType
{
    int integer;
    char chars[5];
};

typedef stxxl::typed_block<BLOCK_SIZE, MyType> block_type;

// forced instantiation
template class stxxl::typed_block<BLOCK_SIZE, MyType>;
template class stxxl::read_write_pool<block_type>;

int main()
{
    stxxl::block_manager* bm = stxxl::block_manager::get_instance();
    STXXL_DEFAULT_ALLOC_STRATEGY alloc;

    {
        STXXL_MSG("Write-After-Write coherence test");
        stxxl::read_write_pool<block_type> pool(2, 10);
        block_type* blk;
        block_type::bid_type bid;

        bm->new_block(alloc, bid);

        // write the block for the first time
        blk = pool.steal();
        (*blk)[0].integer = 42;
        pool.write(blk, bid);

        // read the block
        blk = pool.steal();
        pool.read(blk, bid)->wait();
        delete blk;

        // write the block for the second time
        blk = pool.steal();
        (*blk)[0].integer = 23;
        pool.write(blk, bid);

        // hint the block
        pool.hint(bid); // flush w_pool

        // get the hinted block
        blk = pool.steal();
        pool.read(blk, bid)->wait();

        STXXL_CHECK2((*blk)[0].integer == 23,
                     "WRITE-AFTER-WRITE COHERENCE FAILURE");

        pool.add(blk);
        bm->delete_block(bid);
    }

    {
        STXXL_MSG("Write-After-Hint coherence test #1");
        stxxl::read_write_pool<block_type> pool(1, 1);
        block_type* blk;
        block_type::bid_type bid;

        bm->new_block(alloc, bid);
        blk = pool.steal();
        (*blk)[0].integer = 42;
        pool.write(blk, bid);
        blk = pool.steal(); // flush w_pool

        // hint the block
        pool.hint(bid);

        // update the hinted block
        (*blk)[0].integer = 23;
        pool.write(blk, bid);
        blk = pool.steal(); // flush w_pool

        // get the hinted block
        pool.read(blk, bid)->wait();

        STXXL_CHECK2((*blk)[0].integer == 23,
                     "WRITE-AFTER-HINT COHERENCE FAILURE");

        pool.add(blk);
        bm->delete_block(bid);
    }

    {
        STXXL_MSG("Write-After-Hint coherence test #2");
        stxxl::read_write_pool<block_type> pool(1, 1);
        block_type* blk;
        block_type::bid_type bid;

        bm->new_block(alloc, bid);
        blk = pool.steal();
        (*blk)[0].integer = 42;
        pool.write(blk, bid);

        // hint the block
        pool.hint(bid);

        // update the hinted block
        blk = pool.steal();
        (*blk)[0].integer = 23;
        pool.write(blk, bid);
        blk = pool.steal(); // flush w_pool

        // get the hinted block
        pool.read(blk, bid)->wait();

        STXXL_CHECK2((*blk)[0].integer == 23,
                     "WRITE-AFTER-HINT COHERENCE FAILURE");

        pool.add(blk);
        bm->delete_block(bid);
    }

    return 0;
}
// vim: et:ts=4:sw=4
