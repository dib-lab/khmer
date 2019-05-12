/***************************************************************************
 *  tests/mng/test_pool_pair.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example mng/test_pool_pair.cpp

#include <iostream>
#include <stxxl/mng>
#include <stxxl/bits/mng/prefetch_pool.h>
#include <stxxl/bits/mng/write_pool.h>

#define BLOCK_SIZE (1024 * 512)

struct MyType
{
    int integer;
    char chars[5];
};

typedef stxxl::typed_block<BLOCK_SIZE, MyType> block_type;

// forced instantiation
template class stxxl::typed_block<BLOCK_SIZE, MyType>;
template class stxxl::prefetch_pool<block_type>;
template class stxxl::write_pool<block_type>;

int main()
{
    stxxl::block_manager* bm = stxxl::block_manager::get_instance();
    STXXL_DEFAULT_ALLOC_STRATEGY alloc;

    {
        STXXL_MSG("Write-After-Write coherence test");
        stxxl::prefetch_pool<block_type> p_pool(2);
        stxxl::write_pool<block_type> w_pool(10);
        block_type* blk;
        block_type::bid_type bid;

        bm->new_block(alloc, bid);

        // write the block for the first time
        blk = w_pool.steal();
        (*blk)[0].integer = 42;
        w_pool.write(blk, bid);

        // read the block
        blk = w_pool.steal();
        p_pool.read(blk, bid)->wait();
        delete blk;

        // write the block for the second time
        blk = w_pool.steal();
        (*blk)[0].integer = 23;
        w_pool.write(blk, bid);

        // hint the block
        p_pool.hint(bid, w_pool); // flush w_pool

        // get the hinted block
        blk = w_pool.steal();
        p_pool.read(blk, bid)->wait();

        STXXL_CHECK2((*blk)[0].integer == 23,
                     "WRITE-AFTER-WRITE COHERENCE FAILURE");

        w_pool.add(blk);
        bm->delete_block(bid);
    }

    {
        STXXL_MSG("Write-After-Hint coherence test #1");
        stxxl::prefetch_pool<block_type> p_pool(1);
        stxxl::write_pool<block_type> w_pool(1);
        block_type* blk;
        block_type::bid_type bid;

        bm->new_block(alloc, bid);
        blk = w_pool.steal();
        (*blk)[0].integer = 42;
        w_pool.write(blk, bid);
        blk = w_pool.steal(); // flush w_pool

        // hint the block
        p_pool.hint(bid);

        // update the hinted block
        (*blk)[0].integer = 23;
        w_pool.write(blk, bid);
        p_pool.invalidate(bid);
        blk = w_pool.steal(); // flush w_pool

        // get the hinted block
        p_pool.read(blk, bid)->wait();

        STXXL_CHECK2((*blk)[0].integer == 23,
                     "WRITE-AFTER-HINT COHERENCE FAILURE");

        w_pool.add(blk);
        bm->delete_block(bid);
    }

    {
        STXXL_MSG("Write-After-Hint coherence test #2");
        stxxl::prefetch_pool<block_type> p_pool(1);
        stxxl::write_pool<block_type> w_pool(1);
        block_type* blk;
        block_type::bid_type bid;

        bm->new_block(alloc, bid);
        blk = w_pool.steal();
        (*blk)[0].integer = 42;
        w_pool.write(blk, bid);

        // hint the block
        p_pool.hint(bid, w_pool); // flush w_pool

        // update the hinted block
        blk = w_pool.steal();
        (*blk)[0].integer = 23;
        w_pool.write(blk, bid);
        p_pool.invalidate(bid);
        blk = w_pool.steal(); // flush w_pool

        // get the hinted block
        p_pool.read(blk, bid)->wait();

        STXXL_CHECK2((*blk)[0].integer == 23,
                     "WRITE-AFTER-HINT COHERENCE FAILURE");

        w_pool.add(blk);
        bm->delete_block(bid);
    }

    return 0;
}
// vim: et:ts=4:sw=4
