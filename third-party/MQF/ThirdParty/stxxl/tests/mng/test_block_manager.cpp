/***************************************************************************
 *  tests/mng/test_block_manager.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example mng/test_mng.cpp
//! This is an example of use of completion handlers, \c stxxl::block_manager, and
//! \c stxxl::typed_block

#include <iostream>
#include <stxxl/request>
#include <stxxl/mng>
#include <stxxl/bits/verbose.h>

#define BLOCK_SIZE (1024 * 512)

struct MyType
{
    int integer;
    //char chars[4];
    ~MyType() { }
};

struct my_handler
{
    void operator () (stxxl::request* req)
    {
        STXXL_MSG(req << " done, type=" << req->io_type());
    }
};

typedef stxxl::typed_block<BLOCK_SIZE, MyType> block_type;

int main()
{
    STXXL_MSG(sizeof(MyType) << " " << (BLOCK_SIZE % sizeof(MyType)));
    STXXL_MSG(sizeof(block_type) << " " << BLOCK_SIZE);
    const unsigned nblocks = 2;
    stxxl::BIDArray<BLOCK_SIZE> bids(nblocks);
    std::vector<int> disks(nblocks, 2);
    stxxl::request_ptr* reqs = new stxxl::request_ptr[nblocks];
    stxxl::block_manager* bm = stxxl::block_manager::get_instance();
    bm->new_blocks(stxxl::striping(), bids.begin(), bids.end());

    block_type* block = new block_type[2];
    STXXL_MSG(std::hex);
    STXXL_MSG("Allocated block address    : " << (stxxl::unsigned_type)(block));
    STXXL_MSG("Allocated block address + 1: " << (stxxl::unsigned_type)(block + 1));
    STXXL_MSG(std::dec);
    unsigned i = 0;
    for (i = 0; i < block_type::size; ++i)
    {
        block->elem[i].integer = i;
        //memcpy (block->elem[i].chars, "STXXL", 4);
    }
    for (i = 0; i < nblocks; ++i)
        reqs[i] = block->write(bids[i], my_handler());

    std::cout << "Waiting " << std::endl;
    stxxl::wait_all(reqs, nblocks);

    for (i = 0; i < nblocks; ++i)
    {
        reqs[i] = block->read(bids[i], my_handler());
        reqs[i]->wait();
        for (int j = 0; j < block_type::size; ++j)
        {
            STXXL_CHECK2(j == block->elem[j].integer,
                         "Error in block " << std::hex << i << " pos: " << j
                                           << " value read: " << block->elem[j].integer);
        }
    }

    bm->delete_blocks(bids.begin(), bids.end());

    delete[] reqs;
    delete[] block;

#if 0
    // variable-size blocks, not supported currently

    BIDArray<0> vbids(nblocks);
    for (i = 0; i < nblocks; i++)
        vbids[i].size = 1024 + i;

    bm->new_blocks(striping(), vbids.begin(), vbids.end());

    for (i = 0; i < nblocks; i++)
        STXXL_MSG("Allocated block: offset=" << vbids[i].offset << ", size=" << vbids[i].size);

    bm->delete_blocks(vbids.begin(), vbids.end());
#endif
}
