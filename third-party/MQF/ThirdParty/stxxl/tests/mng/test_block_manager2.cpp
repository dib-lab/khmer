/***************************************************************************
 *  tests/mng/test_block_manager2.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>
#include <stxxl/mng>

#define BLOCK_SIZE (1024 * 1024 * 32)

typedef stxxl::typed_block<BLOCK_SIZE, int> block_type;
template class stxxl::typed_block<BLOCK_SIZE, int>; // forced instantiation

int main()
{
    stxxl::int64 totalsize = 0;
    stxxl::config* config = stxxl::config::get_instance();

    for (size_t i = 0; i < config->disks_number(); ++i)
        totalsize += config->disk_size(i);

    stxxl::unsigned_type totalblocks =
        (stxxl::unsigned_type)(totalsize / block_type::raw_size);

    STXXL_MSG("external memory: " << totalsize << " bytes  ==  " << totalblocks << " blocks");

    stxxl::BIDArray<BLOCK_SIZE> b5a(totalblocks / 5);
    stxxl::BIDArray<BLOCK_SIZE> b5b(totalblocks / 5);
    stxxl::BIDArray<BLOCK_SIZE> b5c(totalblocks / 5);
    stxxl::BIDArray<BLOCK_SIZE> b5d(totalblocks / 5);
    stxxl::BIDArray<BLOCK_SIZE> b2(totalblocks / 2);

    stxxl::block_manager* bm = stxxl::block_manager::get_instance();

    STXXL_MSG("get 4 x " << totalblocks / 5);
    bm->new_blocks(stxxl::striping(), b5a.begin(), b5a.end());
    bm->new_blocks(stxxl::striping(), b5b.begin(), b5b.end());
    bm->new_blocks(stxxl::striping(), b5c.begin(), b5c.end());
    bm->new_blocks(stxxl::striping(), b5d.begin(), b5d.end());

    STXXL_MSG("free 2 x " << totalblocks / 5);
    bm->delete_blocks(b5a.begin(), b5a.end());
    bm->delete_blocks(b5c.begin(), b5c.end());

    // the external memory should now be fragmented enough,
    // s.t. the following request needs to be split into smaller ones
    STXXL_MSG("get 1 x " << totalblocks / 2);
    bm->new_blocks(stxxl::striping(), b2.begin(), b2.end());

    bm->delete_blocks(b5b.begin(), b5b.end());
    bm->delete_blocks(b5d.begin(), b5d.end());

    bm->delete_blocks(b2.begin(), b2.end());
}
