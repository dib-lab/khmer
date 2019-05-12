/***************************************************************************
 *  tests/containers/hash_map/test_hash_map_block_cache.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Markus Westphal <marwes@users.sourceforge.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>

#include <stxxl.h>
#include <stxxl/bits/common/seed.h>
#include <stxxl/bits/containers/hash_map/block_cache.h>

bool test_block_cache()
{
    typedef std::pair<int, int> value_type;

    const unsigned subblock_raw_size = 1024 * 8; // 8KB subblocks
    const unsigned block_size = 128;             // 1MB blocks (=128 subblocks)

    const unsigned num_blocks = 64;              // number of blocks to use for this test
    const unsigned cache_size = 8;               // size of cache in blocks

    typedef stxxl::typed_block<subblock_raw_size, value_type> subblock_type;
    typedef stxxl::typed_block<block_size* sizeof(subblock_type), subblock_type> block_type;

    const unsigned subblock_size = subblock_type::size;          // size in values

    typedef block_type::bid_type bid_type;
    typedef std::vector<bid_type> bid_container_type;

    // prepare test: allocate blocks, fill them with values and write to disk
    bid_container_type bids(num_blocks);
    stxxl::block_manager* bm = stxxl::block_manager::get_instance();
    bm->new_blocks(stxxl::striping(), bids.begin(), bids.end());

    block_type* block = new block_type;
    for (unsigned i_block = 0; i_block < num_blocks; i_block++) {
        for (unsigned i_subblock = 0; i_subblock < block_size; i_subblock++) {
            for (unsigned i_value = 0; i_value < subblock_size; i_value++) {
                int value = i_value + i_subblock * subblock_size + i_block * block_size;
                (*block)[i_subblock][i_value] = value_type(value, value);
            }
        }
        stxxl::request_ptr req = block->write(bids[i_block]);
        req->wait();
    }

    stxxl::random_number32 rand32;

    // create block_cache
    typedef stxxl::hash_map::block_cache<block_type> cache_type;
    cache_type cache(cache_size);

    // load random subblocks and check for values
    int n_runs = cache_size * 10;
    for (int i_run = 0; i_run < n_runs; i_run++) {
        int i_block = rand32() % num_blocks;
        int i_subblock = rand32() % block_size;

        subblock_type* subblock = cache.get_subblock(bids[i_block], i_subblock);

        int expected = i_block * block_size + i_subblock * subblock_size + 1;
        STXXL_CHECK((*subblock)[1].first == expected);
    }

    // do the same again but this time with prefetching
    for (int i_run = 0; i_run < n_runs; i_run++) {
        int i_block = rand32() % num_blocks;
        int i_subblock = rand32() % block_size;

        cache.prefetch_block(bids[i_block]);
        subblock_type* subblock = cache.get_subblock(bids[i_block], i_subblock);
        int expected = i_block * block_size + i_subblock * subblock_size + 1;
        STXXL_CHECK((*subblock)[1].first == expected);
    }

    // load and modify some subblocks; flush cache and check values
    unsigned myseed = stxxl::get_next_seed();
    stxxl::set_seed(myseed);
    for (int i_run = 0; i_run < n_runs; i_run++) {
        int i_block = rand32() % num_blocks;
        int i_subblock = rand32() % block_size;

        subblock_type* subblock = cache.get_subblock(bids[i_block], i_subblock);

        STXXL_CHECK(cache.make_dirty(bids[i_block]));
        (*subblock)[1].first = (*subblock)[1].second + 42;
    }
    stxxl::set_seed(myseed);
    for (int i_run = 0; i_run < n_runs; i_run++) {
        int i_block = rand32() % num_blocks;
        int i_subblock = rand32() % block_size;
        subblock_type* subblock = cache.get_subblock(bids[i_block], i_subblock);

        int expected = i_block * block_size + i_subblock * subblock_size + 1;
        STXXL_CHECK((*subblock)[1].first == expected + 42);
    }

    // test retaining
    cache.clear();

    // not yet cached
    STXXL_CHECK(cache.retain_block(bids[0]) == false);
    cache.prefetch_block(bids[0]);

    // cached, should be retained
    STXXL_CHECK(cache.retain_block(bids[0]) == true);
    // release again
    STXXL_CHECK(cache.release_block(bids[0]) == true);
    // retrain-count should be 0, release fails
    STXXL_CHECK(cache.release_block(bids[0]) == false);

    // cache new block
    subblock_type* kicked_subblock = cache.get_subblock(bids[1], 0);
    // load other blocks, so that kicked_subblock, well, gets kicked
    for (unsigned i = 0; i < cache_size + 5; i++) {
        cache.prefetch_block(bids[i + 3]);
    }
    // load kicked subblock again, should be at a different location
    STXXL_CHECK(cache.get_subblock(bids[1], 0) != kicked_subblock);

    subblock_type* retained_subblock = cache.get_subblock(bids[1], 0);
    // now retain subblock
    STXXL_CHECK(cache.retain_block(bids[1]) == true);
    for (unsigned i = 0; i < cache_size + 5; i++) {
        cache.prefetch_block(bids[i + 3]);
    }
    // retained_subblock should not have been kicked
    STXXL_CHECK(cache.get_subblock(bids[1], 0) == retained_subblock);
    cache.clear();

    // test swapping
    subblock_type* a_subblock = cache.get_subblock(bids[6], 1);
    cache_type cache2(cache_size / 2);
    std::swap(cache, cache2);
    STXXL_CHECK(cache.size() == cache_size / 2);
    STXXL_CHECK(cache2.size() == cache_size);
    STXXL_CHECK(cache2.get_subblock(bids[6], 1) == a_subblock);

    STXXL_MSG("Passed Block-Cache Test");

    return true;
}

int main()
{
    test_block_cache();
    return 0;
}
