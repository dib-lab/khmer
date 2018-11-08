/***************************************************************************
 *  tests/algo/test_sort_all_parameters.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//#define PLAY_WITH_OPT_PREF

#include <stxxl/mng>
#include <stxxl/sort>
#include <stxxl/scan>
#include <stxxl/vector>
#include <stxxl/random>

#include "test_sort_all_parameters.h"

#ifndef RECORD_SIZE
 #define RECORD_SIZE 4
#endif

#define MB (1024 * 1024)

template <typename T, typename alloc_strategy_type, unsigned block_size>
void test(stxxl::uint64 data_mem, unsigned memory_to_use)
{
    stxxl::uint64 records_to_sort = data_mem / sizeof(T);
    typedef stxxl::vector<T, 2, stxxl::lru_pager<8>, block_size, alloc_strategy_type> vector_type;
    vector_type v(records_to_sort);

    size_t ndisks = stxxl::config::get_instance()->disks_number();
    STXXL_MSG("Sorting " << records_to_sort << " records of size " << sizeof(T));
    STXXL_MSG("Total volume " << (records_to_sort * sizeof(T)) / MB << " MiB");
    STXXL_MSG("Using " << memory_to_use / MB << " MiB");
    STXXL_MSG("Using " << ndisks << " disks");
    STXXL_MSG("Using " << alloc_strategy_type::name() << " allocation strategy ");
    STXXL_MSG("Block size " << vector_type::block_type::raw_size / 1024 << " KiB");

    STXXL_MSG("Filling vector...");
    stxxl::generate(v.begin(), v.end(), stxxl::random_number32_r(), 32);

    STXXL_MSG("Sorting vector...");

    stxxl::stats_data before(*stxxl::stats::get_instance());

    stxxl::sort(v.begin(), v.end(), Cmp<T>(), memory_to_use);

    stxxl::stats_data after(*stxxl::stats::get_instance());

    STXXL_MSG("Checking order...");
    STXXL_CHECK(stxxl::is_sorted(v.begin(), v.end(), Cmp<T>()));

    STXXL_MSG("Sorting: " << (after - before));
    STXXL_MSG("Total:   " << *stxxl::stats::get_instance());
}

template <typename T, unsigned block_size>
void test_all_strategies(
    stxxl::uint64 data_mem,
    unsigned memory_to_use,
    int strategy)
{
    switch (strategy)
    {
    case 0:
        test<T, stxxl::striping, block_size>(data_mem, memory_to_use);
        break;
    case 1:
        test<T, stxxl::SR, block_size>(data_mem, memory_to_use);
        break;
    case 2:
        test<T, stxxl::FR, block_size>(data_mem, memory_to_use);
        break;
    case 3:
        test<T, stxxl::RC, block_size>(data_mem, memory_to_use);
        break;
    default:
        STXXL_ERRMSG("Unknown allocation strategy: " << strategy << ", aborting");
        abort();
    }
}

int main(int argc, char* argv[])
{
    if (argc < 6)
    {
        STXXL_ERRMSG("Usage: " << argv[0] <<
                     " <MiB to sort> <MiB to use> <alloc_strategy [0..3]> <blk_size [0..14]> <seed>");
        return -1;
    }

#if STXXL_PARALLEL_MULTIWAY_MERGE
    STXXL_MSG("STXXL_PARALLEL_MULTIWAY_MERGE");
#endif
    stxxl::uint64 data_mem = stxxl::atouint64(argv[1]) * MB;
    int sort_mem = atoi(argv[2]) * MB;
    int strategy = atoi(argv[3]);
    int block_size = atoi(argv[4]);
    stxxl::set_seed((unsigned)strtoul(argv[5], NULL, 10));
    STXXL_MSG("Seed " << stxxl::get_next_seed());
    stxxl::srandom_number32();

    typedef my_type<unsigned, RECORD_SIZE> my_default_type;

    switch (block_size)
    {
    case 0:
        test_all_strategies<my_default_type, 128* 1024>(data_mem, sort_mem, strategy);
        break;
    case 1:
        test_all_strategies<my_default_type, 256* 1024>(data_mem, sort_mem, strategy);
        break;
    case 2:
        test_all_strategies<my_default_type, 512* 1024>(data_mem, sort_mem, strategy);
        break;
    case 3:
        test_all_strategies<my_default_type, 1024* 1024>(data_mem, sort_mem, strategy);
        break;
    case 4:
        test_all_strategies<my_default_type, 2* 1024* 1024>(data_mem, sort_mem, strategy);
        break;
    case 5:
        test_all_strategies<my_default_type, 4* 1024* 1024>(data_mem, sort_mem, strategy);
        break;
    case 6:
        test_all_strategies<my_default_type, 8* 1024* 1024>(data_mem, sort_mem, strategy);
        break;
    case 7:
        test_all_strategies<my_default_type, 16* 1024* 1024>(data_mem, sort_mem, strategy);
        break;
    case 8:
        test_all_strategies<my_default_type, 640* 1024>(data_mem, sort_mem, strategy);
        break;
    case 9:
        test_all_strategies<my_default_type, 768* 1024>(data_mem, sort_mem, strategy);
        break;
    case 10:
        test_all_strategies<my_default_type, 896* 1024>(data_mem, sort_mem, strategy);
        break;
    case 11:
        test_all_strategies<my_type<unsigned, 12>, 2* MB>(data_mem, sort_mem, strategy);
        break;
    case 12:
        test_all_strategies<my_type<unsigned, 12>, 2* MB + 4096>(data_mem, sort_mem, strategy);
        break;
    case 13:
        test_all_strategies<my_type<unsigned, 20>, 2* MB + 4096>(data_mem, sort_mem, strategy);
        break;
    case 14:
        test_all_strategies<my_type<unsigned, 128>, 2* MB>(data_mem, sort_mem, strategy);
        break;
    default:
        STXXL_ERRMSG("Unknown block size: " << block_size << ", aborting");
        abort();
    }

    return 0;
}
