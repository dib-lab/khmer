/***************************************************************************
 *  tests/mng/test_block_alloc_strategy.cpp
 *
 *  instantiate all block allocation strategies
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/mng>
#include <stxxl/bits/mng/block_alloc_interleaved.h>
#include <stxxl/bits/common/error_handling.h>

template <typename strategy>
void test_strategy()
{
    STXXL_MSG(STXXL_PRETTY_FUNCTION_NAME);
    strategy s0;
    strategy s2(1, 3);
    stxxl::offset_allocator<strategy> o(1, s0);
    typedef typename stxxl::interleaved_alloc_traits<strategy>::strategy interleaved;
    interleaved itl(23, strategy(1, 3));
    for (int i = 0; i < 16; ++i)
        STXXL_MSG(s0(i) << " " << s2(i) << " " << o(i) << " " << itl(i));

    int firstdisk = 0;
    int ndisks = 4;
    int nruns = 10;
    int runsize = 15;
    std::cout << "firstdisk=" << firstdisk << "  ndisks=" << ndisks << "  nruns=" << nruns << "  runsize=" << runsize;
    interleaved itl2(nruns, strategy(firstdisk, firstdisk + ndisks));
    for (int i = 0; i < nruns * runsize; ++i) {
        if (i % nruns == 0)
            std::cout << std::endl;
        std::cout << itl2(i) << " ";
    }
    std::cout << std::endl;
}

int main()
{
    stxxl::config* cfg = stxxl::config::get_instance();

    // instantiate the allocation strategies
    STXXL_MSG("Number of disks: " << cfg->disks_number());
    for (unsigned i = 0; i < cfg->disks_number(); ++i)
        STXXL_MSG(cfg->disk_path(i) << ", " << cfg->disk_size(i) << ", " << cfg->disk_io_impl(i));
    test_strategy<stxxl::striping>();
    test_strategy<stxxl::FR>();
    test_strategy<stxxl::SR>();
    test_strategy<stxxl::RC>();
    STXXL_MSG("Regular disks: [" << cfg->regular_disk_range().first << "," << cfg->regular_disk_range().second << ")");
    if (cfg->regular_disk_range().first != cfg->regular_disk_range().second)
        test_strategy<stxxl::RC_disk>();
    STXXL_MSG("Flash devices: [" << cfg->flash_range().first << "," << cfg->flash_range().second << ")");
    if (cfg->flash_range().first != cfg->flash_range().second)
        test_strategy<stxxl::RC_flash>();
    test_strategy<stxxl::single_disk>();
}
