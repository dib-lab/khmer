/***************************************************************************
 *  tests/algo/test_asch.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <cstdlib>
#include <stxxl/bits/algo/async_schedule.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/random>

// Test async schedule algorithm

int main(int argc, char* argv[])
{
    if (argc < 5)
    {
        STXXL_ERRMSG("Usage: " << argv[0] << " D L m seed");
        return -1;
    }
    const int D = atoi(argv[1]);
    const int L = atoi(argv[2]);
    const stxxl::int_type m = atoi(argv[3]);
    stxxl::ran32State = atoi(argv[4]);
    stxxl::int_type* disks = new stxxl::int_type[L];
    stxxl::int_type* prefetch_order = new stxxl::int_type[L];
    int* count = new int[D];

    for (int i = 0; i < D; i++)
        count[i] = 0;

    stxxl::random_number32 rnd;
    for (int i = 0; i < L; i++)
    {
        disks[i] = rnd(D);
        count[disks[i]]++;
    }

    for (int i = 0; i < D; i++)
        std::cout << "Disk " << i << " has " << count[i] << " blocks" << std::endl;

    stxxl::compute_prefetch_schedule(disks, disks + L, prefetch_order, m, D);

    STXXL_MSG("Prefetch order:");
    for (int i = 0; i < L; ++i) {
        STXXL_MSG("request " << prefetch_order[i] << "  on disk " << disks[prefetch_order[i]]);
    }
    STXXL_MSG("Request order:");
    for (int i = 0; i < L; ++i) {
        int j;
        for (j = 0; prefetch_order[j] != i; ++j) ;
        STXXL_MSG("request " << i << "  on disk " << disks[i] << "  scheduled as " << j);
    }

    delete[] count;
    delete[] disks;
    delete[] prefetch_order;
}

// vim: et:ts=4:sw=4
