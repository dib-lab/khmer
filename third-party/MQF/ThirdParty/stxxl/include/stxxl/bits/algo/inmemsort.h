/***************************************************************************
 *  include/stxxl/bits/algo/inmemsort.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_INMEMSORT_HEADER
#define STXXL_ALGO_INMEMSORT_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/common/simple_vector.h>
#include <stxxl/bits/io/request_operations.h>
#include <stxxl/bits/algo/adaptor.h>
#include <stxxl/bits/mng/adaptor.h>
#include <stxxl/bits/parallel.h>

#include <algorithm>

STXXL_BEGIN_NAMESPACE

template <typename ExtIterator, typename StrictWeakOrdering>
void stl_in_memory_sort(ExtIterator first, ExtIterator last, StrictWeakOrdering cmp)
{
    typedef typename ExtIterator::block_type block_type;

    STXXL_VERBOSE("stl_in_memory_sort, range: " << (last - first));
    first.flush();
    unsigned_type nblocks = last.bid() - first.bid() + (last.block_offset() ? 1 : 0);
    simple_vector<block_type> blocks(nblocks);
    simple_vector<request_ptr> reqs(nblocks);
    unsigned_type i;

    for (i = 0; i < nblocks; ++i)
        reqs[i] = blocks[i].read(*(first.bid() + i));

    wait_all(reqs.begin(), nblocks);

    unsigned_type last_block_correction = last.block_offset() ? (block_type::size - last.block_offset()) : 0;
    check_sort_settings();
    potentially_parallel::
    sort(make_element_iterator(blocks.begin(), first.block_offset()),
         make_element_iterator(blocks.begin(), nblocks * block_type::size - last_block_correction),
         cmp);

    for (i = 0; i < nblocks; ++i)
        reqs[i] = blocks[i].write(*(first.bid() + i));

    wait_all(reqs.begin(), nblocks);
}

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_INMEMSORT_HEADER
